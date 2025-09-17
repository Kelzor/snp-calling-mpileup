import pysam
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os


def vcf_to_alignment(vcf_file, fasta_output, min_af, max_af, depth_threshold, deletion_threshold=None):
    """
    Convert a multi-sample VCF file into a SNP alignment in FASTA format,
    applying allele frequency and depth thresholds, and skipping INDELs.
    Optionally filter alignment columns by deletion_threshold (proportion of Ns).
    """

    records = []        # SeqRecord objects for FASTA output
    debug_entries = []  # Debug info per site/sample

    # Open the VCF
    with pysam.VariantFile(vcf_file, 'r') as vcf:
        samples = list(vcf.header.samples)

        # Dictionary of sequences (lists of bases) per sample
        sample_sequences = {sample: [] for sample in samples}
        site_positions = []  # Track positions for filtering

        for record in vcf:
            ref_allele = record.ref
            alt_alleles = record.alts

            # Skip INDELs
            if len(ref_allele) > 1 or any(len(alt) > 1 for alt in alt_alleles):
                continue

            site_bases = []

            for sample in samples:
                if sample in record.samples:
                    genotype = record.samples[sample]['GT']
                    af = record.samples[sample].get('AF')
                    dp = record.samples[sample].get('DP')

                    if isinstance(genotype, tuple) and len(genotype) == 2:
                        allele1, allele2 = genotype

                        if dp is not None and dp >= depth_threshold:
                            if allele1 == 0 and allele2 == 0:
                                base = ref_allele
                            else:
                                if af is not None:
                                    allele_frequency_alt = af
                                    if allele_frequency_alt >= max_af:
                                        if allele1 > 0:
                                            base = alt_alleles[allele1 - 1]
                                        else:
                                            base = alt_alleles[allele2 - 1]
                                    elif allele_frequency_alt <= min_af:
                                        base = ref_allele
                                    else:
                                        base = 'N'
                                else:
                                    base = 'N'
                        else:
                            base = 'N'
                    elif genotype == './.':
                        base = 'N'
                    else:
                        base = 'N'
                else:
                    base = 'N'

                sample_sequences[sample].append(base)
                site_bases.append(base)

                debug_entries.append({
                    "pos": record.pos,
                    "sample": sample,
                    "base": base,
                    "genotype": genotype,
                    "af": af,
                    "dp": dp
                })

            site_positions.append(record.pos)

    # Apply deletion filtering if requested
    if deletion_threshold is not None:
        num_samples = len(samples)
        keep_sites = []

        for i in range(len(site_positions)):
            non_missing = sum(1 for sample in samples if sample_sequences[sample][i] != 'N')
            prop_non_missing = non_missing / num_samples
            if prop_non_missing >= deletion_threshold:
                keep_sites.append(i)

        # Filter sequences
        for sample in samples:
            sample_sequences[sample] = [sample_sequences[sample][i] for i in keep_sites]

        # Filter debug entries
        kept_positions = {site_positions[i] for i in keep_sites}
        debug_entries = [entry for entry in debug_entries if entry["pos"] in kept_positions]

    # Build SeqRecord objects
    for sample in samples:
        sequence = ''.join(sample_sequences[sample])
        record = SeqRecord(Seq(sequence), id=sample, description="")
        records.append(record)

    # Write FASTA alignment
    with open(fasta_output, 'w') as output_file:
        SeqIO.write(records, output_file, "fasta")

    # Name debug file based on output FASTA file
    debug_file_path = fasta_output + ".debug.txt"

    # Write debug log
    with open(debug_file_path, 'w') as debug_file:
        for entry in debug_entries:
            debug_file.write(
                f"POS: {entry['pos']}, Sample: {entry['sample']}, Added: {entry['base']}, "
                f"Genotype: {entry['genotype']}, AF: {entry['af']}, DP: {entry['dp']}\n"
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert VCF to SNP alignment FASTA.")
    parser.add_argument("-i", "--vcf", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("--min_af", type=float, default=0.2, help="Minimum allele frequency threshold (default: 0.2)")
    parser.add_argument("--max_af", type=float, default=0.8, help="Maximum allele frequency threshold (default: 0.8)")
    parser.add_argument("--depth", type=int, default=10, help="Depth threshold (default: 10)")
    parser.add_argument("--deletion", type=float, default=None,
                        help="Optional: filter sites where fewer than this proportion of samples have data (0.1–1.0)")

    args = parser.parse_args()

    vcf_to_alignment(args.vcf, args.output, args.min_af, args.max_af, args.depth, args.deletion)
