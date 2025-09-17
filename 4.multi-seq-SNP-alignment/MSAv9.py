import pysam
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def vcf_to_alignment(vcf_file, fasta_output, min_af, max_af, depth_threshold, deletion_threshold=None):
    records = []

    # Open the VCF file
    with pysam.VariantFile(vcf_file, 'r') as vcf:
        samples = list(vcf.header.samples)
        sample_sequences = {sample: [] for sample in samples}
        debug_entries = []  # store debug info in memory

        for record in vcf:
            ref_allele = record.ref
            alt_alleles = record.alts

            # Skip indels
            if len(ref_allele) > 1 or any(len(alt) > 1 for alt in alt_alleles):
                continue

            site_calls = {}  # store calls for each sample at this site

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
                                        base = alt_alleles[allele1 - 1]
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
                site_calls[sample] = (base, genotype, af, dp)

            # save debug entries for this site
            for sample in samples:
                base, genotype, af, dp = site_calls[sample]
                debug_entries.append(
                    f"POS: {record.pos}, Sample: {sample}, Added to Alignment: {base}, "
                    f"Genotype: {genotype}, Allele Frequency: {af}, Depth: {dp}"
                )

    # Apply deletion filtering if requested
    if deletion_threshold is not None:
        keep_sites = []
        num_sites = len(next(iter(sample_sequences.values())))
        num_samples = len(samples)
        for i in range(num_sites):
            non_missing = sum(seq[i] != 'N' for seq in sample_sequences.values())
            if non_missing / num_samples >= deletion_threshold:
                keep_sites.append(i)

        # filter sequences
        for sample in samples:
            sample_sequences[sample] = ''.join(sample_sequences[sample][i] for i in keep_sites)

        # filter debug entries
        filtered_debug = []
        site_index = 0
        for entry in debug_entries:
            if f"POS:" in entry:
                if site_index in keep_sites:
                    filtered_debug.append(entry)
                if "Sample:" in entry:
                    pass
            site_index += 0  # site_index already tracked above
        debug_entries = [entry for idx, entry in enumerate(debug_entries)
                         if (idx // len(samples)) in keep_sites]

    # Create SeqRecord objects
    records = [SeqRecord(Seq(''.join(sample_sequences[sample])),
                         id=sample, description="") for sample in samples]

    # Write FASTA
    with open(fasta_output, 'w') as output_file:
        SeqIO.write(records, output_file, "fasta")

    # Write debug file (with total SNP count at top)
    total_snps = len(records[0].seq) if records else 0
    debug_file = fasta_output + ".debug.txt"
    with open(debug_file, 'w') as df:
        df.write(f"Total SNPs in final alignment: {total_snps}\n")
        for line in debug_entries:
            df.write(line + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert VCF to SNP alignment FASTA.")
    parser.add_argument("-i", "--vcf", required=True, help="Input VCF file (can be .vcf or .vcf.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA alignment file")
    parser.add_argument("--min_af", type=float, default=0.2, help="Minimum allele frequency threshold")
    parser.add_argument("--max_af", type=float, default=0.8, help="Maximum allele frequency threshold")
    parser.add_argument("--depth", type=int, default=10, help="Minimum depth threshold")
    parser.add_argument("--deletion", type=float, default=None,
                        help="Optional filter: drop sites with fewer than this proportion of non-missing calls")
    args = parser.parse_args()

    vcf_to_alignment(args.vcf, args.output, args.min_af, args.max_af, args.depth, args.deletion)
