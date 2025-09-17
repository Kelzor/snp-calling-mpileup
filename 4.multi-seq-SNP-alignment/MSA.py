#!/usr/bin/env python

import pysam
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def vcf_to_alignment(vcf_file, fasta_output, min_af, max_af, depth_threshold, deletion_threshold=None, no_debug=False):
    """
    Convert a multi-sample VCF file into a SNP alignment FASTA.
    Optional deletion filtering and debug output.
    """

    sample_sequences = {}
    debug_entries = []

    # Open VCF
    with pysam.VariantFile(vcf_file, 'r') as vcf:
        samples = list(vcf.header.samples)
        sample_sequences = {sample: [] for sample in samples}
        site_positions = []

        for record in vcf:
            ref_allele = record.ref
            alt_alleles = record.alts

            # Skip INDELs
            if len(ref_allele) > 1 or any(len(alt) > 1 for alt in alt_alleles):
                continue

            site_calls = {}

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

            # Store debug entries
            for sample in samples:
                base, genotype, af, dp = site_calls[sample]
                debug_entries.append(
                    f"POS: {record.pos}, Sample: {sample}, Added to Alignment: {base}, "
                    f"Genotype: {genotype}, Allele Frequency: {af}, Depth: {dp}"
                )

            site_positions.append(record.pos)

    # Deletion filtering
    if deletion_threshold is not None:
        num_samples = len(samples)
        keep_sites = []
        num_sites = len(next(iter(sample_sequences.values())))
        for i in range(num_sites):
            non_missing = sum(seq[i] != 'N' for seq in sample_sequences.values())
            if non_missing / num_samples >= deletion_threshold:
                keep_sites.append(i)

        # Filter sequences
        for sample in samples:
            sample_sequences[sample] = [sample_sequences[sample][i] for i in keep_sites]

        # Filter debug entries
        kept_positions = {site_positions[i] for i in keep_sites}
        debug_entries = [entry for entry in debug_entries if f"POS: {entry.split('POS: ')[1].split(',')[0]}" in map(str, kept_positions)]

    # Build SeqRecords
    records = [SeqRecord(Seq(''.join(sample_sequences[sample])),
                         id=sample, description="") for sample in samples]

    # Write FASTA
    with open(fasta_output, 'w') as output_file:
        SeqIO.write(records, output_file, "fasta")

    # Total SNPs in final alignment
    total_snps = len(records[0].seq) if records else 0

    if not no_debug:
        # Write full debug file
        debug_file = fasta_output + ".debug.txt"
        with open(debug_file, 'w') as df:
            df.write(f"Total SNPs in final alignment: {total_snps}\n")
            for line in debug_entries:
                df.write(line + "\n")
    else:
        # Write a small summary file
        summary_file = fasta_output + ".stats.txt"
        with open(summary_file, 'w') as sf:
            sf.write(f"Total SNPs in final alignment: {total_snps}\n")
            sf.write(f"Number of samples: {len(samples)}\n")

        # Also print to console
        print(f"Total SNPs in final alignment: {total_snps}")
        print(f"Number of samples: {len(samples)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert VCF to SNP alignment FASTA.")
    parser.add_argument("-i", "--vcf", required=True, help="Input VCF file (can be .vcf or .vcf.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA alignment file")
    parser.add_argument("--min_af", type=float, default=0.2, help="Minimum allele frequency threshold")
    parser.add_argument("--max_af", type=float, default=0.8, help="Maximum allele frequency threshold")
    parser.add_argument("--depth", type=int, default=10, help="Minimum depth threshold")
    parser.add_argument("--deletion", type=float, default=None,
                        help="Optional filter: drop sites with fewer than this proportion of non-missing calls")
    parser.add_argument("--no-debug", action='store_true', help="Do not write full debug file; write only summary stats and print to console")

    args = parser.parse_args()

    vcf_to_alignment(args.vcf, args.output, args.min_af, args.max_af, args.depth, args.deletion, args.no_debug)
