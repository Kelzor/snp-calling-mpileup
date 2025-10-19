#!/usr/bin/env python
# This line allows the script to be executed directly on UNIX-like systems

import pysam                   # For reading and processing VCF files efficiently
import argparse                # For handling command-line arguments
from Bio import SeqIO          # For writing FASTA output
from Bio.SeqRecord import SeqRecord  # Represents a sequence entry for FASTA output
from Bio.Seq import Seq        # Represents DNA or RNA sequences

# -------------------------------------------------------------------------
# MAIN FUNCTION
# -------------------------------------------------------------------------

def vcf_to_alignment(vcf_file, fasta_output, min_af, max_af, depth_threshold,
                     deletion_threshold=None, keep_monomorphic=False, no_debug=False):
    """
    Convert a multi-sample VCF file into a SNP alignment in FASTA format.

    Each VCF site corresponds to one column in the alignment (if it passes filters).
    Monomorphic sites are excluded by default (can be kept with --keep-monomorphic).
    """

    # Open the VCF file using pysam (handles both .vcf and .vcf.gz)
    with pysam.VariantFile(vcf_file, 'r') as vcf:

        # Extract list of sample IDs from the VCF header
        samples = list(vcf.header.samples)

        # Prepare a dictionary mapping each sample → list of base calls
        # Each list will later be concatenated into a full sequence
        sample_sequences = {sample: [] for sample in samples}

        # Counter for sites retained in the final alignment
        total_snps = 0

        # Prepare debug or summary output files depending on the flag
        debug_file = None
        summary_file = None
        if not no_debug:
            debug_file = open(fasta_output + ".debug.txt", 'w')
        else:
            summary_file = open(fasta_output + ".stats.txt", 'w')

        # -----------------------------------------------------------------
        # Iterate through each variant site in the VCF
        # -----------------------------------------------------------------
        for idx, record in enumerate(vcf):

            # Get reference and alternate alleles (strings)
            ref_allele = record.ref
            alt_alleles = record.alts

            # Skip multi-base variants (INDELs, MNPs) – only keep single-nucleotide polymorphisms
            if len(ref_allele) > 1 or any(len(alt) > 1 for alt in alt_alleles):
                continue

            # Store per-sample information for this site
            site_calls = {}
            non_missing = 0  # Counter for samples with a valid (non-"N") base

            # -------------------------------------------------------------
            # Process each sample at this site
            # -------------------------------------------------------------
            for sample in samples:

                if sample in record.samples:
                    # Extract genotype (GT), allele frequency (AF), and read depth (DP)
                    genotype = record.samples[sample]['GT']  # e.g., (0, 1) or (1, 1)
                    af = record.samples[sample].get('AF')    # float, allele frequency of ALT in this sample (if present)
                    dp = record.samples[sample].get('DP')    # int, read depth

                    # -----------------------------------------------------
                    # Check that genotype is valid and diploid (two alleles)
                    # -----------------------------------------------------
                    if isinstance(genotype, tuple) and len(genotype) == 2:
                        allele1, allele2 = genotype

                        # Only proceed if read depth is high enough
                        if dp is not None and dp >= depth_threshold:

                            # Case 1: Homozygous reference (0/0)
                            # → Assign reference base directly
                            if allele1 == 0 and allele2 == 0:
                                base = ref_allele

                            # Case 2: At least one ALT allele present (e.g., 0/1, 1/1)
                            else:
                                # -------------------------------------------------
                                # The allele frequency (AF) helps interpret the call:
                                # - AF is the proportion of reads supporting the ALT allele.
                                # - AF ~ 0.0 → nearly all reads are REF
                                # - AF ~ 1.0 → nearly all reads are ALT
                                # - AF around 0.5 → heterozygous or uncertain
                                #
                                # We use user-provided min_af / max_af thresholds to decide:
                                #   - AF >= max_af → trust ALT base
                                #   - AF <= min_af → trust REF base
                                #   - Otherwise → ambiguous (too intermediate), mark as 'N'
                                # -------------------------------------------------
                                if af is not None:
                                    allele_frequency_alt = af

                                    # High confidence ALT (e.g., AF ≥ 0.8)
                                    if allele_frequency_alt >= max_af:
                                        # allele1 - 1 converts genotype index (1-based for ALT) to 0-based Python index
                                        base = alt_alleles[allele1 - 1]

                                    # High confidence REF (e.g., AF ≤ 0.2)
                                    elif allele_frequency_alt <= min_af:
                                        base = ref_allele

                                    # Uncertain (intermediate frequency, likely heterozygous)
                                    else:
                                        base = 'N'
                                else:
                                    # If AF field is missing, mark as ambiguous
                                    base = 'N'
                        else:
                            # Not enough coverage → missing
                            base = 'N'

                    # Missing genotype field (e.g., "./.") → no call
                    elif genotype == './.':
                        base = 'N'

                    # Any unexpected case → mark missing
                    else:
                        base = 'N'

                else:
                    # Sample not found in record → missing
                    base = 'N'

                # Count valid bases (not "N")
                if base != 'N':
                    non_missing += 1

                # Store all details for this site and sample
                site_calls[sample] = (base, genotype, af, dp)

            # -------------------------------------------------------------
            # Deletion filter: require a minimum proportion of non-missing data
            # -------------------------------------------------------------
            if deletion_threshold is not None:
                if non_missing / len(samples) < deletion_threshold:
                    continue  # Skip site if too sparse

            # -------------------------------------------------------------
            # Skip monomorphic sites unless user allows keeping them
            # -------------------------------------------------------------
            if not keep_monomorphic:
                # Collect only called (non-"N") bases
                called_bases = [base for base, _, _, _ in site_calls.values() if base != 'N']

                # If all samples share the same base → not variable → skip
                if len(set(called_bases)) <= 1:
                    continue

            # -------------------------------------------------------------
            # Site passes all filters — keep it
            # -------------------------------------------------------------
            total_snps += 1  # Increment total site counter

            # Append this site’s base to each sample’s running sequence
            for sample in samples:
                base, genotype, af, dp = site_calls[sample]
                sample_sequences[sample].append(base)

                # Optionally log debugging information for every sample
                if not no_debug:
                    debug_file.write(
                        f"POS: {record.pos}, Sample: {sample}, Added to Alignment: {base}, "
                        f"Genotype: {genotype}, Allele Frequency: {af}, Depth: {dp}\n"
                    )

        # -----------------------------------------------------------------
        # After all sites processed: build final FASTA records
        # -----------------------------------------------------------------
        records = [
            SeqRecord(
                Seq(''.join(sample_sequences[sample])),  # Join list of bases → string sequence
                id=sample,                               # FASTA header ID
                description=""                           # Empty description line
            )
            for sample in samples
        ]

        # Write all sequences to FASTA output file
        with open(fasta_output, 'w') as output_file:
            SeqIO.write(records, output_file, "fasta")

        # Write summary or debug stats
        if not no_debug:
            debug_file.write(f"Total SNPs in final alignment: {total_snps}\n")
            debug_file.close()
        else:
            summary_file.write(f"Total SNPs in final alignment: {total_snps}\n")
            summary_file.write(f"Number of samples: {len(samples)}\n")
            summary_file.close()
            print(f"Total SNPs in final alignment: {total_snps}")
            print(f"Number of samples: {len(samples)}")

# -------------------------------------------------------------------------
# COMMAND-LINE INTERFACE
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Define and parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert VCF to SNP alignment FASTA.")

    # Input and output paths
    parser.add_argument("-i", "--vcf", required=True,
                        help="Input VCF file (can be .vcf or .vcf.gz)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA alignment file")

    # Allele frequency thresholds
    parser.add_argument("--min_af", type=float, default=0.2,
                        help="Minimum allele frequency threshold below which REF is assumed")
    parser.add_argument("--max_af", type=float, default=0.8,
                        help="Maximum allele frequency threshold above which ALT is assumed")

    # Depth and missing-data filters
    parser.add_argument("--depth", type=int, default=10,
                        help="Minimum read depth required for a valid genotype call")
    parser.add_argument("--deletion", type=float, default=None,
                        help="Drop sites with fewer than this proportion of non-missing calls")

    # Option to retain monomorphic sites
    parser.add_argument("--keep-monomorphic", action='store_true',
                        help="Include monomorphic sites (default: skip them)")

    # Option to suppress detailed debugging logs
    parser.add_argument("--no-debug", action='store_true',
                        help="Do not write per-site debug file; only summary stats")

    # Parse all arguments into a namespace object
    args = parser.parse_args()

    # Call main function with parsed parameters
    vcf_to_alignment(
        args.vcf,
        args.output,
        args.min_af,
        args.max_af,
        args.depth,
        args.deletion,
        args.keep_monomorphic,
        args.no_debug
    )