import argparse
import pysam

def parse_args():
    # Parse command line arguments to configure filtering options
    parser = argparse.ArgumentParser(description='Filter VCF file based on DP, AF, GT, SNP type, and SNP depth')
    parser.add_argument('vcf_file', help='input VCF file')  # Input VCF file
    parser.add_argument('-s', '--sample', help='comma-separated list of sample names to filter')  # Optional list of samples to filter
    parser.add_argument('-d', '--dp', type=int, required=True, help='minimum read depth (DP) to retain a variant')  # Minimum read depth for a sample to be retained
    parser.add_argument('-n', '--snpdp', type=int, required=True, help='minimum SNP depth (alternate allele depth in AD) to retain a variant')  # Minimum SNP depth to retain a variant
    parser.add_argument('-a', '--af', type=float, required=True, help='minimum allele frequency (AF) to retain a variant')  # Minimum allele frequency to retain a variant
    parser.add_argument('-b', '--bed', help='input BED file to specify regions to filter')  # Optional BED file for region-based filtering
    parser.add_argument('-o', '--output', required=True, help='output VCF file')  # Output VCF file
    parser.add_argument('-v', '--snp', action='store_true', help='only include SNPs (exclude indels)')  # Flag to include only SNPs, excluding indels
    args = parser.parse_args()

    # If sample names are provided, split them into a list; otherwise, use an empty list
    if args.sample:
        args.sample = args.sample.split(',')
    else:
        args.sample = []  # Set an empty list if no sample names are provided

    return args

def load_bed_regions(bed_file):
    # Load the regions from a BED file and return as a set of (chrom, position) tuples
    bed_regions = set()
    with open(bed_file, 'r') as bed:
        for line in bed:
            parts = line.strip().split('\t')  # Split the line by tabs
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])  # Extract chromosome, start, and end positions
            for position in range(start, end + 1):
                bed_regions.add((chrom, position))  # Add each position in the range to the set
    return bed_regions

def is_snp(record):
    # Check if the record is a SNP (not an indel)
    # SNPs do not have the key 'INDEL' in their INFO field
    return 'INDEL' not in record.info

def filter_record(record, args, bed_regions):
    # Filter out the record based on various conditions (DP, AF, GT, AD, SNP type, and BED regions)

    # Filter out indels if --snp is specified by the user
    if args.snp and not is_snp(record):
        return False  # If it's an indel, and only SNPs are requested, return False to exclude this record

    # Check if the variant's position is within any of the specified BED regions
    if (record.chrom, record.pos) in bed_regions:
        return False  # Exclude the variant if it falls within the BED regions

    # Loop through all samples in the VCF record
    for sample_name in record.samples.keys():
        # Skip the user-specified samples (those listed in args.sample)
        if sample_name in args.sample:
            continue  # Move on to the next sample if the current sample is specified for exclusion

        # Get the data for the current sample
        sample = record.samples[sample_name]

        # Check if the sample passes the DP filter (minimum read depth)
        if sample['DP'] is not None and sample['DP'] >= args.dp:
            # Ensure that AD (Alternate Depth) is not None and contains at least 2 values (ref and alt alleles)
            if sample['AD'] is not None and len(sample['AD']) > 1:
                # Get the depth of the alternate allele (AD[1])
                alt_reads = sample['AD'][1]  # AD[1] is the count of alternate allele reads
                # Ensure that alt_reads is not None and meets the minimum SNP depth (args.snpdp)
                if alt_reads is not None and alt_reads >= args.snpdp:
                    # Check if the sample has a valid genotype (GT) and if the genotype is heterozygous or homozygous alternate
                    gt_tuple = tuple(sample['GT'])  # Convert the genotype (GT) to a tuple for easy comparison
                    if gt_tuple == (1, 0) or gt_tuple == (0, 1) or gt_tuple == (1, 1):
                        # Check if the sample has an allele frequency (AF) that meets the minimum requirement (args.af)
                        if sample['AF'] is not None and sample['AF'] >= args.af:
                            return True  # All conditions are satisfied, so this record should be included

    # If none of the samples meet the criteria, return False to exclude the record
    return False

def main():
    # Parse arguments
    args = parse_args()

    # Load BED regions (if specified)
    bed_regions = set()
    if args.bed:
        bed_regions = load_bed_regions(args.bed)

    # Open input VCF file and output VCF file
    input_vcf = pysam.VariantFile(args.vcf_file, 'r')  # Read input VCF file
    output_vcf = pysam.VariantFile(args.output, 'w', header=input_vcf.header)  # Write output VCF file

    # Iterate over each record in the input VCF and filter them
    for record in input_vcf:
        # If the record passes the filter, write it to the output VCF
        if filter_record(record, args, bed_regions):
            output_vcf.write(record)

    # Close both input and output VCF files
    input_vcf.close()
    output_vcf.close()

# If this script is being run directly, invoke the main function
if __name__ == '__main__':
    main()