import argparse
import pysam

def parse_args():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Filter VCF file based on DP, AF, and GT')
    parser.add_argument('vcf_file', help='input VCF file')
    parser.add_argument('-s', '--sample', help='comma-separated list of sample names to filter')
    parser.add_argument('-d', '--dp', type=int, required=True, help='minimum read depth (DP) to retain a variant')
    parser.add_argument('-a', '--af', type=float, required=True, help='minimum allele frequency (AF) to retain a variant')
    parser.add_argument('-b', '--bed', help='input BED file to specify regions to filter')
    parser.add_argument('-o', '--output', required=True, help='output VCF file')
    args = parser.parse_args()

    # If sample names are provided, split comma-separated sample names into a list
    if args.sample:
        args.sample = args.sample.split(',')
    else:
        args.sample = []  # Set an empty list if sample names are not provided

    return args

def load_bed_regions(bed_file):
    bed_regions = set()
    with open(bed_file, 'r') as bed:
        for line in bed:
            parts = line.strip().split('\t')
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            for position in range(start, end + 1):
                bed_regions.add((chrom, position))
    return bed_regions

def filter_record(record, args, bed_regions):
    # Check if the variant's position is in the BED regions
    if (record.chrom, record.pos) in bed_regions:
        return False

    # Loop through all samples in the record
    for sample_name in record.samples.keys():
        # Skip the user-specified samples
        if sample_name in args.sample:
            continue
        sample = record.samples[sample_name]
        # Check if the sample passes the DP and AF filters
        if sample['DP'] is not None and sample['DP'] >= args.dp:
            # Check if the sample passes the GT filter
            gt_tuple = tuple(sample['GT'])
            if gt_tuple == (1, 0) or gt_tuple == (0, 1) or gt_tuple == (1, 1):
                if sample['AF'] is not None and sample['AF'] >= args.af:
                    return True
    return False

def main():
    args = parse_args()

    # Load BED regions
    bed_regions = set()
    if args.bed:
        bed_regions = load_bed_regions(args.bed)

    # Open input and output VCF files
    input_vcf = pysam.VariantFile(args.vcf_file, 'r')
    output_vcf = pysam.VariantFile(args.output, 'w', header=input_vcf.header)

    # Filter and write output VCF file
    for record in input_vcf:
        if filter_record(record, args, bed_regions):
            output_vcf.write(record)

    # Close input and output VCF files
    input_vcf.close()
    output_vcf.close()

if __name__ == '__main__':
    main()

