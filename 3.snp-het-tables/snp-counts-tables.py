import argparse
import pysam
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze SNPs from a multi-sample VCF file.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file.")
    parser.add_argument("-o", "--output", type=str, default="snp_summary.tsv", help="Output file for the SNP summary table.")
    parser.add_argument("-d", "--dp", type=int, required=True, help="Minimum depth of coverage (DP).")
    parser.add_argument("-ad", "--allele_depth", type=int, required=True, help="Minimum allele depth (AD).")
    parser.add_argument("-a", "--af", type=float, required=True, help="Minimum allele frequency (AF) to classify as a SNP.")
    parser.add_argument("-homa", "--homo_af", type=float, required=True, help="Minimum allele frequency for homozygous SNPs.")
    parser.add_argument("-s", "--sample", nargs="*", default=[], help="Samples to exclude from analysis.")
    return parser.parse_args()


def process_vcf(vcf_file, min_depth, min_allele_depth, min_af, homo_af, exclude_samples):
    sample_data = {}
    vcf = pysam.VariantFile(vcf_file)

    # Initialize the sample data
    for sample_name in vcf.header.samples:
        if sample_name not in exclude_samples:
            sample_data[sample_name] = {
                "snps": 0,
                "homozygous_snps": 0,
                "heterozygous_snps": 0,
                "positions": []  # To store SNP positions, AF, DP, and AD for each sample
            }

    # Process each record in the VCF file
    for record in vcf:
        for sample_name in record.samples.keys():
            if sample_name in exclude_samples:
                continue

            sample = record.samples[sample_name]

            # Apply DP filter
            if sample.get("DP") is not None and sample["DP"] >= min_depth:
                # Apply AD filter (minimum allele depth)
                ad_values = sample.get("AD", None)
                if ad_values is not None and isinstance(ad_values, (list, tuple)) and len(ad_values) > 1:
                    allele_depth = sum(x for x in ad_values[1:] if x is not None)  # Filter out None values
                    if allele_depth < min_allele_depth:
                        continue  # Skip if allele depth is below the threshold

                    # Apply GT and AF filters
                    gt_tuple = tuple(sample.get("GT", []))
                    af = sample.get("AF", None)

                    if af is not None and isinstance(af, float) and gt_tuple in [(1, 0), (0, 1), (1, 1)]:
                        if af >= min_af:
                            sample_data[sample_name]["snps"] += 1
                            if af >= homo_af:
                                sample_data[sample_name]["homozygous_snps"] += 1
                            elif min_af <= af < homo_af:
                                sample_data[sample_name]["heterozygous_snps"] += 1

                        # Save SNP position, AF (formatted), DP, and AD for each SNP
                        sample_data[sample_name]["positions"].append((record.pos, format(af, ".4f"), sample["DP"], allele_depth))

    return sample_data


def write_output(output_file, results, min_af, homo_af, min_depth, min_allele_depth):
    with open(output_file, "w") as out:
        # Write the summary table headers
        out.write("Sample\tTotal SNPs\tHomozygous SNPs\tHeterozygous SNPs\n")
        
        # Write the results for each sample
        for sample, counts in results.items():
            out.write(f"{sample}\t{counts['snps']}\t{counts['homozygous_snps']}\t{counts['heterozygous_snps']}\n")
        
        # Add a summary of the allele frequency thresholds and SNP depth used
        out.write(f"\nAllele frequency thresholds used:\n")
        out.write(f"SNP AF >= {min_af:.4f}, Homozygous SNP AF >= {homo_af:.4f}, Heterozygous SNP {min_af:.4f} <= AF < {homo_af:.4f}\n")
        out.write(f"Read Depth used: {min_depth}, Minimum Allele Depth used: {min_allele_depth}\n")


def write_sample_snp_details(sample_data, output_folder):
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Write individual tables for each sample
    for sample, data in sample_data.items():
        snp_details_file = f"{output_folder}/{sample}_snp_details.tsv"
        with open(snp_details_file, "w") as snp_file:
            snp_file.write("Position\tAF\tDepth\tAllele Depth\n")
            for pos, af, dp, ad in data["positions"]:
                snp_file.write(f"{pos}\t{af}\t{dp}\t{ad}\n")
        print(f"SNP details for {sample} written to {snp_details_file}")


def main():
    args = parse_args()

    # Process the VCF file
    results = process_vcf(
        args.input, 
        args.dp, 
        args.allele_depth,
        args.af, 
        args.homo_af, 
        args.sample
    )

    # Write the results to the summary output file
    write_output(args.output, results, args.af, args.homo_af, args.dp, args.allele_depth)

    # Write individual SNP tables for each sample
    write_sample_snp_details(results, output_folder="snp_details")

    print(f"SNP summary written to {args.output}")


if __name__ == "__main__":
    main()
