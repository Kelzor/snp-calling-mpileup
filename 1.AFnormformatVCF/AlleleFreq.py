import argparse
import pysam


def add_af_field(input_vcf, output_vcf):
    # Read the input file
    myvcf = pysam.VariantFile(input_vcf, "r")

    # Add the AF field to header.
    myvcf.header.formats.add('AF', '1', 'Float', 'Allele frequency calculated from AD/DP')

    # Create an object of the new vcf file and open it to write data.
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=myvcf.header)

    for variant in myvcf:
        for sample in variant.samples:
            print(variant.samples[sample]['GT'])
            if variant.samples[sample]['GT'] == (0, 1) or variant.samples[sample]['GT'] == (1, 0) or variant.samples[sample]['GT'] == (1, 1):
                ad = variant.samples[sample]['AD']
                #print(variant.pos)
                print(ad)
                ad = ad[1]
                dp = variant.samples[sample]['DP']
                af = ad / dp
                variant.samples[sample]["AF"] = af
            elif variant.samples[sample]['GT'] == (0, 2) or variant.samples[sample]['GT'] == (2, 2):
                ad = variant.samples[sample]['AD']
                #print(variant.pos)
                print(ad)
                ad = ad[2]
                dp = variant.samples[sample]['DP']
                af = ad / dp
                variant.samples[sample]["AF"] = af
        vcf_out.write(variant)

    vcf_out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add AF field to VCF file.')
    parser.add_argument('input_vcf', help='Input VCF file')
    parser.add_argument('output_vcf', help='Output VCF file with added AF field')

    args = parser.parse_args()

    add_af_field(args.input_vcf, args.output_vcf)

