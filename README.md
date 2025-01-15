1. snakemake pipeline for variant calling in ancient pathogen samples using mpileup
   produces multivcf
   adds custom field called [AF] allele frequency
2. python script for filtering multivcf output

exp python BED-OG-variable.py AFmergenormformat.vcf.gz -d 5 -a .10 -n 3 -v -b /data/stonelab/Kelly_TB/mpileup/MTBC_regions_to_exclude.bed -o SNP5.13TLTpreBlast.vcf

'-s', '--sample', help='comma-separated list of sample names to filter; used to ignore sites only variable in outgroup'
'-d', '--dp', type=int, required=True, help='minimum overall read depth (DP) to retain a variant'
'-n', '--snpdp', type=int, required=True, help='minimum reads supporting the SNP (AD) to retain a variant'
'-a', '--af', type=float, required=True, help='minimum allele frequency (AF) to retain a variant'
'-b', '--bed', help='input BED file to specify regions to filter/mask'
'-o', '--output', required=True, help='output VCF file'
'-v', '--snp', action='store_true', help='only include SNPs (exclude indels)'
