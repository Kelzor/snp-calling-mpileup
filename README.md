### 1. Snakemake Pipeline for Variant Calling in Ancient Pathogen Samples Using mpileup
This pipeline produces a **multivcf** file and adds a custom field called `[AF]` for allele frequency.

---

### 2. Python Script for Filtering Multivcf Output
Use the following Python script to filter the output:

```bash
python BED-OG-variable.py AFmergenormformat.vcf.gz -d 5 -a .10 -n 3 -v -b /data/stonelab/Kelly_TB/mpileup/MTBC_regions_to_exclude.bed -o SNP5.13TLTpreBlast.vcf
```

### Available Options:

- **`-s`** `--sample`: A comma-separated list of sample names to filter. This option helps ignore sites that are only variable in the outgroup.
- **`-d`** `--dp` (type=int): Minimum overall read depth (DP) to retain a variant.
- **`-n`** `--snpdp` (type=int, required): Minimum reads supporting the SNP (AD) to retain a variant.
- **`-a`** `--af` (type=float, required): Minimum allele frequency (AF) to retain a variant.
- **`-b`** `--bed`: Input BED file to specify regions to filter/mask.
- **`-o`** `--output` (required): Output VCF file.
- **`-v`** `--snp`: Only include SNPs (exclude indels).
