## 1. Snakemake Pipeline for Variant Calling in Ancient Pathogen Samples Using mpileup
This pipeline produces a **multivcf** file and adds a custom field called `[AF]` for allele frequency.

---

## 2. Python Script for Filtering Multivcf Output

Parses multivcf and outputs only variant sites as defined by user options.

**Important**: if one sample passes the criteria, all sample records are kept, even if they do not themselves pass the filters. 

May add functionality to convert these "failures" to N/As?

Use the following Python script to filter the output:

### Command-Line Arguments  

| Argument        | Short | Type    | Required | Default | Description |
|----------------|-------|---------|----------|---------|-------------|
| `vcf_file`     | N/A   | File    | ✅ Yes   | N/A     | Input VCF file to be filtered. |
| `--sample`     | `-s`  | String  | ❌ No    | N/A     | Comma-separated list of sample names to filter out. This option ignores sites that are ONLY VARIANLE in the outgroup. |
| `--dp`         | `-d`  | Integer | ✅ Yes   | N/A     | Minimum read depth (DP) to retain a variant. |
| `--snpdp`      | `-n`  | Integer | ✅ Yes   | N/A     | Minimum SNP depth (alternate allele depth in AD) to retain a variant. |
| `--af`         | `-a`  | Float   | ✅ Yes   | N/A     | Minimum allele frequency (AF) to retain a variant. |
| `--bed`        | `-b`  | File    | ❌ No    | N/A     | Input BED file specifying regions for filtering. |
| `--output`     | `-o`  | File    | ✅ Yes   | N/A     | Output VCF file after filtering. |
| `--snp`        | `-v`  | Flag    | ❌ No    | `False` | Only include SNPs (exclude indels). |

### Example Usage  
```bash
python BED-OG-variable.py AFmergenormformat.vcf.gz -d 5 -a .10 -n 3 -v -b /data/stonelab/Kelly_TB/mpileup/MTBC_regions_to_exclude.bed -o SNP5.13TLTpreBlast.vcf

```
---

## 3. Python Script for counting SNPs (overall, homozygous, and het) based on user defined critera

### Command-Line Arguments  

| Argument        | Short | Type    | Required | Default             | Description |
|----------------|-------|---------|----------|---------------------|-------------|
| `--input`      | `-i`  | File    | ✅ Yes   | N/A                 | Input VCF file. |
| `--output`     | `-o`  | File    | ❌ No    | `snp_summary.tsv`   | Output file for the SNP summary table. |
| `--dp`         | `-d`  | Integer | ✅ Yes   | N/A                 | Minimum depth of coverage (DP) required for a record to be considered. |
| `--allele_depth` | `-ad` | Integer | ✅ Yes   | N/A                 | Minimum allele depth (AD) required for a record to be counted as a SNP. |
| `--af`         | `-a`  | Float   | ✅ Yes   | N/A                 | Minimum allele frequency (AF) to classify a variant as a SNP. |
| `--homo_af`    | `-homa` | Float | ✅ Yes   | N/A                 | Minimum allele frequency for homozygous SNP classification. |
| `--sample`     | `-s`  | List    | ❌ No    | N/A   | List of samples to exclude from analysis. |

---

### Example Usage  

```bash
python snp-counts.py -i input.vcf -o output.tsv -d 10 -ad 5 -a 0.2 -homa 0.8 -s sample1 sample2
