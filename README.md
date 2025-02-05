### 1. Snakemake Pipeline for Variant Calling in Ancient Pathogen Samples Using mpileup
This pipeline produces a **multivcf** file and adds a custom field called `[AF]` for allele frequency.

---

### 2. Python Script for Filtering Multivcf Output

Parses multivcf and outputs only variant sites as defined by user options.
**Important**: if one sample passes the criteria, all sample records are kept, even if they do not themselves pass the filters. 

May add functionality to convert these "failures" to N/As?

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

### 3. Python Script for counting SNPs (overall, homozygous, and het) based on user defined critera
## Command-Line Arguments  

| Argument        | Short | Type    | Required | Default             | Description |
|----------------|-------|---------|----------|---------------------|-------------|
| `--input`      | `-i`  | File    | ✅ Yes   | N/A                 | Input VCF file. |
| `--output`     | `-o`  | File    | ❌ No    | `snp_summary.tsv`   | Output file for the SNP summary table. |
| `--dp`         | `-d`  | Integer | ✅ Yes   | N/A                 | Minimum depth of coverage (DP) required for a record to be considered. |
| `--allele_depth` | `-ad` | Integer | ✅ Yes   | N/A                 | Minimum allele depth (AD) required for a record to be counted as a SNP. |
| `--af`         | `-a`  | Float   | ✅ Yes   | N/A                 | Minimum allele frequency (AF) to classify a variant as a SNP. |
| `--homo_af`    | `-homa` | Float | ✅ Yes   | N/A                 | Minimum allele frequency for homozygous SNP classification. |
| `--sample`     | `-s`  | List    | ❌ No    | `[]` (Empty List)   | List of samples to exclude from analysis. |

---

### Example Usage  

```bash
python snp-counts.py -i input.vcf -o output.tsv -d 10 -ad 5 -a 0.2 -homa 0.8 -s sample1 sample2
