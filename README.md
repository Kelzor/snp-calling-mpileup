## 1. Snakemake Pipeline for Variant Calling in Ancient Pathogen Samples Using mpileup
This pipeline produces a **multivcf** file and adds a custom field called `[AF]` for allele frequency.

The produced vcf files have a record for every site in the reference genome. This is important logic for combining ancient DNA VCFs. Ancient DNA data are low coverage. It therefore cannot be assumed that a position is a reference call if a variant is not called there.

---

## 2. Python Script for Filtering Multivcf Output

Parses multivcf and outputs only variant sites as defined by user options.

**Important**: if one sample passes the criteria, all sample records are kept, even if they do not themselves pass the filters. 

May add functionality to convert these "failures" to Ns?

Use the following Python script to filter the multiVCF file produced by mpileup:

### Command-Line Arguments  

| Argument        | Short | Type    | Required | Default | Description |
|----------------|-------|---------|----------|---------|-------------|
| `vcf_file`     | N/A   | File    | ✅ Yes   | N/A     | Input VCF file to be filtered. |
| `--sample`     | `-s`  | String  | ❌ No    | N/A     | Comma-separated list of sample names to filter out. This option ignores sites that are ONLY VARIABLE in the outgroup. |
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

```

## 4. Python Script for converting multi-sample VCF into SNP alignment FASTA

### Command-Line Arguments

| Argument       | Short | Type    | Required | Default   | Description |
|----------------|-------|---------|----------|-----------|-------------|
| `--vcf`        | `-i`  | File    | ✅ Yes   | N/A       | Input multi-sample VCF file. |
| `--output`     | `-o`  | File    | ✅ Yes   | N/A       | Output FASTA file containing the SNP alignment. A debug file with the same name and `.debug.txt` extension will also be created. |
| `--min_af`     |       | Float   | ❌ No    | `0.2`     | Minimum allele frequency threshold below which the reference allele is chosen. |
| `--max_af`     |       | Float   | ❌ No    | `0.8`     | Maximum allele frequency threshold above which the alternate allele is chosen. |
| `--depth`      |       | Integer | ❌ No    | `10`      | Minimum depth threshold required for a base call. |
| `--deletion`   |       | Float   | ❌ No    | `None`    | Post-alignment filter: removes sites where fewer than this proportion of samples have data (non-`N`). Value must be between `0.1` and `1.0`. |
| `--no-debug`   |       | Flag   | ❌ No    | `False`    | If set, the full debug file is not written. Instead, a small summary file with total SNPs and number of samples is generated, and the information is printed to the console. This is to save time for big alignments.|

The debug file is a tab-delimited text file that logs decisions made for each variant and sample.
Each line corresponds to one sample at one SNP site.

If `--deletion` is specified, only positions retained in the final filtered alignment are written to the debug file.

INDELs are skipped — only SNPs are included.
Positions with insufficient depth (DP < --depth) or ambiguous allele frequency (--min_af < AF < --max_af) are called as N.
Missing genotypes (./.) are also represented as N.

---

### Example Usage  

Convert a VCF to a SNP alignment, keeping only sites where at least 90% of samples are present (non-missing):  

```bash
python MSAv8.py -i SNP5.13TLTpostMalt.vcf -o v8alignment9PD.fasta --min_af 0.1 --max_af 0.8 --depth 3 --deletion 0.9
```
or
```bash
python MSAv8.py -i SNP5.13TLTpostMalt.vcf -o v8alignment9PD.fasta --min_af 0.1 --max_af 0.8 --depth 3
```
