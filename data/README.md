# Input Data Directory

This directory should contain your input data files for GWAS analysis.

## Required Files

### 1. Genotype Data
Place your genotype file(s) here in one of the supported formats:

- **HapMap format** (`.hmp.txt`) - Recommended
- **VCF format** (`.vcf`)
- **Numeric format** (`.txt`)
- **Plink format** (`.ped` and `.map`)

### 2. Phenotype Data
Place your phenotype file(s) here:

- Tab-delimited (`.txt`) or CSV (`.csv`) format
- First column: Sample/Taxa names (must match genotype file)
- Subsequent columns: Trait values (numeric)

## Example File Structure

```
data/
├── genotype.hmp.txt          # Genotype data in HapMap format
├── phenotype.txt             # Phenotype data with trait values
└── README.md                 # This file
```

## Data Format Details

### HapMap Format Example
```
rs#    alleles    chrom    pos    strand    assembly#    Sample1    Sample2    Sample3
SNP1   A/G        1        1234   +         NA           AA         AG         GG
SNP2   C/T        1        5678   +         NA           CC         CT         TT
```

### Phenotype Format Example
```
Taxa       EarHeight    EarDiameter    DaysToPollen
Sample1    45.2         3.5            65
Sample2    52.1         4.1            68
Sample3    48.7         3.8            67
```

## Notes

- Ensure sample names are consistent between genotype and phenotype files
- Remove missing data or code as `NA`
- For best results, filter SNPs with MAF < 0.05 before analysis
- Phenotype data should be normally distributed (transform if needed)
