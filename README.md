# GWAS Analysis with GAPIT

A comprehensive repository for performing Genome-Wide Association Studies (GWAS) using the GAPIT (Genome Association and Prediction Integrated Tool) R package.

## 📋 Overview

**Genome-Wide Association Studies (GWAS)** are a powerful approach to identify genetic variants associated with traits of interest by scanning markers across complete genomes of many individuals. This repository contains scripts, results, and documentation for conducting GWAS analyses using GAPIT.

**GAPIT** is a comprehensive R package that provides a user-friendly interface for performing GWAS with multiple statistical models, including:
- General Linear Model (GLM)
- Mixed Linear Model (MLM)
- BLINK (Bayesian-information and Linkage-disequilibrium Iteratively Nested Keyway)
- FarmCPU (Fixed and random model Circulating Probability Unification)

## 🗂️ Repository Structure

```
GWAS_GAPIT/
├── README.md                    # This file
├── scripts/                     # R scripts for GWAS analysis
│   ├── GAPIT3_IIRR.R           # Main analysis script
│   └── GAPIT.R                 # GAPIT package functions
├── data/                        # Input data directory (genotype & phenotype)
├── results/                     # Analysis outputs
│   ├── association/            # Association test results
│   ├── manhattan_plots/        # Manhattan plot visualizations
│   ├── pca/                    # Principal Component Analysis results
│   ├── qq_plots/               # QQ plot visualizations
│   ├── genotype/               # Genotype analysis (MAF, kinship, LD)
│   └── phenotype/              # Phenotype distribution plots
└── docs/                        # Additional documentation
    └── output_guide.md         # Detailed guide to GAPIT outputs
```

## 🔧 Installation & Requirements

### R Version
- R >= 3.5.0 (recommended: R >= 4.0.0)

### Required R Packages

```r
# Install required packages
install.packages(c("multtest", "gplots", "LDheatmap", "genetics", 
                   "ape", "EMMREML", "scatterplot3d"))

# Install GAPIT from GitHub
install.packages("devtools")
devtools::install_github("jiabowang/GAPIT3", force=TRUE)
```

### Core Dependencies
- **multtest**: Multiple hypothesis testing
- **gplots**: Enhanced plotting functions
- **LDheatmap**: Linkage disequilibrium visualization
- **genetics**: Population genetics tools
- **ape**: Phylogenetic analysis
- **EMMREML**: Mixed model analysis
- **scatterplot3d**: 3D visualization

## 📊 Input Data Format

### Genotype Data
GAPIT accepts multiple genotype formats:
- **HapMap format** (recommended): Tab-delimited text file with SNP information
- **VCF format**: Variant Call Format
- **Numeric format**: 0/1/2 coding for SNP genotypes
- **Plink format**: .ped and .map files

Example HapMap format:
```
rs#    alleles    chrom    pos    strand    assembly#    center    protLSID    assayLSID    panelLSID    QCcode    Sample1    Sample2    ...
SNP1   A/G        1        1234   +         NA           NA        NA          NA           NA           NA        AA         AG         ...
```

### Phenotype Data
Tab-delimited or CSV file with:
- **Column 1**: Taxa/Sample names (must match genotype file)
- **Column 2+**: Trait values (numeric)

Example:
```
Taxa       EarHT    EarDia    dpoll
Sample1    45.2     3.5       65
Sample2    52.1     4.1       68
```

## 🚀 Usage

### Basic GWAS Analysis

```r
# Load GAPIT
library(GAPIT3)

# Set working directory
setwd("path/to/GWAS_GAPIT")

# Load your data
myY <- read.table("data/phenotype.txt", head = TRUE)
myG <- read.table("data/genotype.hmp.txt", head = FALSE)

# Run GAPIT with GLM model
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  model = "GLM",
  PCA.total = 3
)

# Run GAPIT with MLM model (accounts for population structure)
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  model = "MLM",
  PCA.total = 3
)
```

### Advanced Options

```r
# Multiple models with kinship matrix
myGAPIT <- GAPIT(
  Y = myY,
  G = myG,
  model = c("GLM", "MLM", "BLINK", "FarmCPU"),
  PCA.total = 5,
  kinship.algorithm = "VanRaden",
  Multiple_analysis = TRUE
)
```

## 📈 Understanding GAPIT Outputs

### Association Results (`results/association/`)

| File | Description |
|------|-------------|
| `GWAS_Results.*.csv` | SNP-trait associations with p-values, effect sizes, and statistics |
| `GWAS_StdErr.*.csv` | Standard errors for effect estimates |
| `Prediction_results.*.csv` | Genomic prediction accuracy and breeding values |
| `PVE.*.csv` | Proportion of Variance Explained by significant markers |
| `Filter_GWAS_results.csv` | Significant SNPs passing threshold filters |

**Key columns in GWAS_Results:**
- **SNP**: Marker identifier
- **Chromosome**: Chromosome number
- **Position**: Physical position (bp)
- **P.value**: Statistical significance
- **maf**: Minor allele frequency
- **effect**: Allele effect size

### Visualization Outputs

#### Manhattan Plots (`results/manhattan_plots/`)
- Display -log10(p-values) across the genome
- Peaks indicate significant trait-associated regions
- Horizontal line shows significance threshold (typically p < 0.05 after Bonferroni correction)

#### QQ Plots (`results/qq_plots/`)
- Compare observed vs. expected p-value distribution
- Assess model fit and population structure control
- Deviation from diagonal indicates inflation/deflation

#### PCA Results (`results/pca/`)
- Principal components for population structure
- Used as covariates in association models
- Visualize genetic relationships among samples

### Genotype Analysis (`results/genotype/`)

| File | Description |
|------|-------------|
| `Kin_Zhang.csv` | Kinship matrix (genetic relatedness) |
| `Frequency_MAF.csv` | Minor allele frequency distribution |
| `Distance.Rsquare.csv` | Linkage disequilibrium decay |

### Phenotype Analysis (`results/phenotype/`)
- Distribution plots for each trait
- Identifies outliers and data quality issues
- Shows trait variation across samples

## 🔬 Example Workflow

This repository contains results from analyzing three maize traits:
1. **EarHT**: Ear height (cm)
2. **EarDia**: Ear diameter (cm)
3. **dpoll**: Days to pollen

### Analysis Steps

1. **Data Preparation**
   - Quality control on genotype data (MAF > 0.05)
   - Phenotype normalization and outlier detection

2. **Population Structure Analysis**
   - PCA to capture population stratification
   - Kinship matrix calculation

3. **Association Testing**
   - GLM: Fast, suitable for unrelated individuals
   - MLM: Accounts for relatedness and population structure

4. **Results Interpretation**
   - Identify significant SNPs (p < 1e-5)
   - Examine Manhattan and QQ plots
   - Validate top associations

## 📚 Key Concepts

### P-value
Statistical significance of SNP-trait association. Lower values indicate stronger evidence.

### MAF (Minor Allele Frequency)
Frequency of the less common allele. Rare variants (MAF < 0.05) often excluded due to low power.

### Effect Size
Magnitude of phenotypic change per allele substitution.

### PVE (Proportion of Variance Explained)
Percentage of trait variation explained by a SNP or set of SNPs.

### Kinship Matrix
Measures genetic relatedness between individuals. Controls for family structure in MLM.

### Linkage Disequilibrium (LD)
Non-random association between alleles at different loci. Affects resolution of GWAS.

## 📖 References

### GAPIT Publications
- **GAPIT Version 1**: Lipka et al. (2012). *Bioinformatics* 28(18):2397-2399
- **GAPIT Version 2**: Tang et al. (2016). *The Plant Genome* 9(2)
- **GAPIT Version 3**: Wang & Zhang (2021). *Genomics, Proteomics & Bioinformatics* 19(4):629-640

### Useful Resources
- [GAPIT GitHub Repository](https://github.com/jiabowang/GAPIT)
- [GAPIT Tutorial](http://zzlab.net/GAPIT/)
- [GWAS Catalog](https://www.ebi.ac.uk/gwas/)

## 🤝 Contributing

For questions or contributions, please open an issue or submit a pull request.

## 📄 License

This project is open source and available for research and educational purposes.

---

**Last Updated**: January 2026
