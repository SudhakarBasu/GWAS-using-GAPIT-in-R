# GAPIT Output Guide

This document provides detailed explanations of all GAPIT output files.

## Association Results (`results/association/`)

### GWAS Results Files
**Format**: `GAPIT.Association.GWAS_Results.[MODEL].[TRAIT].csv`

Contains the main association test results for each SNP.

**Key Columns**:
- `SNP`: Marker identifier
- `Chromosome`: Chromosome number
- `Position`: Physical position in base pairs
- `P.value`: Statistical significance (lower = stronger association)
- `maf`: Minor allele frequency
- `nobs`: Number of observations
- `Rsquare.of.Model.without.SNP`: Model R² without this SNP
- `Rsquare.of.Model.with.SNP`: Model R² with this SNP
- `effect`: Estimated allele effect size

**Interpretation**: SNPs with P.value < 1e-5 are typically considered significant. The `effect` column shows the phenotypic change per allele substitution.

### Standard Error Files
**Format**: `GAPIT.Association.GWAS_StdErr.[MODEL].[TRAIT].csv`

Standard errors for effect size estimates. Used to calculate confidence intervals.

### Prediction Results
**Format**: `GAPIT.Association.Prediction_results.[MODEL].[TRAIT].csv`

Genomic prediction accuracy and predicted breeding values for each individual.

**Key Columns**:
- `Taxa`: Sample identifier
- `Observed`: Actual phenotype value
- `Predicted`: Genomic predicted value
- `Residual`: Observed - Predicted

**Interpretation**: Correlation between Observed and Predicted indicates prediction accuracy.

### PVE (Proportion of Variance Explained)
**Format**: `GAPIT.Association.PVE.[MODEL].[TRAIT].csv`

Shows the percentage of phenotypic variance explained by significant markers.

### Filtered Results
**Format**: `GAPIT.Association.Filter_GWAS_results.csv`

Contains only SNPs passing significance thresholds, useful for downstream analysis.

---

## Visualization Outputs

### Manhattan Plots (`results/manhattan_plots/`)

**Files**:
- `Manhattan_Chro.[MODEL].[TRAIT].pdf`: Chromosome-based layout
- `Manhattan_Geno.[MODEL].[TRAIT].pdf`: Genome-wide view
- `Manhattans_Circular.pdf`: Circular genome representation
- `Manhattans_Symphysic.pdf`: Multi-trait comparison

**How to Read**:
- X-axis: Genomic position (chromosome and position)
- Y-axis: -log10(P-value)
- Horizontal line: Significance threshold (typically -log10(5e-8) ≈ 7.3)
- Peaks above threshold: Significant associations

**What to Look For**:
- Sharp peaks indicate strong associations
- Multiple peaks in same region suggest linkage disequilibrium
- Peaks across genome suggest polygenic trait

### QQ Plots (`results/qq_plots/`)

**Files**: `GAPIT.Association.QQplot.[MODEL].[TRAIT].pdf`

**How to Read**:
- X-axis: Expected -log10(P-value) under null hypothesis
- Y-axis: Observed -log10(P-value)
- Diagonal line: Perfect fit (no association)
- Deviation from diagonal: Indicates true associations or model issues

**Interpretation**:
- Points on diagonal: Good model fit, population structure controlled
- Early deviation: Inflation due to population structure (use MLM instead of GLM)
- Late deviation only: True associations present

### PCA Plots (`results/pca/`)

**Files**: 
- `GAPIT.Genotype.PCA.csv`: Principal component scores
- `GAPIT.Genotype.PCA_2D.pdf`: 2D scatter plots
- `GAPIT.Genotype.PCA_3D.pdf`: 3D visualization

**Purpose**: Visualize population structure and genetic relationships.

**Interpretation**: Clusters indicate subpopulations. PCs are used as covariates in association models.

---

## Genotype Analysis (`results/genotype/`)

### Kinship Matrix
**Files**: 
- `GAPIT.Genotype.Kin_Zhang.csv`: Kinship coefficients
- `GAPIT.Genotype.Kin_Zhang.pdf`: Heatmap visualization

**Purpose**: Measures genetic relatedness between individuals.

**Values**:
- 1.0: Identical (same individual or clone)
- 0.5: Parent-offspring or full siblings
- 0.25: Half-siblings or grandparent-grandchild
- 0.0: Unrelated

### MAF Distribution
**Files**:
- `GAPIT.Genotype.Frequency_MAF.csv`: MAF for each SNP
- `GAPIT.Genotype.Frequency.pdf`: Histogram
- `GAPIT.Genotype.MAF_Heterozosity.pdf`: MAF vs heterozygosity

**Purpose**: Quality control and filtering.

**Interpretation**: Most SNPs should have MAF > 0.05 for reliable association testing.

### Linkage Disequilibrium
**Files**:
- `GAPIT.Genotype.Distance.Rsquare.csv`: LD decay statistics
- `GAPIT.Genotype.Distance_R_Chro.pdf`: LD decay plot

**Purpose**: Understand LD patterns, affects GWAS resolution.

**Interpretation**: Faster LD decay = higher mapping resolution.

---

## Phenotype Analysis (`results/phenotype/`)

### Distribution Plots
**Files**: 
- `GAPIT.Phenotype.View.[TRAIT].pdf`: Histogram and summary statistics
- `GAPIT.Phenotype.Distribution_Significantmarkers.[MODEL].[TRAIT].pdf`: Phenotype distribution by genotype at significant SNPs

**Purpose**: 
- Check for normality (important for GLM/MLM)
- Identify outliers
- Visualize trait variation

**What to Look For**:
- Normal distribution (bell curve)
- No extreme outliers
- Sufficient variation for association testing

---

## Model Comparison

### GLM vs MLM

**GLM (General Linear Model)**:
- Faster computation
- Suitable for unrelated individuals
- May have false positives with population structure

**MLM (Mixed Linear Model)**:
- Accounts for kinship and population structure
- More conservative (fewer false positives)
- Slower computation
- Recommended for most GWAS

**When to Use Each**:
- Use GLM first for quick exploration
- Use MLM for final analysis and publication
- Compare QQ plots to assess model fit

---

## Tips for Interpretation

1. **Always check QQ plots first** - Ensures model is appropriate
2. **Look for consistency across models** - True associations appear in both GLM and MLM
3. **Examine LD around peaks** - Significant SNPs cluster due to LD
4. **Consider biological context** - Are associations near known candidate genes?
5. **Validate top hits** - Replicate in independent populations if possible

---

## Common Issues

**Problem**: QQ plot shows early deviation from diagonal  
**Solution**: Use MLM instead of GLM, or add more PCs

**Problem**: No significant associations  
**Solution**: Check sample size, trait heritability, marker density

**Problem**: Too many significant SNPs  
**Solution**: May indicate population structure; use stricter threshold or MLM

**Problem**: Kinship matrix has negative values  
**Solution**: Normal for VanRaden method; indicates less related than average

---

For more details, see the main [README.md](../README.md) or visit the [GAPIT documentation](http://zzlab.net/GAPIT/).
