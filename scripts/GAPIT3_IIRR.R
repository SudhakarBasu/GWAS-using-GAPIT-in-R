################################################################################
# GWAS Analysis with GAPIT3
#
# This script performs Genome-Wide Association Studies (GWAS) using the GAPIT3
# R package. It includes installation instructions, data loading, and analysis
# examples with different models and parameters.
#
# Author: [Your Name]
# Date: January 2026
# GAPIT Version: 3.x
################################################################################

# ===============================================================================
# SECTION 1: PACKAGE INSTALLATION
# ===============================================================================

# GAPIT3 can be installed in multiple ways. Choose ONE of the following methods:

# METHOD 1: Install from GitHub (Recommended)
# This installs the latest version from the official repository
if (!require("remotes")) install.packages("remotes")
remotes::install_github("jiabowang/GAPIT3", force = TRUE)
library(GAPIT3)

# METHOD 2: Source from ZZlab website
# Alternative method that loads functions directly from the web
# source("http://zzlab.net/GAPIT/GAPIT.library.R")
# source("http://zzlab.net/GAPIT/gapit_functions.txt")

# METHOD 3: Install from local archive file
# If you have downloaded a .tar.gz or .zip file
# install.packages("GAPIT3_3.1.0.9000.tar.gz", repos = NULL, type = "source")

# ===============================================================================
# SECTION 2: SET WORKING DIRECTORY
# ===============================================================================

# Set the working directory to where your data files are located
# IMPORTANT: Update this path to match your data location
setwd("path/to/your/data") # Change this to your actual data directory

# Alternative: Use RStudio's GUI
# Session -> Set Working Directory -> Choose Directory

# Verify current working directory
getwd()

# ===============================================================================
# SECTION 3: LOAD INPUT DATA
# ===============================================================================

# GAPIT requires two main input files:
# 1. Genotype file (HapMap format, VCF, or numeric)
# 2. Phenotype file (tab-delimited with taxa names and trait values)

# Load genotype data (HapMap format)
# HapMap format: header = FALSE
# First 11 columns: SNP info (rs#, alleles, chrom, pos, etc.)
# Remaining columns: Genotype calls for each sample
myGenotype <- read.table("genotype.hmp.txt", header = FALSE)

# Load phenotype data
# Column 1: Taxa/Sample names (must match genotype file)
# Column 2+: Trait values (numeric)
myPhenotype <- read.table("phenotype.txt", header = TRUE)

# ===============================================================================
# SECTION 4: DATA QUALITY CHECKS
# ===============================================================================

# Examine phenotype data structure
str(myPhenotype) # Check data types and structure
head(myPhenotype) # View first few rows
summary(myPhenotype) # Summary statistics

# Check for missing data
colSums(is.na(myPhenotype)) # Count NAs per column

# Visualize trait distribution (replace 'TraitName' with your actual trait column)
# hist(myPhenotype$TraitName,
#      main = "Trait Distribution",
#      xlab = "Trait Value",
#      col = "lightblue")

# Basic statistics for a trait
# mean(myPhenotype$TraitName, na.rm = TRUE)
# sd(myPhenotype$TraitName, na.rm = TRUE)
# range(myPhenotype$TraitName, na.rm = TRUE)

# ===============================================================================
# SECTION 5: BASIC GWAS ANALYSIS - GLM MODEL
# ===============================================================================

# General Linear Model (GLM)
# - Fast computation
# - Suitable for populations with minimal structure
# - Does not account for kinship/relatedness
# - Good for initial exploratory analysis

myResults_GLM <- GAPIT(
  Y = myPhenotype, # Phenotype data
  G = myGenotype, # Genotype data
  model = "GLM", # Model type
  PCA.total = 3, # Number of principal components to use
  SNP.impute = "Major", # Impute missing SNPs with major allele
  Major.allele.zero = TRUE # Code major allele as 0
)

# ===============================================================================
# SECTION 6: MIXED LINEAR MODEL (MLM) ANALYSIS
# ===============================================================================

# Mixed Linear Model (MLM)
# - Accounts for population structure (via PCs)
# - Accounts for kinship/relatedness (via kinship matrix)
# - More conservative than GLM (fewer false positives)
# - Recommended for final analysis and publication

myResults_MLM <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = "MLM", # Mixed Linear Model
  PCA.total = 3, # Number of PCs as fixed effects
  SNP.impute = "Major",
  Major.allele.zero = TRUE
)

# ===============================================================================
# SECTION 7: MULTIPLE MODEL COMPARISON
# ===============================================================================

# Run multiple models simultaneously to compare results
# Available models: GLM, MLM, SUPER, MLMM, FarmCPU, BLINK

myResults_Multi <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = c("GLM", "MLM", "FarmCPU", "BLINK"), # Multiple models
  PCA.total = 3,
  SNP.impute = "Major",
  Major.allele.zero = TRUE
)

# ===============================================================================
# SECTION 8: ADVANCED OPTIONS - KINSHIP COMPRESSION
# ===============================================================================

# Kinship compression groups similar individuals to speed up computation
# Useful for large populations (>1000 individuals)
# Default: group.by = 10 (groups of ~10 individuals)

# WITHOUT compression (each individual is its own group)
# Use when population size is small (<500) or you want maximum precision
myResults_NoCompression <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = "MLM",
  PCA.total = 3,
  SNP.impute = "Major",
  Major.allele.zero = TRUE,
  group.from = nrow(myPhenotype), # Set to population size
  group.to = nrow(myPhenotype),
  group.by = 1 # Each individual = 1 group
)

# WITH compression (default settings)
# GAPIT automatically determines optimal grouping
myResults_WithCompression <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = "MLM",
  PCA.total = 3,
  SNP.impute = "Major",
  Major.allele.zero = TRUE
  # group.by defaults to 10
)

# ===============================================================================
# SECTION 9: MODEL SELECTION AND OPTIMIZATION
# ===============================================================================

# Automatic model selection using Bayesian Information Criterion (BIC)
# - Determines optimal number of PCs for each trait
# - Tests different kinship clustering methods
# - Tests different kinship grouping strategies

myResults_Optimized <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = c("GLM", "MLM", "FarmCPU", "BLINK"),
  SNP.impute = "Major",
  PCA.total = 5, # Test up to 5 PCs
  kinship.cluster = c("complete", "ward.D"), # Clustering methods
  kinship.group = c("Mean", "Max", "Median"), # Grouping strategies
  Major.allele.zero = TRUE,
  Model.selection = TRUE # Enable BIC-based selection
)

# ===============================================================================
# SECTION 10: USING CUSTOM KINSHIP MATRIX AND COVARIATES
# ===============================================================================

# If you have pre-computed kinship matrix or custom covariates

# Load custom files
# myKinship <- read.table("kinship_matrix.txt", header = FALSE)
# myCovariates <- read.table("covariates.txt", header = TRUE)

# Run GAPIT with custom inputs
# myResults_Custom <- GAPIT(
#   Y = myPhenotype,
#   G = myGenotype,
#   KI = myKinship,        # Custom kinship matrix
#   CV = myCovariates,     # Custom covariates (e.g., environmental factors)
#   model = "MLM"
# )

# ===============================================================================
# SECTION 11: CONVERT HAPMAP TO NUMERICAL FORMAT
# ===============================================================================

# Some software requires genotype data in numerical format (0/1/2)
# GAPIT can convert HapMap format to numerical

# myConversion <- GAPIT(
#   G = myGenotype,
#   output.numerical = TRUE
# )
#
# myGD <- myConversion$GD  # Numerical genotype data
# myGM <- myConversion$GM  # SNP map information
#
# # Save numerical format
# write.table(myGD, "genotype_numerical.txt", quote = FALSE, row.names = FALSE)
# write.table(myGM, "genotype_map.txt", quote = FALSE, row.names = FALSE)

# ===============================================================================
# SECTION 12: INTERPRETING RESULTS
# ===============================================================================

# GAPIT automatically generates output files in the working directory:
#
# ASSOCIATION RESULTS:
# - GAPIT.Association.GWAS_Results.[MODEL].[TRAIT].csv
#   Contains SNP positions, p-values, effect sizes, MAF
#
# VISUALIZATIONS:
# - Manhattan plots: Show -log10(p-value) across genome
# - QQ plots: Assess model fit and population structure control
# - PCA plots: Visualize population structure
#
# GENOTYPE ANALYSIS:
# - Kinship matrix: Genetic relatedness between individuals
# - MAF distribution: Minor allele frequency
# - LD decay: Linkage disequilibrium patterns
#
# PREDICTIONS:
# - Genomic predicted breeding values for each individual
# - Prediction accuracy (correlation between observed and predicted)

# For detailed interpretation, see: docs/output_guide.md

# ===============================================================================
# SECTION 13: TIPS AND BEST PRACTICES
# ===============================================================================

# 1. Always check QQ plots first to assess model fit
# 2. Use GLM for quick exploration, MLM for final analysis
# 3. If QQ plot shows early deviation, increase PCA.total or use MLM
# 4. For large populations (>1000), use compression or FarmCPU/BLINK
# 5. Compare results across multiple models for robustness
# 6. Filter SNPs with MAF < 0.05 before analysis for better power
# 7. Ensure phenotype data is approximately normally distributed
# 8. Check for batch effects or confounding factors in your data

# ===============================================================================
# SECTION 14: TROUBLESHOOTING
# ===============================================================================

# Problem: "Matrix is singular" error
# Solution: Try different compression levels or kinship clustering methods

# Problem: No significant associations
# Solution: Check sample size, trait heritability, marker density

# Problem: Too many significant SNPs
# Solution: Use MLM instead of GLM, or increase stringency threshold

# Problem: Installation fails
# Solution: Update R to latest version, install dependencies manually

# ===============================================================================
# END OF SCRIPT
# ===============================================================================

# For more information:
# - GAPIT Manual: http://zzlab.net/GAPIT/
# - GitHub: https://github.com/jiabowang/GAPIT3
# - Publications: See README.md for citations
