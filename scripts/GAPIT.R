################################################################################
# Basic GWAS Analysis with GAPIT3
#
# This is a simplified script for performing GWAS analysis using GAPIT3.
# For more advanced options, see GAPIT3_IIRR.R
#
# Author: [Your Name]
# Date: January 2026
################################################################################

# ===============================================================================
# 1. INSTALL AND LOAD GAPIT3
# ===============================================================================

# Install GAPIT3 from GitHub (run once)
if (!require("remotes")) install.packages("remotes")
remotes::install_github("jiabowang/GAPIT3", force = TRUE)

# Load the library
library(GAPIT3)

# ===============================================================================
# 2. SET WORKING DIRECTORY
# ===============================================================================

# Set to your data directory
# IMPORTANT: Update this path to match your actual data location
setwd("path/to/your/data")

# Verify the directory
getwd()

# ===============================================================================
# 3. LOAD DATA FILES
# ===============================================================================

# Load phenotype data
# Format: Tab-delimited file with header
# Column 1: Taxa/Sample names
# Column 2+: Trait values (numeric)
myPhenotype <- read.table("phenotype.txt", header = TRUE)

# Load genotype data (HapMap format)
# Format: Tab-delimited file without header
# First 11 columns: SNP information
# Remaining columns: Genotype calls
myGenotype <- read.table("genotype.hmp.txt", header = FALSE)

# ===============================================================================
# 4. EXAMINE DATA
# ===============================================================================

# Check phenotype structure
str(myPhenotype)
head(myPhenotype)

# Check for missing values
colSums(is.na(myPhenotype))

# ===============================================================================
# 5. RUN GWAS ANALYSIS
# ===============================================================================

# OPTION A: Single Model Analysis (GLM - Fast)
# Good for initial exploration
myResults_GLM <- GAPIT(
  Y = myPhenotype, # Phenotype data
  G = myGenotype, # Genotype data
  model = "GLM", # General Linear Model
  PCA.total = 3 # Use 3 principal components
)

# OPTION B: Single Model Analysis (MLM - Recommended)
# Accounts for population structure and kinship
# More conservative, fewer false positives
myResults_MLM <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = "MLM", # Mixed Linear Model
  PCA.total = 3
)

# OPTION C: Multiple Model Comparison
# Run several models and compare results
# Models: GLM, MLM, SUPER, MLMM, FarmCPU, BLINK
myResults_Multi <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = c("GLM", "MLM", "FarmCPU", "BLINK"),
  PCA.total = 3
)

# ===============================================================================
# 6. RESULTS INTERPRETATION
# ===============================================================================

# GAPIT automatically generates output files:
#
# Key Output Files:
# - GAPIT.Association.GWAS_Results.[MODEL].[TRAIT].csv
#   → Main results: SNP positions, p-values, effect sizes
#
# - Manhattan plots (PDF)
#   → Visual representation of associations across genome
#   → Peaks indicate significant associations
#
# - QQ plots (PDF)
#   → Check model fit and population structure control
#   → Points should follow diagonal line
#
# - PCA results
#   → Population structure visualization
#
# For detailed interpretation, see: ../docs/output_guide.md

# ===============================================================================
# 7. ADDITIONAL PARAMETERS (OPTIONAL)
# ===============================================================================

# Advanced analysis with more options
myResults_Advanced <- GAPIT(
  Y = myPhenotype,
  G = myGenotype,
  model = "MLM",
  PCA.total = 5, # Test more PCs
  SNP.impute = "Major", # Impute missing SNPs
  Major.allele.zero = TRUE, # Code major allele as 0
  Model.selection = TRUE, # Auto-select optimal # of PCs
  kinship.cluster = "ward.D", # Kinship clustering method
  kinship.group = "Mean" # Kinship grouping strategy
)

# ===============================================================================
# QUICK REFERENCE
# ===============================================================================

# Common GAPIT Parameters:
# - Y: Phenotype data frame
# - G: Genotype data (HapMap, VCF, or numeric)
# - model: "GLM", "MLM", "SUPER", "MLMM", "FarmCPU", "BLINK"
# - PCA.total: Number of principal components (typically 3-5)
# - SNP.impute: "Major", "Minor", or "Middle"
# - Major.allele.zero: TRUE or FALSE
# - Model.selection: TRUE (auto-select PCs) or FALSE
# - KI: Custom kinship matrix
# - CV: Custom covariates

# For help:
# ?GAPIT
# help(package = "GAPIT3")

################################################################################
# END OF SCRIPT
################################################################################
