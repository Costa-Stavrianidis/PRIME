#### Script to fit PRIME model to gnomAD v4 variant data aggregated by CATH/Gene3D ####
#### domain annotations and exon annotations ####

library(tidyverse)
library(rjags)

set.seed(100)

### Domain model
## Preprocess the variant data
# Read in variant data
variant_counts <- read_tsv("../data/domain_variant_counts.txt")

# Filter out genes that have no variants
genes_no_variants <- variant_counts %>%
  group_by(gene_id) %>%
  summarise(sum_total_variants = sum(total_variants)) %>%
  filter(sum_total_variants == 0)
variant_counts <- variant_counts %>%
  filter(! gene_id %in% genes_no_variants$gene_id)

# Factorize gene ID column
variant_counts <- variant_counts %>% 
  mutate(gene_id = as.factor(gene_id))

## Run MCMC
# Import Bayesian hierarchical model
source("PRIME_model_wrapper.R")

# Generate the MCMC chain and save
mcmcCoda = genMCMC(data = variant_counts,
                   xName = "missense", NName = "total_variants", sName = "subregion",
                   gName = "gene_id", numSavedSteps = 20000, thinSteps = 1, 
                   saveName = "../data/PRIME_domain_")

# Generate summary info and save
summaryInfo <- smryMCMC(mcmcCoda, compVal = NULL, saveName = "../data/PRIME_domain_")

### Exon model
## Preprocess the variant data
# Read in variant data
variant_counts <- read_tsv("../data/exon_variant_counts.txt")

# Filter out genes that have no variants
genes_no_variants <- variant_counts %>%
  group_by(gene_id) %>%
  summarise(sum_total_variants = sum(total_variants)) %>%
  filter(sum_total_variants == 0)
variant_counts <- variant_counts %>%
  filter(! gene_id %in% genes_no_variants$gene_id)

# Factorize gene ID column
variant_counts <- variant_counts %>% 
  mutate(gene_id = as.factor(gene_id))

## Run MCMC
# Generate the MCMC chain and save
mcmcCoda = genMCMC(data = myData,
                   xName = "missense", NName = "total_variants", sName = "subregion",
                   gName = "gene_id", numSavedSteps = 20000, thinSteps = 1, 
                   saveName = "../data/PRIME_exon_")

# Generate summary info and save
summaryInfo <- smryMCMC(mcmcCoda, compVal = NULL, saveName = "../data/PRIME_exon_")