#### Script to generate final, cleaned Excel files of regional estimates for ####
#### domain annotations and exon annotations ####

library(tidyverse)
library(writexl)

## Domain results
# Read in domain data
domains <- read_tsv("../data/domain_variant_counts.txt", show_col_types = FALSE)

# Read in MCMC summary info
results <- read.csv("../data/PRIME_domain_SummaryInfo.csv", row.names = 1)

# Filter out genes that have no variants
genes_no_variants <- domains %>%
  group_by(gene_id) %>%
  summarise(sum_total_variants = sum(total_variants)) %>%
  filter(sum_total_variants == 0)
domains <- domains %>%
  filter(! gene_id %in% genes_no_variants$gene_id)

# Extract theta posterior modes and highest density intervals
theta_modes <- results[grep("theta", rownames(results)), ]$Mode
theta_lower <- results[grep("theta", rownames(results)), ]$HDIlow
theta_upper <- results[grep("theta", rownames(results)), ]$HDIhigh

# Add theta modes and highest density intervals to domain data
domains$theta <- theta_modes
domains$theta_lower <- theta_lower
domains$theta_upper <- theta_upper

# Parameter names in MCMC
domains <- domains %>% 
  mutate(mcmc_id = paste0("theta", "[", as.character(row_number()), "]"))

# Extract genes
genes <- domains %>% 
  dplyr::select(gene_id, gene_name) %>% 
  unique() %>% 
  arrange(gene_id)

# Extract omega posterior modes, highest density intervals, and overall omega
omega_modes <- results[grep("omega", rownames(results)), ]$Mode
omega_modes <- omega_modes[-length(omega_modes)]  # don't include overall omega

omega_lower <- results[grep("omega", rownames(results)), ]$HDIlow
omega_lower <- omega_lower[-length(omega_lower)]

omega_upper <- results[grep("omega", rownames(results)), ]$HDIhigh
omega_upper <- omega_upper[-length(omega_upper)]

overall_omega_mode <- results[grep("omegaO", rownames(results)), ]$Mode

# Extract kappa posterior modes and highest density intervals
kappa_modes <- results[grep("kappa", rownames(results)), ]$Mode

kappa_lower <- results[grep("kappa", rownames(results)), ]$HDIlow

kappa_upper <- results[grep("kappa", rownames(results)), ]$HDIhigh

# Add gene-level information
genes$omega <- omega_modes
genes$omega_lower <- omega_lower
genes$omega_upper <- omega_upper

genes$kappa <- kappa_modes
genes$kappa_lower <- kappa_lower
genes$kappa_upper <- kappa_upper

# Parameter names in MCMC
genes <- genes %>% 
  mutate(mcmc_id_omega = paste0("omega", "[", as.character(row_number()), "]"),
         mcmc_id_kappa = paste0("kappa", "[", as.character(row_number()), "]"))

## Exon results
# Read in exon data
exons <- read_tsv("../data/exon_variant_counts.txt", show_col_types = FALSE)

# Read in MCMC summary info
results <- read.csv("../data/PRIME_exon_SummaryInfo.csv", row.names = 1)

# Filter out genes that have no variants
genes_no_variants <- exons %>%
  group_by(gene_id) %>%
  summarise(sum_total_variants = sum(total_variants)) %>%
  filter(sum_total_variants == 0)
exons <- exons %>%
  filter(! gene_id %in% genes_no_variants$gene_id)

# Extract theta posterior modes and highest density intervals
theta_modes <- results[grep("theta", rownames(results)), ]$Mode
theta_lower <- results[grep("theta", rownames(results)), ]$HDIlow
theta_upper <- results[grep("theta", rownames(results)), ]$HDIhigh

# Add theta modes and highest density intervals to exon data
exons$theta <- theta_modes
exons$theta_lower <- theta_lower
exons$theta_upper <- theta_upper

# Parameter names in MCMC
exons <- exons %>% 
  mutate(mcmc_id = paste0("theta", "[", as.character(row_number()), "]"))

# Extract genes
genes <- exons %>% 
  dplyr::select(gene_id, gene_name) %>% 
  unique() %>% 
  arrange(gene_id)

# Extract omega posterior modes, highest density intervals, and overall omega
omega_modes <- results[grep("omega", rownames(results)), ]$Mode
omega_modes <- omega_modes[-length(omega_modes)]  # don't include overall omega

omega_lower <- results[grep("omega", rownames(results)), ]$HDIlow
omega_lower <- omega_lower[-length(omega_lower)]

omega_upper <- results[grep("omega", rownames(results)), ]$HDIhigh
omega_upper <- omega_upper[-length(omega_upper)]

overall_omega_mode <- results[grep("omegaO", rownames(results)), ]$Mode

# Extract kappa posterior modes and highest density intervals
kappa_modes <- results[grep("kappa", rownames(results)), ]$Mode

kappa_lower <- results[grep("kappa", rownames(results)), ]$HDIlow

kappa_upper <- results[grep("kappa", rownames(results)), ]$HDIhigh

# Add gene-level information
genes$omega <- omega_modes
genes$omega_lower <- omega_lower
genes$omega_upper <- omega_upper

genes$kappa <- kappa_modes
genes$kappa_lower <- kappa_lower
genes$kappa_upper <- kappa_upper

# Parameter names in MCMC
genes <- genes %>% 
  mutate(mcmc_id_omega = paste0("omega", "[", as.character(row_number()), "]"),
         mcmc_id_kappa = paste0("kappa", "[", as.character(row_number()), "]"))

## Clean files
cleaned_file_domains <- domains %>% select(chromosome, strand, min_start, max_end, coding_region, exon_number, amino_start, amino_end, gene_name, gene_id, transcript_id, gene_region_identifier, homologous_superfamily, missense, lof, synonymous, total_variants, theta, theta_lower, theta_upper, pathogenic_counts) %>% 
  rename(theta_HDI_low = theta_lower, theta_HDI_high = theta_upper)

write_xlsx(cleaned_file_domains, "../data/PRIME_domain_estimates.xlsx")

cleaned_file_exons <- exons %>% select(chromosome, strand, start, end, exon_number, gene_name, gene_id, transcript_id, missense, lof, synonymous, total_variants, theta, theta_lower, theta_upper, pathogenic_counts) %>% 
  rename(theta_HDI_low = theta_lower, theta_HDI_high = theta_upper)

write_xlsx(cleaned_file_exons, "../data/PRIME_exon_estimates.xlsx")