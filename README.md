# PRIME

PRIME (Probabilistic Regional Intolerance to Missense Estimation) is a Bayesian hierarchical modeling approach for measuring intolerance, or constraint, at sub-genic resolution across coding regions of the genome. The model estimates the probability that a variant observed in the population is missense within a genomic region using gnomAD v4 exome and genome sequencing data. Regions under stronger purifying selection are expected to have lower probabilities, reflecting a depletion of functional variants in the population.

Through its beta–binomial hierarchical structure, PRIME induces shrinkage on sub-genic estimates of intolerance by allowing them to be jointly informed by other regions within the same gene and by genes across the genome. This borrowing of information stabilizes estimates in regions with sparse variant counts. The Bayesian framework also yields full posterior distributions for each region, enabling principled quantification of uncertainty in the intolerance estimates.

---

## Repository contents

Final processed summary files are provided in the `data/` directory. Two sets of estimates are included:

- **Domain-based regions** (defined by CATH/Gene3D boundaries):
  `PRIME_domain_estimates.xlsx`

- **Exon-based regions** (defined by Ensembl exon boundaries):
  `PRIME_exon_estimates.xlsx`

These files contain posterior summary statistics for all modeled regions.

---

## Reproducing the analysis

To reproduce the PRIME analysis, first clone the repository:

git clone https://github.com/Costa-Stavrianidis/PRIME.git

cd PRIME

The analysis is performed in R and requires JAGS to be installed. The pipeline should be run in the following order:

1. Fit the PRIME model

   Rscript code/fit_PRIME_model.R

   This script fits the Bayesian hierarchical model using JAGS and generates raw MCMC output along with summary statistics for genome, gene, and region-level intolerance.

2. Generate final summary files

   Rscript code/generate_final_files.R

   This script post-processes the MCMC output and produces the cleaned .xlsx files containing posterior summary statistics for domain and exon-defined regions.
