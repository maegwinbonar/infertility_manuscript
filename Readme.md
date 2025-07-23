# Infertility_manuscript

This repository contains pipelines for data processing and analysis for the manuscript titled "Genome-wide analysis reveals loci associated with infertility under positive selection in human populations".

## Project Overview

We compiled a list of 31 infertility-associated loci (IAL) (28 loci associated with female infertility and 3 loci associated with male infertility) identified by large GWAS studies and tested for signatures of positive selection and archaic introgression using 1073 high coverage whole-genome sequencing data from 69 geographically diverse human populations including 177 newly-generated genomes from 23 populations in Near Oceania. We used two complementary test statistics for positive selection: PBS (population branch statistic) and XP-EHH (cross-population extended haplotype homozygosity) and reported evidence of selection using Fisher’s combined scores.

## Repository Structure

This repository was made for reproducibility purposes, the full script for each analysis is provided in an .Rmd file and corresponding html.

Description of the scripts in order of execution:

- **01_constraint** - Tests for genetic constraint on 25 genes containing IAL
- **02_ai_deserts** - Overlap IAL with archaic deserts from Sankararaman et al. 2014
- **selection scans** - PBS and XP-EHH methods
- **03_positive_selection** - analysis of PBS and XP-EHH, generating heatmaps
- **04_allele_frequency** - Allele frequencies of risk alleles
- **05_ld_calculation** - Calculating r2 for IAL

- **scripts** - Additional R Scripts and notes

## Contact

- **Author**: Maegwin Bonar
- **Email**: maegwinbonar@gmail.com
- **Institution**: Yale University

---

**Note**: This repository is part of ongoing research. Please contact the authors before using the data or code for publication purposes.
