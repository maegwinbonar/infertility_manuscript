# Infertility_manuscript

This repository contains pipelines for data processing and analysis for the manuscript titled "Genome-wide analysis reveals loci associated with infertility under positive selection in human populations".

## Project Overview

We compiled a list of 31 infertility-associated loci (IAL) (28 loci associated with female infertility and 3 loci associated with male infertility) identified by large GWAS studies and tested for signatures of positive selection and archaic introgression using 1073 high coverage whole-genome sequencing data from 69 geographically diverse human populations including 177 newly-generated genomes from 23 populations in Near Oceania. We used two complementary test statistics for positive selection: PBS (population branch statistic) and XP-EHH (cross-population extended haplotype homozygosity) and reported evidence of selection using Fisher’s combined scores.

## Repository Structure

The repository is organized into the following sections:

- **genetic constraint analysis** - Tests for genetic constraint on 25 genes containing IAL
- **overlap with archaic deserts** - Overlap IAL with archaic deserts from Sankararaman et al. 2014
- **selection scans** - PBS and XP-EHH methods
- **linkage disequilibrium** - Calculating r2 for IAL
- **allele frequency** - Allele frequencies of risk alleles
- **documentation** - Additional documentation and notes

### Prerequisites
```
- R (version 4.0+)
- Python (version 3.8+)
```

## Documentation

Please see the README.md file in each folder for detailed information about the corresponding analysis:

- [Genetic Constraint](./src/genetic_constraint/README.md)
- [Archaic Deserts](./src/ai_deserts/README.md)
- [Selection Scans](./src/selection_scans/README.md)
- [Linkage Disequilibrium](./linkage_disequilibrium/README.md)
- [Allele Frequency](./allele_frequency/README.md)

## Contact

- **Author**: Maegwin Bonar
- **Email**: maegwin.bonar@yale.edu
- **Institution**: Yale University

---

**Note**: This repository is part of ongoing research. Please contact the authors before using the data or code for publication purposes.
