# Multiomic investigation of shared genetic pathways in paediatric congenital heart disease and neurodevelopmental disorders


## Overview

This repository contains the computational analysis code for our multiomic study investigating the genetic and epigenetic factors contributing to the co-occurrence of neurodevelopmental disorders (NDD) and congenital heart disease (CHD).

## Study Design

- **Cohort**: 14 trios and 1 duo (probands with NDD and/or CHD plus biological parents)
- **Data types**: 
  - Whole genome sequencing (WGS)
  - RNA sequencing (subset of 10 trios)
  - DNA methylation arrays (15 trios)
- **Key findings**:
  - 57% P/LP detection rate from genetic variants
  - Epigenetic age acceleration in CHD patients
  - Methylation dysregulation in specific patients

## Analysis Scripts

### DNA Methylation Analysis
- **`dmr_analysis.r`** - Identifies differentially methylated probes in study probands compared to control samples from GEO
- **`epigenetic_age_analysis.r`** - Calculates epigenetic age using Hannum clock and identifies age acceleration
- **`methylation_pathway_analysis.r`** - Collapses CpG sites into differentially methylated regions (DMRs) and performs GO enrichment analysis

### Genetic Variant Analysis
- **`mae_analysis.sh`** - Monoallelic expression analysis pipeline for identifying parent-of-origin expression bias
- **`prs_pipeline_pgs-calc.sh`** - Calculates polygenic risk scores for neurodevelopmental and psychiatric traits using pgsc_calc

## Citation

If you use this code, please cite:

Thompson, J. M., Gao, Y., Iwasawa, E., Das, D., Rath, E., Troup, M., Humphreys, D. T., Heydarian, H., Anixt, J., Kasparian, N. A., Froehlich, T. E., Tchieu, J., Weaver, K. N., Congenital Heart Disease Synergy Study Group, Kirk, E. P., Dale, R., Dunwoodie, S. L., Winlaw, D. S., & Giannoulatou, E. (2026). Multiomic investigation of shared genetic pathways in paediatric congenital heart disease and neurodevelopmental disorders. Human Mutation, 2026, 7869246. https://doi.org/10.1155/humu/7869246

## License

This project is licensed under the CC License - see the LICENSE file for details.

## Contact

For questions about the code or analysis:
- Jamie-Lee Thompson: j.thompson@victorchang.edu.au

For questions about the study or data access:
- David Winlaw: david.winlaw@northwestern.edu
- Eleni Giannoulatou: e.giannoulatou@victorchang.edu.au

## Acknowledgments

- Victor Chang Cardiac Research Institute
- Cincinnati Children's Hospital Medical Center
- The families who participated in this study
