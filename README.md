# Multi-omic analysis achieves high diagnostic yield and identifies accelerated biological aging in paediatric congenital heart disease and neurodevelopmental disorders


## Overview

This repository contains the computational analysis code for our multi-omic study investigating the genetic and epigenetic factors contributing to the co-occurrence of neurodevelopmental disorders (NDD) and congenital heart disease (CHD).

## Study Design

- **Cohort**: 14 trios and 1 duo (probands with NDD and/or CHD plus biological parents)
- **Data types**: 
  - Whole genome sequencing (WGS)
  - RNA sequencing (subset of 10 trios)
  - DNA methylation arrays (all 15 families)
- **Key findings**:
  - 64% diagnostic yield from genetic variants
  - Dramatic epigenetic age acceleration in CHD patients (~17 years)
  - Extreme methylation dysregulation in specific patients

## Analysis Scripts

### DNA Methylation Analysis
- **`denovo_methylation_sites.r`** - Identifies de novo methylation changes in trios using effect-size criteria (Δβ > 0.2 from both parents)
- **`epigenetic_age_analysis.r`** - Calculates epigenetic age using Hannum clock and identifies age acceleration
- **`inheritance_analysis.r`** - Analyzes parent-of-origin methylation inheritance patterns in trios
- **`methylation_pathway_analysis.r`** - Collapses CpG sites into differentially methylated regions (DMRs) and performs GO enrichment analysis

### Genetic Variant Analysis
- **`splicing_validation.py`** - Validates splice-affecting variants using RNA-seq junction analysis with inheritance awareness
- **`mae_analysis.sh`** - Monoallelic expression (MAE) analysis pipeline for identifying parent-of-origin expression bias
- **`prs_pipeline_pgs-calc.sh`** - Calculates polygenic risk scores for neurodevelopmental and psychiatric traits using pgsc_calc

## Requirements

**R packages:**
- minfi, methylclock, dplyr, ggplot2, GenomicRanges, clusterProfiler, lme4, lmerTest

**Python:**
- Python 3.6+ (standard library only)

**Command-line tools:**
- bcftools, samtools, GATK, nextflow (for shell scripts)

## Usage

Each script contains usage examples in comments. Input file formats and parameters are documented within each file.

## Citation

If you use this code, please cite:

Thompson JM, Gao Y, Iwasawa E, et al. Multi-omic analysis achieves high diagnostic yield and identifies accelerated biological aging in paediatric congenital heart disease and neurodevelopmental disorders. [Journal] 2025; [details to be added upon publication]

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
- - Cincinnati Children's Hospital Medical Center
- The families who participated in this study
