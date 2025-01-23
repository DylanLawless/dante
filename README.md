# Heracles for genetic report generation

This repository contains R code for generating clinical genetic reports from the output of a Whole Genome Sequencing (WGS) pipeline. The main functionality transforms raw pipeline data into structured reports (PDF/HTML) that provide variant interpretation and annotations for clinical analysis.

<img src="images/logo_heracles.webp" style="width: 80%;" alt="Heracles logo"/>

## Overview

### Key Features
- **Automated report generation**: Generates PDF reports based on patient-specific data.
- **Customisable input parameters**: Allows parameterised execution for single-case analysis.
- **Data visualisation**: Includes PCA plots and tabular summaries of variants.
- **Standards compliant**: Adheres to ACMG guidelines for variant interpretation.

## Requirements
### Software
- R (>= 4.0)
- LaTeX (e.g., TeX Live, MikTeX) with `xelatex`

### R Packages
- `rmarkdown`
- `kableExtra`
- `dplyr`
- `digest`
- `stringr`
- `ggplot2`

Install packages using:
```R
install.packages(c("rmarkdown", "kableExtra", "dplyr", "digest", "stringr", "ggplot2"))
```

## Usage
### Input Data
Ensure the following input files are available:
- WGS pipeline output files from ACMGuru/GuRu (i.e. .Rds).
- Phenotypic and clinical metadata.
- Supporting annotation files (e.g., population reference data).

### Execution Steps
1. **Prepare input data**: Place the required data files in the appropriate directories.
2. **Set parameters**: Update parameters in the `report.Rmd` file or pass them dynamically via the R script.
3. **Run Report Generation**:
   ```R
   source("run_reports.R")
   ```
   This script processes the input data and generates reports for each case.

### Output
- Individual PDF reports, named by sample ID and rank (e.g., `report_priority_1_sample_ABC123.pdf`).

## References
- [Pipedev documentation](https://swisspedhealth-pipelinedev.github.io/docs/)

## Support
For issues, please create a GitHub issue or contact the repository maintainer.


