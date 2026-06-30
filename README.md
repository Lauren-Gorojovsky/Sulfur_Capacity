# Solubility and speciation of sulfur under crustal conditions

This repository contains the supplementary materials for:

**L. R. Gorojovsky and B. J. Wood (2026), "Solubility and speciation of sulfur in silicate melts under crustal conditions", Earth and Planetary Science Letters.**

The files provide the calibration data, model-fitting workflow, sulfur speciation calculator, and MAGEMin exports used alongside the manuscript. The R workflow fits LASSO regression models for sulfide and sulfate capacity in silicate melts, estimates prediction intervals, and propagates uncertainty into sulfur speciation calculations.

A rendered version of the LASSO notebook is available here:

[View the Quarto notebook](https://lauren-gorojovsky.github.io/Solubility-and-speciation-of-sulfur-under-crustal-conditions-Supplementary-materials/S2_LASSO_regression.html)

## Files

| File | Description |
| --- | --- |
| `S1_Calibration_data_and_fitting_results.xlsx` | Calibration data for the `logCS2-` and `logCS6+` models, together with fitting-result summaries. |
| `S2_LASSO_regression.R` | R script used to fit the LASSO models, estimate prediction intervals, and reproduce the measured-vs-predicted plots. |
| `S2_LASSO_regression.qmd` | Quarto notebook version of the LASSO workflow. This is intended as a more readable, step-by-step version of the R script. |
| `S3_Sulfur_fO2_calculator.xlsx` | Spreadsheet calculator for sulfur solubility, sulfur speciation, and oxygen fugacity calculations. |
| `S4_MAGEMin_Exports.xlsx` | MAGEMin model outputs for arc and MORB compositions at different QFM offsets. |
| `Supplementary_Files.Rproj` | RStudio project file for running the R workflow locally. |

## Running the R workflow

The easiest way to run the analysis locally is to open `Supplementary_Files.Rproj` in RStudio and then run `S2_LASSO_regression.R` from top to bottom. The numbered sections of the script can also be run one at a time.

The same workflow is provided as a Quarto notebook in `S2_LASSO_regression.qmd`. This notebook can be opened in RStudio or any Quarto-capable editor and run interactively.

The workflow uses the following R packages:

- `readxl`
- `glmnet`
- `tibble`
- `ggplot2`
- `dplyr`
- `patchwork`
- `ggprism`
- `scales`
- `conformalInference`

Most packages can be installed from CRAN:

```r
install.packages(c(
  "readxl",
  "glmnet",
  "tibble",
  "ggplot2",
  "dplyr",
  "patchwork",
  "ggprism",
  "scales"
))
```

The `conformalInference` package can be installed from GitHub:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("ryantibs/conformal", subdir = "conformalInference")
```

## Analysis overview

The LASSO workflow:

1. Imports the sulfide and sulfate capacity calibration datasets from `S1_Calibration_data_and_fitting_results.xlsx`.
2. Builds model matrices for `logCS2-` and `logCS6+`.
3. Runs repeated Monte Carlo cross-validation using `glmnet`.
4. Fits final sparse LASSO models using the mean `lambda.1se` value.
5. Estimates 95% prediction intervals using split-conformal prediction.
6. Recreates the measured-vs-predicted plots shown in the manuscript.
7. Provides an example Monte Carlo calculation for propagating uncertainty into `S6+/Stot`.

The script sets a random seed before the repeated cross-validation and Monte Carlo sections. Small numerical differences may still occur between R or package versions.

## Reusing the workflow

The script can be adapted to fit the models with another calibration dataset. To do this, replace the workbook and sheet names in block (2) of `S2_LASSO_regression.R`.

The replacement data should contain the same response variables and model structure used in the supplied calibration files:

- `logCS2` for the sulfide capacity model
- `logCS6` for the sulfate capacity model
- `Study` and `T`, which are included in the data file but excluded from the model matrix

## Citation

If you use these materials, please cite:

Gorojovsky, L. R. and Wood, B. J. (2026). Solubility and speciation of sulfur in silicate melts under crustal conditions. *Earth and Planetary Science Letters*.

## Contact

For questions about the supplementary files or workflow, please contact Lauren Gorojovsky: Lauren.gorojovsky@earth.ox.ac.uk.
