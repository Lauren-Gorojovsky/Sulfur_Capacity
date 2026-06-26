# Solubility and speciation of sulfur under crustal conditions

This repository contains the supplementary R project for:

**L. R. Gorojovsky and B. J. Wood (2026), "Solubility and speciation of sulfur in silicate melts under crustal conditions", Earth and Planetary Science Letters.**

The project provides calibration data and an R workflow for fitting LASSO regression models for sulfide and sulfate capacity in silicate melts, estimating prediction intervals, and propagating uncertainty into sulfur speciation calculations. Go to page: https://lauren-gorojovsky.github.io/Solubility-and-speciation-of-sulfur-under-crustal-conditions-Supplementary-materials/S2_LASSO_regression.html 

## Repository contents

| File | Description |
| --- | --- |
| `S1_Calibration_data_and_fitting_results.xlsx` | Supplementary calibration workbook. Includes fitting data for `logCS2-` and `logCS6+`, plus fitting-result summaries. |
| `S2_LASSO_regression.R` | R script for data import, repeated cross-validation, final LASSO fits, split-conformal prediction intervals, measured-vs-predicted plots, and example Monte Carlo propagation for `S6+/Stot`. |
| `S2_LASSO_regression.qmd` | Quarto notebook version of the LASSO workflow, intended for readable, executable use in RStudio, Quarto, GitHub Codespaces, or Binder-style cloud environments. |
| `S3_Sulfur_fO2_calculator.xlsx` | Spreadsheet calculator for sulfur solubility, speciation, and `logfO2`, using major-element composition, temperature, pressure, water content, and oxygen fugacity inputs. |
| `S4_MAGEMin_Exports.xlsx` | MAGEMin model export workbook for arc and MORB compositions at different QFM offsets. |
| `Supplementary_Files.Rproj` | RStudio project file. Opening this file sets the project working directory for the analysis. |

## Workbook sheets

The calibration workbook, `S1_Calibration_data_and_fitting_results.xlsx`, contains:

| Sheet | Contents |
| --- | --- |
| `Info` | Brief description of the supplementary workbook. |
| `S1_logCS2_data_for_fitting` | Sulfide-capacity calibration data used to fit the `logCS2-` model. |
| `S2_logCS6_data_for_fitting` | Sulfate-capacity calibration data used to fit the `logCS6+` model. |
| `S3-6_Fitting_results_summary` | Summary tables for fitted model covariates and results. |

The calculator workbook, `S3_Sulfur_fO2_calculator.xlsx`, contains an `Info` sheet and calculator sheets for entering melt composition, temperature, pressure, `H2O`, and `logfO2` relative to QFM.

The MAGEMin workbook, `S4_MAGEMin_Exports.xlsx`, contains export sheets for `MAGEMin_Arc_QFM`, `MAGEMin_Arc_QFM+0.5`, `MAGEMin_Arc_QFM+1`, and `MAGEMin_MORB_QFM`.

## Requirements

The workflow was written for R and is easiest to run from RStudio. It uses the following R packages:

- `readxl`
- `glmnet`
- `tibble`
- `ggplot2`
- `dplyr`
- `patchwork`
- `ggprism`
- `scales`
- `conformalInference`

Install the CRAN packages with:

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

## Running the analysis

1. Clone or download this repository.
2. Open `Supplementary_Files.Rproj` in RStudio.
3. Confirm that the supplementary workbooks and `S2_LASSO_regression.R` are in the project root.
4. Open `S2_LASSO_regression.R`.
5. Run the script from top to bottom, or run each numbered block in order.

Most results print to the R console. In RStudio, fitted models, coefficient tables, prediction intervals, and plotting objects are also available in the Environment pane.

Alternatively, open `S2_LASSO_regression.qmd` in RStudio or another Quarto-capable editor and run the notebook chunks interactively. To render the notebook to HTML from the command line after installing Quarto, run:

```bash
quarto render S2_LASSO_regression.qmd
```

## Rendering with GitHub Actions

This repository includes a GitHub Actions workflow at `.github/workflows/quarto-render.yml`. When pushed to GitHub, the workflow installs Quarto and the required R packages, renders `S2_LASSO_regression.qmd`, and deploys the rendered HTML through GitHub Pages.

To enable it on GitHub:

1. Open the repository on GitHub.
2. Go to `Settings > Pages`.
3. Under `Build and deployment`, set `Source` to `GitHub Actions`.
4. Push changes to the `main` branch, or run `Render Quarto notebook` manually from the `Actions` tab.

After the workflow succeeds, the rendered notebook will be available from the repository's GitHub Pages site.

## Analysis overview

The script performs the following steps:

1. Imports the `logCS2-` and `logCS6+` calibration datasets from the supplementary workbook.
2. Builds design matrices for LASSO regression.
3. Runs repeated Monte Carlo cross-validation using `glmnet`.
4. Fits final sparse LASSO models at the mean `lambda.1se` value.
5. Estimates 95% prediction intervals using split-conformal prediction.
6. Recreates measured-vs-predicted plots corresponding to Fig. 3a,b of the manuscript.
7. Provides an example Monte Carlo workflow for propagating uncertainty in `logCS2-`, `logCS6+`, and `logfO2` into `S6+/Stot`.

The script sets a random seed before the repeated cross-validation and Monte Carlo sections. Small numerical differences may still occur across R, package, or platform versions.

## Reusing the workflow

To fit the models with a different calibration dataset, replace the workbook and sheet names in block (2) of `S2_LASSO_regression.R`. The replacement data should preserve the expected response variables and covariate structure used by the model matrices:

- `logCS2` for the sulfide-capacity model
- `logCS6` for the sulfate-capacity model
- `Study` and `T`, which are excluded from the model matrix in the supplied workflow

## Citation

If you use this repository, please cite the associated publication:

Gorojovsky, L. R. and Wood, B. J. (2026). Solubility and speciation of sulfur in silicate melts under crustal conditions. *Earth and Planetary Science Letters*.

## Contact

For questions about the supplementary files or workflow, contact Lauren Gorojovsky: Lauren.gorojovsky@earth.ox.ac.uk.

## License

No license file is currently included. Before making the repository public, add a license that reflects how you would like others to reuse the code and data.
