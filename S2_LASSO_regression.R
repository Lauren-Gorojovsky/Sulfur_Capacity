

# 
#     LASSO REGRESSION & UNCERTAINTY ESTIMATES FOR LOGCS2- & LOGCS6+
#     Written by L.R. Gorojovsky, University of Oxford, October 2025
# 



# This code accompanies the publication:
# L.R. Gorojovsky & B.J. Wood, 2026, "Solubility and speciation of sulfur 
# in silicate melts under crustal conditions" 



# ─────────────────────────────────────────────────────────────────────────────
#                                   (0) 
#                         HOW TO USE THIS PROGRAM
# ─────────────────────────────────────────────────────────────────────────────

# This is a self-contained R project, designed this way to make importing the
# provided calibration data fairly straightforward for those who are 
# unfamiliar with R. 
# 
# Heres some info to get started:
# 1) Install R (Free) - from CRAN (search “CRAN download R”)
# 2) Install Rstudio desktop (this is the IDE for R - free from Posit).
# 3) Running the code should be simple enough if the following files 
#     are kept in the same folder that should be named "Supplementary_Files"
#       - S2_LASSO_regression.R
#       - Supplementary_Files.Rproj
#       - S1_Calibration_data_and_fitting_results
#     To ensure it works as intended, open the .Rproj file first.
#     In the case you would like to train the model using another data set, 
#     then this script can be used as a standalone (you dont need the .Rproj file), 
#     although you will want to modify the code in block (2), lines 88-95.
# 
# 3) The script is organised into blocks that need to be run in order, either 
#     run the whole script, or block by block. 
# 4) Most of the results will print to the console as you run the script. 
#     However, if you're using RStudio the results can also be accessed in the 
#     'Environment' tab. Any result can be exported by opening the variable 
#     from the tab and copying and pasting into a spreadsheet/text file.
# 
# Please feel free to contact Lauren for help or inquiries:
# Lauren.gorojovsky@earth.ox.ac.uk

# ─────────────────────────────────────────────────────────────────────────────
#                                   (1) 
#                         INSTALL & LOAD LIBRARIES
# ─────────────────────────────────────────────────────────────────────────────

# Uncomment the following lines to install packages if you haven't already 

# install.packages(readxl)       # read_excel()
# install.packages(glmnet)       # lasso / elastic net
# install.packages(tibble)       # tidy coefficient tibbles
# install.packages(ggplot2)      # plotting
# install.packages(dplyr)        # data manipulation
# install.packages(patchwork)    # combine ggplots
# install.packages(ggprism)      # annotation_ticks()

# The ConformalInference library is fiddly to install, I recommend doing it this way:

# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("ryantibs/conformal", subdir = "conformalInference")


# Load packagaes
library(readxl)               # read_excel()
library(glmnet)               # lasso / elastic net
library(tibble)               # tidy coefficient tibbles
library(ggplot2)              # plotting
library(dplyr)                # data manipulation
library(patchwork)            # combine ggplots
library(ggprism)              # annotation_ticks()
library(conformalInference)   # split conformal prediction



# ─────────────────────────────────────────────────────────────────────────────
#                                   (2) 
#                           DATA IMPORT & PREP
# ─────────────────────────────────────────────────────────────────────────────

# Import single-oxygen mole fractions / T to fit the form:
#               logCS = A0 + B/T + Σ(XA/T) (eq. 7)

#  import sulfate capacity df
df_S6 <- read_excel("S1_Calibration_data_and_fitting_results.xlsx", 
                   sheet = "S2_logCS6_data_for_fitting", 
                   skip = 1)

#  import sulfide capacity df 
df_S2 <- read_excel("S1_Calibration_data_and_fitting_results.xlsx", 
                   sheet = "S1_logCS2_data_for_fitting", 
                   skip = 1)


# Design matrix: no intercept, drop Study & Temp
x1 <- model.matrix(logCS6 ~ . - 1 - Study - T, df_S6)
y1 <- df_S6$logCS6

x2 <- model.matrix(logCS2 ~ . - 1 - Study - T, df_S2)
y2 <- df_S2$logCS2


# ─────────────────────────────────────────────────────────────────────────────
#                                   (3) 
#    REPEATED CV FOR LASSO (100× Monte-Carlo CV at RMSE, 20 folds each)
#    We summarise λ.1se and the distribution of selected coefficients.
# ─────────────────────────────────────────────────────────────────────────────

# fix sequence for random sampling 
set.seed(80085)

B_reps <- 100               # number of Monte-Carlo CV repeats
nfolds <- 20                # CV folds

# Helper function to run repeated cv.glmnet() and collect summaries
run_mc_cv <- function(x, y, B = 100, nfolds = 20) {
  lambda_1se_vals <- numeric(B)
  rmse_1se_vals   <- numeric(B)
  coeff_list      <- vector("list", B)
  
  for (b in seq_len(B)) {
    cv <- cv.glmnet(x, y, alpha = 1, nfolds = nfolds, type.measure = "mse")
    
    lambda_1se_vals[b] <- cv$lambda.1se
    rmse_1se_vals[b]   <- sqrt(cv$cvm[ cv$index["1se", 1] ])
    
    beta_1se <- coef(cv, s = "lambda.1se")
    coeff_list[[b]] <- tibble(
      term  = beta_1se@Dimnames[[1]][ beta_1se@i + 1 ],
      value = beta_1se@x
    )
  }
  
  lambda_1se_mean <- mean(lambda_1se_vals)
  lambda_1se_sd   <- sd(lambda_1se_vals)
  rmse_1se_mean   <- mean(rmse_1se_vals)
  
  coef_tbl <- bind_rows(coeff_list, .id = "iter") %>%
    group_by(term) %>%
    summarise(
      mean   = mean(value),
      sd     = sd(value),
      n_non0 = n(),              # how often term is nonzero across B repeats
      .groups = "drop"
    ) %>%
    mutate(select_freq = n_non0 / B)
  
  list(
    lambda_1se_mean = lambda_1se_mean,
    lambda_1se_sd   = lambda_1se_sd,
    rmse_1se_mean   = rmse_1se_mean,
    coef_tbl        = coef_tbl
  )
}

# -- Run for sulfate (S6+)
cv_S6 <- run_mc_cv(x1, y1, B = B_reps, nfolds = nfolds)
cat("\nSulfate λ-1SE (mean of 100 CVs) =", round(cv_S6$lambda_1se_mean, 5),
    "; mean RMSE =", round(cv_S6$rmse_1se_mean, 4), "\n")

# -- Run for sulfide (S2–)
cv_S2 <- run_mc_cv(x2, y2, B = B_reps, nfolds = nfolds)
cat("Sulfide λ-1SE (mean of 100 CVs) =", round(cv_S2$lambda_1se_mean, 5),
    "; mean RMSE =", round(cv_S2$rmse_1se_mean, 4), "\n\n")


# ─────────────────────────────────────────────────────────────────────────────
#                                   (4) 
#           FINAL LASSO FITS at λ-1SE (fit on ALL calibration data)
# ─────────────────────────────────────────────────────────────────────────────

fit_S6 <- glmnet(x1, y1, alpha = 1, lambda = cv_S6$lambda_1se_mean)
beta_S6 <- as.matrix(coef(fit_S6))
beta_S6_sparse <- beta_S6[ beta_S6[, 1] != 0, , drop = FALSE ]

cat("Sulfate: Sparse coefficient vector (λ-1SE on full data):\n")
print(beta_S6_sparse)

fit_S2 <- glmnet(x2, y2, alpha = 1, lambda = cv_S2$lambda_1se_mean)
beta_S2 <- as.matrix(coef(fit_S2))
beta_S2_sparse <- beta_S2[ beta_S2[, 1] != 0, , drop = FALSE ]

cat("Sulfide: Sparse coefficient vector (λ-1SE on full data):\n")
print(beta_S2_sparse)


# ─────────────────────────────────────────────────────────────────────────────
#                                   (5) 
#     95% PREDICTION INTERVALS VIA SPLIT-CONFORMAL (λ fixed from above)
#                         as shown in Fig. 3a,b
# ─────────────────────────────────────────────────────────────────────────────

# For S6+
train_S6   <- function(x, y) glmnet(x, y, alpha = 1, lambda = cv_S6$lambda_1se_mean, standardize = TRUE)
predict_S6 <- function(fit, newx) as.numeric(predict(fit, newx, s = cv_S6$lambda_1se_mean))

conf_S6 <- conformal.pred.split(
  x = x1, y = y1, x0 = x1,
  train.fun = train_S6,
  predict.fun = predict_S6,
  alpha = 0.05
)

pred_S6 <- df_S6 %>%
  mutate(pred = conf_S6$pred, lo95 = conf_S6$lo, hi95 = conf_S6$up)

# For S2–
train_S2   <- function(x, y) glmnet(x, y, alpha = 1, lambda = cv_S2$lambda_1se_mean, standardize = TRUE)
predict_S2 <- function(fit, newx) as.numeric(predict(fit, newx, s = cv_S2$lambda_1se_mean))

conf_S2 <- conformal.pred.split(
  x = x2, y = y2, x0 = x2,
  train.fun = train_S2,
  predict.fun = predict_S2,
  alpha = 0.05
)

pred_S2 <- df_S2 %>%
  mutate(pred = conf_S2$pred, lo95 = conf_S2$lo, hi95 = conf_S2$up)


# ─────────────────────────────────────────────────────────────────────────────
#                                   (6) 
#               MEASURED vs PREDICTED PLOTS Fig. 3a,b of manuscript 
#               Adds ±1.96×RMSE conformal bands around the 1:1 line
# ─────────────────────────────────────────────────────────────────────────────

# Order + shapes for Study
study_levels <- c("OM", "BW", "GW")

# Colour by experimental temperature (°C)
temp_range <- c(1000, 1500)  
common_fill <- scale_fill_viridis_c(
  option = "turbo",
  name   = "T (°C)",
  limits = temp_range,
  breaks = seq(1050, 2100, 200),
  oob    = scales::squish
)

# ---- S6+
res_S6  <- pred_S6$logCS6 - pred_S6$pred
rmse_S6 <- sqrt(mean(res_S6^2))
delta_S6 <- 1.96 * rmse_S6

pred_S6 <- pred_S6 %>%
  arrange(factor(Study, levels = study_levels)) %>%
  mutate(Study = factor(Study, levels = study_levels))

pS6 <- ggplot(pred_S6,
              aes(x = logCS6, y = pred, shape = Study, fill = T - 273)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(slope = 1, intercept =  delta_S6, linetype = "dotted") +
  geom_abline(slope = 1, intercept = -delta_S6, linetype = "dotted") +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0, linewidth = 0.2) +
  geom_point(size = 2) +
  scale_shape_manual(values = c(22, 23, 21)) +
  coord_equal() +
  common_fill +
  labs(x = "Measured logCS6+", y = "Predicted logCS6+") +
  theme_minimal() +
  coord_cartesian(clip = "off") +
  annotation_ticks(sides = "bl", type = "both", outside = TRUE) +
  scale_x_continuous(limits = c(2, 13), breaks = seq(2, 13, by = 2)) +
  scale_y_continuous(limits = c(2, 13), breaks = seq(2, 13, by = 2)) +
  theme(
    legend.key.size = unit(10,"point"),
    aspect.ratio       = .8,
    panel.grid.major.y = element_line(color = "grey", size = 0.3, linetype = 2),
    panel.background   = element_rect(colour = "black")
  )

# ---- S2–
res_S2  <- pred_S2$logCS2 - pred_S2$pred
rmse_S2 <- sqrt(mean(res_S2^2))
delta_S2 <- 1.96 * rmse_S2

pred_S2 <- pred_S2 %>%
  arrange(factor(Study, levels = study_levels)) %>%
  mutate(Study = factor(Study, levels = study_levels))

pS2 <- ggplot(pred_S2,
              aes(x = logCS2, y = pred, shape = Study, fill = T - 273)) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0, linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_abline(slope = 1, intercept =  delta_S2, linetype = "dotted") +
  geom_abline(slope = 1, intercept = -delta_S2, linetype = "dotted") +
  geom_point(size = 2) +
  scale_shape_manual(values = c(22, 23, 21)) +
  coord_cartesian(clip = "off") +
  annotation_ticks(sides = "bl", type = "both", outside = TRUE) +
  scale_x_continuous(limits = c(-8, -2), breaks = seq(-8, -2, by = 1)) +
  scale_y_continuous(limits = c(-8, -2), breaks = seq(-8, -2, by = 1)) +
  common_fill +
  labs(x = "Measured logCS2–", y = "Predicted logCS2–") +
  theme_minimal() +
  theme(
    legend.position = "none",
    aspect.ratio       = .8,
    panel.grid.major.y = element_line(color = "grey", size = 0.3, linetype = 2),
    panel.background   = element_rect(colour = "black")
  )

# Combined panel
pS2 + pS6 + plot_layout(ncol = 2)


# ─────────────────────────────────────────────────────────────────────────────
#                                   (7) 
#   MONTE-CARLO PROPAGATION OVER A LOGfO2 GRID (e.g. Figs 4 and 5)
#    Propagate uncertainty from logCS6+ and logCS2- into S6+/Stot.
# ─────────────────────────────────────────────────────────────────────────────

# The following values are just an example, modify as needed:

# Inputs (replace with your capacity+T±P values / SEs)
mu6    <- 10.824          # a value of logCS6+
mu2    <- -6.226          # a value of logCS2–
se6    <- 0.55 / 1.96     # standard error for logCS6+
se2    <- 0.35 / 1.96     # standard error for logCS2-
se_fO2 <- 0.01            # standard error for logfO2 (± 0.1 s.d., n = 100)

set.seed(80085)
R  <- 50000
draw6 <- rnorm(R, mu6, se6)
draw2 <- rnorm(R, mu2, se2)

# vector of logfO2 vals (absolute)
fO2_vec <- seq(-11, -5, by = 0.1)

out <- lapply(fO2_vec, function(lgfO2) {
  # L_draw is log10(K), where K = 10^(logCS6 - logCS2 + 2*log fO2)
  L_draw <- draw6 - draw2 + 2 * (lgfO2 + rnorm(1, 0, se_fO2))
  frac   <- 1 - 1 / (1 + 10^L_draw)  # sulfur as S6+ fraction
  tibble(
    point = 1 - 1 / (1 + 10^(mu6 - mu2 + 2 * lgfO2)),
    lo95  = quantile(frac, 0.025),
    hi95  = quantile(frac, 0.975)
  )
})

grid_df <- bind_rows(out) |>
  mutate(logfO2 = fO2_vec)

# plot S6+/S vs logfO2(FMQ) with error envelope
ggplot(grid_df, aes(x = logfO2 - (-25096.3 / 1323 + 8.735), y = point)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill="purple", alpha = .3) +
  geom_line(linewidth = 1) +
  coord_cartesian(clip = "off") +
  annotation_ticks(sides = "bl", type = "both", outside = TRUE) +
  scale_x_continuous(limits = c(-1, 4), breaks = seq(-1, 4, by = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  labs(x = expression(logfO[2]*Delta*FMQ), y = "S6+ fraction") +
  theme_minimal() +
  theme(
    aspect.ratio       = .8,
    panel.grid.major.y = element_line(color = "grey", size = 0.3, linetype = 2),
    panel.background   = element_rect(colour = "black")
  )
