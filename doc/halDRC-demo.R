## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(ggplot2)
library(hal9001)
library(glmnet)
library(here)
library(halDRC)

## -----------------------------------------------------------------------------
set.seed(145)

generate_data <- function(n, a = NA){
  U_W <- rnorm(n)
  U_A <- rnorm(n, sd = 2)
  U_Y <- runif(n)
  W <- U_W

  if (is.na(a)) {
    A <- 2 - 0.5 * W + U_A
    A <- pmin(pmax(A, 0), 5)
  } else {
    A <- rep(a, n)
  }

  Y <- as.numeric(U_Y < plogis(-3 + 0.5*W + 1.25*A - 0.5*W*A))
  data.frame(W, A, Y)
}

obs <- generate_data(n = 500)
head(obs)

## -----------------------------------------------------------------------------
DRC_fit_UAdaptive <- fit_UHAL_DRC(dat = obs, y_var_name = "Y", trt_var_name = "A", family = "binomial")

## -----------------------------------------------------------------------------
df_ests <- DRC_fit_UAdaptive$curve_est

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ci_lwr, ymax = ci_upr), alpha = 0.3) +
  labs(title = "Estimated Dose-Response (Binary Outcome)",
       x = "Treatment (A)",
       y = "E[Y | do(A = a)]") +
  theme_bw()

## -----------------------------------------------------------------------------
a.vals <- seq(0, 5, length.out = 20)
psi0_a <- sapply(a.vals, function(a) mean(generate_data(n = 1e6, a = a)$Y))
psi0 <- data.frame(a = a.vals, psi0 = psi0_a)

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lwr, ymax = ci_upr), alpha = 0.3) +
  geom_point(data = psi0, aes(x = a, y = psi0), color = "red") +
  geom_line(data = psi0, aes(x = a, y = psi0), color = "red") +
  labs(title = "Estimated vs. True Dose-Response (Binary Y)",
       y = "E[Y | do(A = a)]") +
  theme_bw()

## -----------------------------------------------------------------------------
generate_data_2 <- function(n, a = NA) {
  U_W <- rnorm(n)
  U_A <- rnorm(n, sd = 0.8)
  U_Y <- rnorm(n, sd = 0.5)

  W <- U_W
  A <- if (is.na(a)) 0.5 * W + U_A else rep(a, n)
  logit_Y <- 0.5 * A^2 + A + 0.3 * W + U_Y
  Y <- 1 / (1 + exp(-logit_Y))
  data.frame(W, A, Y)
}

obs <- generate_data_2(n = 500)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2))
plot(obs$W, obs$A, main = "W vs A", pch = 19, col = rgb(0, 0, 1, 0.4))
plot(obs$A, obs$Y, main = "A vs Y", pch = 19, col = rgb(1, 0, 0, 0.4))

## -----------------------------------------------------------------------------
DRC_fit_UAdaptive <- fit_UHAL_DRC(obs, "Y", "A", "gaussian")
df_ests <- DRC_fit_UAdaptive$curve_est

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ci_lwr, ymax = ci_upr), alpha = 0.3) +
  labs(title = "Estimated Dose-Response (Continuous Outcome)",
       x = "Treatment (A)",
       y = "E[Y | do(A = a)]") +
  theme_bw()

## -----------------------------------------------------------------------------
a.vals <- seq(min(obs$A), max(obs$A), length.out = 100)
psi0 <- data.frame(
  a = a.vals,
  psi0 = sapply(a.vals, function(a) mean(generate_data_2(n = 1e6, a = a)$Y))
)

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lwr, ymax = ci_upr), alpha = 0.3) +
  geom_line(data = psi0, aes(x = a, y = psi0), color = "red") +
  labs(title = "Estimated vs. True Dose-Response (Continuous Y)",
       y = "E[Y | do(A = a)]") +
  theme_bw()

## ----get_real_data------------------------------------------------------------
# Load raw dataset
ObsData <- read.csv("https://hbiostat.org/data/repo/rhc.csv", header = TRUE)

# Create outcome Y: length of stay
ObsData$Y <- ObsData$dschdte - ObsData$sadmdte
ObsData$Y[is.na(ObsData$Y)] <- ObsData$dthdte[is.na(ObsData$Y)] - ObsData$sadmdte[is.na(ObsData$Y)]

# Add synthetic ProcedureCount (used as A)
set.seed(145)
ObsData$ProcedureCount <- round(runif(nrow(ObsData), min = 0, max = 10), 0)

# Rename variables needed for analysis
ObsData <- ObsData %>%
  dplyr::rename(
    APACHE.score = aps1,
    DASIndex = das2d3pc,
    Creatinine = crea1,
    Bilirubin = bili1,
    WBC = wblc1,
    Heart.rate = hrt1
  )

# Subset to required columns for analysis
analysis_df <- ObsData %>%
  dplyr::select(
    Y,
    A = ProcedureCount,
    APACHE.score,
    DASIndex,
    Creatinine,
    Bilirubin,
    WBC,
    Heart.rate
  ) %>%
  na.omit()

# ℹ️ For computational efficiency in this demo, we randomly selected 1,000 observations from the full cleaned dataset for this HAL analysis.
set.seed(2025)
analysis_df <- analysis_df %>% sample_n(1000)

head(analysis_df)

## ----halDRC_realdata----------------------------------------------------------
DRC_fit_UAdaptive <- fit_UHAL_DRC(dat = analysis_df, y_var_name = "Y", trt_var_name = "A", family = "gaussian")
df_ests <- DRC_fit_UAdaptive$curve_est

head(df_ests)

ggplot(df_ests, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = ci_lwr, ymax = ci_upr), alpha = 0.3) +
  labs(
    title = "Estimated Marginal Dose-Response Curve",
    x = "Number of Procedures (A)",
    y = "Expected Length of Stay (Y)"
  ) +
  theme_bw()

