# halDRC <img src="man/figures/logo.png" align="right" height="120" />

[![GitHub](https://img.shields.io/github/last-commit/SeraphinaShi/halDRC?label=last%20update)](https://github.com/SeraphinaShi/halDRC)
[![R-CMD-check](https://github.com/SeraphinaShi/halDRC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SeraphinaShi/halDRC/actions)

> Highly Adaptive Lasso-Based Estimation of Causal Dose-Response Curves

`halDRC` is an R package for estimating **causal marginal dose-response relationships** using the **Highly Adaptive Lasso (HAL)**. It supports both binary and continuous outcomes, with optional undersmoothing and smoothness-order adaptation.

The package provides a plugin estimator that is interpretable and flexible, suitable for both simulation studies and real-world observational data.

---

## 📦 Installation

To install the development version from GitHub:

```r
# Install the devtools package if you haven't already
install.packages("devtools")

# Install halDRC
devtools::install_github("SeraphinaShi/halDRC", build_vignettes = TRUE)
```

Or clone and install locally:

```r
git clone https://github.com/SeraphinaShi/halDRC.git
setwd("halDRC")
devtools::install(build_vignettes = TRUE)
```

To install dependencies: `hal9001` and `sl3`

```r
install.packages("hal9001")
devtools::install_github("tlverse/sl3")
```


---

## 🚀 Quick Start

```r
library(halDRC)

# Simulate data
set.seed(123)

simulate_dose_response <- function(n, a = NA) {
  W <- rnorm(500)
  A <- if (is.na(a)) 0.5 * W + rnorm(n, sd = 0.8) else rep(a, n)
  logit_Y <- 0.5 * A^2 + A + 0.3 * W + rnorm(n, sd = 0.5)
  Y <- 1 / (1 + exp(-logit_Y))

  data.frame(W, A, Y)
}
data <- simulate_dose_response(n = 100)   # Replace with your generator


# Fit dose-response using HAL
library(dplyr)
fit <- fit_UHAL_DRC(
  dat = data,
  y_var_name = "Y",
  trt_var_name = "A",
  family = "gaussian"
)

# Plot results
library(ggplot2)
ggplot(fit$curve_est, aes(x = a, y = y_hat)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lwr, ymax = ci_upr), alpha = 0.3) +
  theme_minimal() +
  labs(title = "Estimated Dose-Response Curve", x = "Treatment (A)", y = "Outcome (Y)")
```

---

## 📘 Vignette Example: Real RHC Dataset

We provide a full demo using the RHC (Right Heart Catheterization) dataset from the TMLE workshop.

### 🔍 View the Vignette

```r
# After installation
browseVignettes("halDRC")
# Then click 'halDRC-demo' to view the full analysis
```

Or view the [online vignette here](https://github.com/SeraphinaShi/halDRC/blob/main/doc/halDRC-demo.html).

This vignette demonstrates the full analysis pipeline using real data:
- Loading and preprocessing the RHC dataset
- Choosing treatment and confounder variables
- Fitting the HAL plugin estimator
- Visualizing the estimated dose-response curve

---

## 📂 Included Functions

- `fit_UHAL_DRC()` – Estimate dose-response curve using HAL

---

## 📚 References

- Shi, Junming, et al. "Hal-based plugin estimation of the causal dose-response curve." arXiv preprint arXiv:2406.05607 (2024).

---

## 🧪 Development

If you want to contribute, feel free to submit an issue or pull request!

---

## License

MIT © Junming Shi
