% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_DRC.R
\name{fit_UHAL_DRC}
\alias{fit_UHAL_DRC}
\title{Fit HAL and Estimate the Dose-Response Curve with Confidence Intervals}
\usage{
fit_UHAL_DRC(
  dat,
  y_var_name,
  trt_var_name,
  family,
  curvePoints = NULL,
  Unsersmoothing = TRUE,
  smoothOrderAdapt = FALSE,
  smoothOrder = 1,
  baseNumKnots = NA,
  boundResults = TRUE
)
}
\arguments{
\item{y_var_name}{A \code{character} string specifying the column name of the
outcome variable in the data.}

\item{trt_var_name}{A \code{character} string specifying the column name of the
treatment (exposure) variable in the data.}

\item{family}{A \code{character} string or a \code{\link[stats]{family}} object
indicating the error distribution and link function to be used in the model.
Must be compatible with the HAL fitting procedure.}

\item{curvePoints}{A \code{numeric vector} containing the treatment values at
which the dose-response curve will be evaluated.}

\item{smoothOrderAdapt}{A \code{logical} flag. If \code{TRUE}, the function will
attempt to data-adaptively choose the smoothness order and number of knots.}

\item{smoothOrder}{A \code{numeric} value specifying the smoothness order for the HAL
basis functions (used only if \code{smoothOrderAdapt = FALSE}).}

\item{baseNumKnots}{A \code{numeric} value indicating the base number of knots to use
in HAL fitting (used only if \code{smoothOrderAdapt = FALSE}).}

\item{boundResults}{A \code{logical} flag. If \code{TRUE}, the function will constrain
the predicted dose-response estimates and confidence intervals to fall within the
observed range of the outcome variable.}

\item{data}{A \code{data.frame} containing the observations. Must include
all relevant covariates, the treatment variable, and the outcome variable.}
}
\value{
A list or data.frame (depending on implementation) containing:
\itemize{
  \item The predicted dose-response curve evaluated at \code{curvePoints}.
  \item Optional confidence intervals for each point on the curve.
  \item Possibly other diagnostic or metadata components used during HAL fitting.
}
}
\description{
This function fits an undersmoothed Highly Adaptive Lasso (HAL) model to the
input data using the specified fitting hyperparameters. It then leverages
the fitted HAL model to estimate the dose-response curve across specified
treatment levels, optionally returning confidence intervals via the delta method.
}
\details{
The core steps include:
1. Fitting a HAL model using the outcome and covariates, excluding the treatment.
2. Predicting counterfactual outcomes for each observation across a range of
   treatment values supplied by `curvePoints`.
3. Estimating the average potential outcome at each treatment level to construct
   the dose-response curve.
4. Optionally computing delta-method-based confidence intervals for the estimated
   means.
5. Optionally bounding the predicted means and confidence intervals within the
   observed outcome range for interpretability or robustness.

The function allows for flexible control of HAL fitting behavior, including
options for data-adaptive smoothing or manually setting smoothness order and
number of knots. It is particularly useful in causal inference or personalized
medicine settings where estimating treatment effect heterogeneity is of interest.
}
