% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/undersmooth_hal.R
\name{get_maxscore}
\alias{get_maxscore}
\title{undersoomthed HAL helper function for global criterion}
\usage{
get_maxscore(basis_mat, resid_mat, sd_est, Nlam, us_fit)
}
\arguments{
\item{basis_mat}{The selected basis matrix from initial fit for undersmoothing,
obtained from the output list of \code{undersmooth_init}.}

\item{resid_mat}{The residual matrix with each column the residuals correspongding to a lambda.}

\item{sd_est}{A numeric vector containing the sd of each column of \code{basis_mat}.}

\item{Nlam}{Number of lambda candidates.}

\item{us_fit}{The \code{glmnet} fit of the sequence of candidate lambdas.}
}
\description{
undersoomthed HAL helper function for global criterion
}
\details{
For each candidate lambda, do:
    1). standardize the score formed by each basis.
    2). calculate the mean of the standardized scores for each basis.
Select the max of the mean.
}
