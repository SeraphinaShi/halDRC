###############################################################################
#'  Calculate the influence curve based SE with given fitted HAL object and evaluation points
#'
#' @details The procedure calculates the SE with the following steps:
#'     1). get all the selected basis functions and corresponding from hal_fit 
#'     2). calculate influence curves for coefficients 
#'     3). use delta method to calculate influence curves for points on the curve
#' @param X An input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that will be used to derive the design matrix of basis
#'  functions.
#' @param Y A \code{numeric} vector of observations of the outcome variable.
#' @param hal_fit The HAL fit object.
#' @param family A \code{character} or a \code{\link[stats]{family}} object
#'  (supported by \code{\link[glmnet]{glmnet}}) specifying the error/link
#'  family for a generalized linear model. See more details in \code{fit_hal}.
#' @param X_unpenalized TAn input \code{matrix} with dimensions number of observations -by-
#'  number of covariates that are not penalized when fit HAL.
#'  
###############################################################################
#'  

IC_based_se <- function(X, Y, hal_fit, eval_points, family = "binomial", X_unpenalized = NULL ){
  
  n = length(Y)
  coef <- hal_fit$coefs
  basis_mat <- cbind(1, as.matrix(hal_fit$x_basis))
  
  nonzero_idx <- which(coef != 0)
  
  if(length(nonzero_idx) > 0) {
    coef_nonzero <- coef[nonzero_idx]
    basis_mat_nonzero <- as.matrix(basis_mat[, nonzero_idx])
    
    if(family == "binomial"){
      if(is.null(X_unpenalized)){
        Y_hat = predict(hal_fit, new_data = X, type = "response")
      } else {
        Y_hat = predict(hal_fit, new_data = X, new_X_unpenalized = as.matrix(X_unpenalized), type = "response")
      }
      

    } else {
      if(is.null(X_unpenalized)){
        Y_hat = predict(hal_fit, new_data = X)
      } else {
        Y_hat = predict(hal_fit, new_X_unpenalized = as.matrix(X_unpenalized), new_data = X)
      }
      
    }
    
    IC_beta <- cal_IC_for_beta(X = basis_mat_nonzero, 
                               Y = Y, 
                               Y_hat =Y_hat,
                               beta_n = coef_nonzero,
                               family = family,
                               X_unpenalized = X_unpenalized
    )
    
    if(any(! is.na(IC_beta))){
      se <- c()
      
      for (i in 1:length(eval_points)) {
        X_new <- X
        X_new$A = eval_points[i]
        
        # efficient influence curve
        x_basis_a <- hal9001::make_design_matrix(as.matrix(X_new), hal_fit$basis_list, p_reserve = 0.75)
        if(is.null(X_unpenalized)){
          x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a)[, nonzero_idx])
        } else {
          x_basis_a_nonzero <- as.matrix(cbind(1, x_basis_a, as.matrix(X_unpenalized))[, nonzero_idx])
        }
        
        IC_EY <- cal_IC_for_EY(X_new = x_basis_a_nonzero, 
                               beta_n = coef_nonzero,
                               IC_beta = IC_beta,
                               family = family)
        
        # empirical SE
        se[i] <- sqrt(var(IC_EY)/n)
        
      }
    } else {
      se <- NA
    }
    
  } else {
    se <- NA
  }
  
  return(se)
}



###############################################################################
# calculating efficient influence curves

cal_IC_for_beta <- function(X, Y, Y_hat, beta_n, family = 'binomial', X_unpenalized = NULL){
  n <- dim(X)[1] 
  p <- length(beta_n)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # 1. calculate score: X'(Y - phi(X))
  res <- Y-Y_hat
  score <- sweep(t(X), 2, res, `*`)
  
  # 2. calculate the derivative of phi:
  if(family == 'binomial'){
    d_phi_scaler <- as.vector(exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2)) # exp(- beta_n %*% t(X)) / ((1 + exp(- beta_n %*% t(X)))^2))
    d_phi <- sweep(X, 1, d_phi_scaler, `*`)
  } else {
    d_phi = - X
  }
  
  # 3. -E_{P_n}(X d_phi)^(-1)
  tmat <- t(X) %*% d_phi / n
  if(! is.matrix(try(solve(tmat), silent = TRUE))){
    return(NA)
  }
  tmat <- -solve(tmat)
  
  # 4. calculate influence curves
  IC <- tmat %*% score
  
  return(IC)
}


cal_IC_for_EY <- function(X_new, beta_n, IC_beta, family = 'binomial'){
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  if(family == 'binomial'){
    d_phi_scaler_new <- as.vector(exp(- beta_n %*% t(X_new)) / ((1 + exp(- beta_n %*% t(X_new)))^2))
    d_phi_new <- sweep(X_new, 1, d_phi_scaler_new, `*`)
  } else {
    d_phi_new = X_new
  }
  
  IC = colMeans(d_phi_new %*% IC_beta)
  
  return(IC)
}


cal_IC_for_ATE <- function(X_new_a, X_new_0, beta_n, IC_beta, family = 'binomial'){
  
  if (!is.matrix(X_new_a)) X_new_a <- as.matrix(X_new_a)
  if (!is.matrix(X_new_0)) X_new_0 <- as.matrix(X_new_0)
  
  if (family == 'binomial') {
    d_phi_scaler_new_a <- as.vector(exp(- beta_n %*% t(X_new_a)) / ((1 + exp(- beta_n %*% t(X_new_a)))^2))
    d_phi_new_a <- sweep(X_new_a, 1, d_phi_scaler_new_a, `*`)
    
    d_phi_scaler_new_0 <- as.vector(exp(- beta_n %*% t(X_new_0)) / ((1 + exp(- beta_n %*% t(X_new_0)))^2))
    d_phi_new_0 <- sweep(X_new_0, 1, d_phi_scaler_new_0, `*`)
    
    d_phi_new <- d_phi_new_a - d_phi_new_0
  } else {
    d_phi_new <- X_new_a - X_new_0
  }
  
  IC = colMeans(d_phi_new %*% IC_beta)
  
  return(IC)
}

