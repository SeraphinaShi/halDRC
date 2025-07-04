
summary_hal9001_unpenalized <- function(object,
                                        lambda = NULL,
                                        only_nonzero_coefs = TRUE,
                                        include_redundant_terms = FALSE,
                                        round_cutoffs = 3,
                                        ...) {
  abs_coef <- basis_list_idx <- coef_idx <- dup <- NULL
  
  # retain coefficients corresponding to lambda
  if (!is.null(lambda)) {
    if (length(lambda) > 1) {
      stop("Cannot summarize over multiple values of lambda.")
    }
    if (lambda != object$lambda_star) {
      if (is.null(object$lasso_fit)) {
        stop(
          "Coefficients for specified lambda do not exist, or are not ",
          "accessible since the fit of the lasso model was not returned ",
          "(i.e., return_lasso was set to FALSE in `hal_fit()`)."
        )
      } else {
        if (!(lambda %in% object$lasso_fit$lambda)) {
          stop("Coefficients for the specified lambda do not exist.")
        } else {
          lambda_idx <- which(object$lasso_fit$lambda == lambda)
          coefs <- object$lasso_fit$glmnet.fit$beta[, lambda_idx]
        }
      }
    } else {
      lambda_idx <- which(object$lambda_star == lambda)
      coefs <- object$coefs[, lambda_idx]
    }
  }
  
  if (is.null(lambda)) {
    lambda <- object$lambda_star
    coefs <- object$coefs
    if (length(lambda) > 1) {
      warning(
        "Coefficients for many lambda exist --\n",
        "Summarizing coefficients corresponding to minimum lambda."
      )
      lambda_idx <- which.min(lambda)
      coefs <- object$coefs[, lambda_idx]
    }
  }
  
  # cox model has no intercept
  if (object$family != "cox") {
    coefs_no_intercept <- coefs[-1]
  } else {
    coefs_no_intercept <- coefs
  }
  
  # subset to non-zero coefficients
  if (only_nonzero_coefs) {
    coef_idxs <- which(coefs_no_intercept != 0)
  } else {
    coef_idxs <- seq_along(coefs_no_intercept)
  }
  copy_map <- object$copy_map[coef_idxs]
  
  if (object$unpenalized_covariates > 0){
    copy_map <- copy_map[1:(length(copy_map) - object$unpenalized_covariates)]
  }
  
  # summarize coefficients with respect to basis list
  coefs_summ <- data.table::rbindlist(
    lapply(seq_along(copy_map), function(map_idx) {
      coef_idx <- coef_idxs[map_idx]
      coef <- coefs_no_intercept[coef_idx]
      
      basis_list_idxs <- copy_map[[map_idx]] # indices of duplicates
      basis_dups <- object$basis_list[basis_list_idxs]
      
      data.table::rbindlist(
        lapply(seq_along(basis_dups), function(i) {
          coef_idx <- ifelse(object$family != "cox", coef_idx + 1, coef_idx)
          dt <- data.table::data.table(
            coef_idx = coef_idx, # coefficient index
            coef, # coefficient
            basis_list_idx = basis_list_idxs[i], # basis list index
            col_idx = basis_dups[[i]]$cols, # column idx in X
            col_cutoff = basis_dups[[i]]$cutoffs, # cutoff
            col_order = basis_dups[[i]]$orders # smoothness order
          )
          return(dt)
        })
      )
    })
  )
  
  if (!include_redundant_terms) {
    coef_idxs <- unique(coefs_summ$coef_idx)
    coefs_summ <- data.table::rbindlist(lapply(coef_idxs, function(idx) {
      # subset to matching coefficient index
      coef_summ <- coefs_summ[coef_idx == idx]
      
      # label duplicates (i.e. basis functions with identical col & cutoff)
      dups_tbl <- coef_summ[, c("col_idx", "col_cutoff", "col_order")]
      if (!anyDuplicated(dups_tbl)) {
        return(coef_summ)
      } else {
        # add col indicating whether or not there is a duplicate
        coef_summ[, dup := (duplicated(dups_tbl) |
                              duplicated(dups_tbl, fromLast = TRUE))]
        
        # if basis_list_idx contains redundant duplicates, remove them
        redundant_dups <- coef_summ[dup == TRUE, "basis_list_idx"]
        if (nrow(redundant_dups) > 1) {
          # keep the redundant duplicate term that has the shortest length
          retain_idx <- which.min(apply(redundant_dups, 1, function(idx) {
            nrow(coef_summ[basis_list_idx == idx])
          }))
          idx_keep <- unname(unlist(redundant_dups[retain_idx]))
          coef_summ <- coef_summ[basis_list_idx == idx_keep]
        }
        return(coef_summ[, -"dup"])
      }
    }))
  }
  
  # summarize with respect to x column names:
  x_names <- data.table::data.table(
    col_idx = 1:length(object$X_colnames),
    col_names = object$X_colnames
  )
  summ <- merge(coefs_summ, x_names, by = "col_idx", all.x = TRUE)
  
  # combine name, cutoff into 0-order basis function (may include interaction)
  summ$zero_term <- paste0(
    "I(", summ$col_names, " >= ", round(summ$col_cutoff, round_cutoffs), ")"
  )
  summ$higher_term <- ifelse(
    summ$col_order == 0, "",
    paste0(
      "(", summ$col_names, " - ",
      round(summ$col_cutoff, round_cutoffs), ")"
    )
  )
  summ$higher_term <- ifelse(
    summ$col_order < 1, summ$higher_term,
    paste0(summ$higher_term, "^", summ$col_order)
  )
  summ$term <- ifelse(
    summ$col_order == 0,
    paste0("[ ", summ$zero_term, " ]"),
    paste0("[ ", summ$zero_term, "*", summ$higher_term, " ]")
  )
  
  term_tbl <- data.table::as.data.table(stats::aggregate(
    term ~ basis_list_idx,
    data = summ, paste, collapse = " * "
  ))
  
  # no longer need the columns or rows that were incorporated in the term
  redundant <- c(
    "term", "col_cutoff", "col_names", "col_idx", "col_order", "zero_term",
    "higher_term"
  )
  summ <- summ[, -..redundant]
  summ_unique <- unique(summ)
  summ <- merge(
    term_tbl, summ_unique,
    by = "basis_list_idx", all.x = TRUE, all.y = FALSE
  )
  
  # summarize in a list
  coefs_list <- lapply(unique(summ$coef_idx), function(this_coef_idx) {
    coef_terms <- summ[coef_idx == this_coef_idx]
    list(coef = unique(coef_terms$coef), term = t(coef_terms$term))
  })
  
  # summarize in a table
  coefs_tbl <- data.table::as.data.table(stats::aggregate(
    term ~ coef_idx,
    data = summ, FUN = paste, collapse = "  OR  "
  ))
  redundant <- c("term", "basis_list_idx")
  summ_unique_coefs <- unique(summ[, -..redundant])
  coefs_tbl <- data.table::data.table(merge(
    summ_unique_coefs, coefs_tbl,
    by = "coef_idx", all = TRUE
  ))
  coefs_tbl[, "abs_coef" := abs(coef)]
  coefs_tbl <- data.table::setorder(coefs_tbl[, -"coef_idx"], -abs_coef)
  coefs_tbl <- coefs_tbl[, -"abs_coef", with = FALSE]
  
  # incorporate intercept
  if (object$family != "cox") {
    intercept <- list(data.table::data.table(
      coef = coefs[1], term = "(Intercept)"
    ))
    
    
    if (object$unpenalized_covariates > 0){
      coefs_unpenalized_tbl <- list(data.table::data.table(
        coef = object$coefs[(nrow(object$coefs) - object$unpenalized_covariates + 1) : nrow(object$coefs),], 
        term = names(object$coefs[(nrow(object$coefs) - object$unpenalized_covariates + 1) : nrow(object$coefs),])
      ))
      
      coefs_tbl <- data.table::rbindlist(
        c(intercept, list(coefs_tbl), c(coefs_unpenalized_tbl)),
        fill = TRUE
      )
      
      intercept <- list(coef = coefs[1], term = "(Intercept)")
      coefs_list <- c(list(intercept), list(coefs_unpenalized_tbl), coefs_list)
    } else {
      coefs_tbl <- data.table::rbindlist(
        c(intercept, list(coefs_tbl)),
        fill = TRUE
      )
      intercept <- list(coef = coefs[1], term = "(Intercept)")
      coefs_list <- c(list(intercept), coefs_list)
    }
    
  }
  
  out <- list(
    table = coefs_tbl,
    list = coefs_list,
    lambda = lambda,
    only_nonzero_coefs = only_nonzero_coefs
  )
  class(out) <- "summary.hal9001"
  return(out)
}

