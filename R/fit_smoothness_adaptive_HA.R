
fit_SL_smoothness_adaptive_HAL <- function(obs, X, Y, x_names, y_name, y_type){
  
  # 1. --- SL
  task <- make_sl3_Task(
    data = obs,
    outcome = y_name,
    covariates = x_names
  )
  
  # num_knots = c(200, 100,  50)
  lrn_hal0 <- Lrnr_hal9001$new(smoothness_orders = 0, family = y_type,
                               return_x_basis = TRUE)
  lrn_hal1 <- Lrnr_hal9001$new(smoothness_orders = 1, family = y_type,
                               return_x_basis = TRUE)
  lrn_hal2 <- Lrnr_hal9001$new(smoothness_orders = 2, family = y_type,
                               return_x_basis = TRUE)
  lrn_hal3 <- Lrnr_hal9001$new(smoothness_orders = 3, family = y_type,
                               return_x_basis = TRUE)
  
  # num_knots = c(20,10,5)
  lrn_hal0_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 0,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 0,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal1_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 1,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 1,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal2_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 2,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 2,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  lrn_hal3_smallerKnots <- Lrnr_hal9001$new(smoothness_orders = 3,
                                            family = y_type, 
                                            return_x_basis = TRUE, 
                                            num_knots = hal9001:::num_knots_generator(
                                              max_degree = ifelse(ncol(X) >= 20, 2, 3), 
                                              smoothness_orders = 3,
                                              base_num_knots_0 = 20,
                                              base_num_knots_1 = 20  
                                            ))
  
  learners <- c(lrn_hal0, lrn_hal0_smallerKnots, lrn_hal1, lrn_hal1_smallerKnots, 
                lrn_hal2, lrn_hal2_smallerKnots, lrn_hal3, lrn_hal3_smallerKnots)
  names(learners) <- c("halcv0", "halcv0_sKnots", "halcv1", "halcv1_sKnots", 
                       "halcv2", "halcv2_sKnots", "halcv3", "halcv3_sKnots")
  stack <- make_learner(Stack, learners)

  cv_selector <- Lrnr_cv_selector$new(eval_function = loss_loglik_binomial) 
  # https://tlverse.org/sl3/reference/loss_functions.html
  dSL <- Lrnr_sl$new(learners = stack, metalearner = cv_selector)

  dSL_fit <- dSL$train(task)

  return(list(dSL_fit = dSL_fit,
              smooth_orders_candidates = c(0, 0, 1, 1, 2, 2, 3, 3),
              num_knots_candidates = c("default", "20", "default", "20",
                                       "default", "20", "default", "20")))
  
}
