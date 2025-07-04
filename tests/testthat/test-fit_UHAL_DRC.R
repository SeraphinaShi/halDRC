test_that("fit_UHAL_DRC runs without error", {
  dat <- data.frame(
    Y = rbinom(100, 1, 0.5),
    A = runif(100),
    X1 = rnorm(100),
    X2 = rnorm(100)
  )
  res <- fit_UHAL_DRC(dat, y_var_name = "Y", trt_var_name = "A", family = "binomial")
  expect_true("curve_est" %in% names(res))
})
