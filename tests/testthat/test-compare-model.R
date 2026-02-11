test_that("model_compare returns summary vectors for fullcompare", {
  set.seed(3)
  arr1 <- array(rnorm(5 * 4 * 3), dim = c(5, 4, 3))
  arr2 <- arr1 + rnorm(length(arr1), sd = 0.05)
  out <- model_compare(arr1, arr2, comparison = "fullcompare_compthenavg")
  expect_true(all(c("corr_vals", "R2_vals", "mae_vals") %in% names(out)))
  expect_length(out$corr_vals, 3)
})

test_that("model_compare supports conditionwise and nodewise compare-then-average", {
  set.seed(4)
  target <- array(rnorm(6 * 5 * 4), dim = c(6, 5, 4))
  pred <- target + rnorm(length(target), sd = 0.1)

  out_cond <- model_compare(target, pred, comparison = "conditionwise_compthenavg")
  expect_equal(dim(out_cond$corr_vals), c(6, 4))
  expect_equal(dim(out_cond$R2_vals), c(6, 4))
  expect_equal(dim(out_cond$mae_vals), c(6, 4))

  out_node <- model_compare(target, pred, comparison = "nodewise_compthenavg")
  expect_equal(dim(out_node$corr_vals), c(5, 4))
  expect_equal(dim(out_node$R2_vals), c(5, 4))
  expect_equal(dim(out_node$mae_vals), c(5, 4))
})

test_that("model_compare supports average-then-compare modes", {
  set.seed(5)
  target <- array(rnorm(7 * 4 * 3), dim = c(7, 4, 3))
  pred <- target + rnorm(length(target), sd = 0.1)

  out_cond <- model_compare(target, pred, comparison = "conditionwise_avgthencomp")
  expect_true(all(c(
    "corr_conditionwise_avgthencomp_bynode",
    "R2_conditionwise_avgthencomp_bynode",
    "maeAcc_bynode_avgthencomp"
  ) %in% names(out_cond)))
  expect_length(out_cond$corr_conditionwise_avgthencomp_bynode, 7)

  out_node <- model_compare(target, pred, comparison = "nodewise_avgthencomp")
  expect_true(all(c(
    "corr_nodewise_avgthencomp_bycond",
    "R2_nodewise_avgthencomp_bycond",
    "maeAcc_nodewise_avgthencomp_bycond"
  ) %in% names(out_node)))
  expect_length(out_node$corr_nodewise_avgthencomp_bycond, 4)
})

test_that("model_compare appends model2 outputs with suffix", {
  set.seed(6)
  target <- array(rnorm(5 * 4 * 3), dim = c(5, 4, 3))
  model1 <- target + rnorm(length(target), sd = 0.05)
  model2 <- target + rnorm(length(target), sd = 0.15)

  out <- model_compare(target, model1, model2 = model2, comparison = "fullcompare_compthenavg")
  expect_true(all(c("corr_vals", "corr_vals_model2", "R2_vals", "R2_vals_model2") %in% names(out)))
  expect_length(out$corr_vals_model2, 3)
})
