test_that("noise_ceiling delegates to model_compare", {
  set.seed(200)
  a1 <- array(rnorm(5 * 4 * 3), dim = c(5, 4, 3))
  a2 <- a1 + rnorm(length(a1), sd = 0.1)

  out <- noise_ceiling(a1, a2, comparison = "fullcompare_compthenavg")
  expect_true(all(c("corr_vals", "R2_vals", "mae_vals") %in% names(out)))
  expect_length(out$corr_vals, 3)
})

test_that("%||% returns fallback only for NULL", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(0 %||% 5, 0)
  expect_equal(FALSE %||% TRUE, FALSE)
})

test_that(".onLoad hook is callable and invisible", {
  expect_invisible(.onLoad("actflower", "actflower"))
})
