test_that("validate_matrix enforces numeric matrix", {
  expect_error(validate_matrix(1:3), "numeric matrix")
  expect_invisible(validate_matrix(matrix(1, 2, 2)))
})
