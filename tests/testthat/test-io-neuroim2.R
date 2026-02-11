test_that("is_neuro_object detects supported neuro classes", {
  x <- structure(list(), class = "NeuroVec")
  y <- matrix(1, 2, 2)
  expect_true(is_neuro_object(x))
  expect_false(is_neuro_object(y))
})

test_that("as_nodes_by_time handles supported neuro objects via matrix fallback", {
  x <- structure(list(X = matrix(1:6, nrow = 2, ncol = 3)), class = "NeuroVec")
  out <- as_nodes_by_time(x, orientation = "nodes_by_time")
  expect_equal(out, matrix(1:6, nrow = 2, ncol = 3))
})

test_that("as_nodes_by_time transposes supported neuro objects when requested", {
  x <- structure(list(X = matrix(1:6, nrow = 2, ncol = 3)), class = "NeuroVec")
  out <- as_nodes_by_time(x, orientation = "time_by_nodes")
  expect_equal(out, t(matrix(1:6, nrow = 2, ncol = 3)))
})

test_that("as_nodes_by_time errors on unsupported non-matrix object", {
  x <- structure(list(foo = 1), class = "not_neuro")
  expect_error(as_nodes_by_time(x), "numeric matrix or a supported neuro object")
})
