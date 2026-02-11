test_that("is_neuroim2_object detects declared neuroim2 classes", {
  x <- structure(list(), class = "NeuroVec")
  y <- matrix(1, 2, 2)
  expect_true(is_neuroim2_object(x))
  expect_false(is_neuroim2_object(y))
})

test_that("as_nodes_by_time_neuroim2 fails clearly without neuroim2", {
  skip_if(has_neuroim2(), "neuroim2 installed; this test targets missing-dependency behavior")
  x <- structure(list(), class = "NeuroVec")
  expect_error(as_nodes_by_time_neuroim2(x), "Install 'neuroim2'")
})
