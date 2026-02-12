test_that(".af_normalize_orientation handles aliases and rejects invalid values", {
  expect_equal(.af_normalize_orientation("nodes_by_time"), "nodes_by_time")
  expect_equal(.af_normalize_orientation("nodes_x_time"), "nodes_by_time")
  expect_equal(.af_normalize_orientation("time_x_nodes"), "time_by_nodes")

  expect_error(.af_normalize_orientation(""), "non-empty character scalar")
  expect_error(.af_normalize_orientation("bogus"), "Unknown `orientation`")
})

test_that("as_nodes_by_time validates numeric matrix-like inputs", {
  expect_error(as_nodes_by_time(1:3), "numeric matrix or a supported neuro object")
  expect_error(as_nodes_by_time(matrix("a", 2, 2)), "numeric matrix or a supported neuro object")
})
