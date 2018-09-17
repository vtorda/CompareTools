context("LadderizeTree")

test_that("ladderizing a phylo object", {
  tree1 <- LadderizeTree(fungi_tree)
  expect_equal(tree1, fungi_tree_lad)
})
