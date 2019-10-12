context("shoe_mask")


test_that("shoe_mask works as expected", {
  expect_equal(dim(shoe_mask("Adidas", size = 10, foot = "L")), c(1998, 3906))
  expect_equal(dim(shoe_mask("Nike", size = 8, foot = "L", ppi = 200)), c(1332, 2616))
  expect_error(shoe_mask("Adidas", size = 8, foot = "L"))
  expect_error(shoe_mask("Nike", size = 7, foot = "L"))
})
