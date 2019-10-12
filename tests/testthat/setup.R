
expect_all_lte <- function(x, y) {
  purrr::walk2(x, y, ~eval(bquote(expect_lte(.x, .y))))
}

expect_all_equal <- function(x, y) {
  purrr::walk2(x, y, ~eval(bquote(expect_equal(.x, .y))))
}
