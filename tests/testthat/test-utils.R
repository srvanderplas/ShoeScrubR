
test_img <- matrix(0, nrow = 10, ncol = 8)
test_white_px <- cbind(c(2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 7, 7, 7, 8, 8, 8, 8, 9, 9),
                       c(5, 6, 7, 2, 3, 4, 5, 3, 4, 5, 4, 5, 6, 7, 5, 2, 5, 7, 4, 5, 6, 7, 3, 4))
test_img[test_white_px] <- 1
test_img <- test_img %>% EBImage::Image()


test_that("image_to_df works", {
  df <- image_to_df(test_img) %>%
    dplyr::arrange(row, col)

  expect_equal(df$row, test_white_px[,1])
  expect_equal(df$col, test_white_px[,2])
  expect_equal(df$frame, rep(1, nrow(test_white_px)))
  expect_equal(df$value, df$frame)
  expect_equal(attr(df, "operation"),
               list(list(type = "convert to df", filter_val = 0)))
})

test_that("img_pad_to_center works", {
  pad12 <- img_pad_to_center(test_img, c(1, 2), .5)
  expect_equal(EBImage::imageData(pad12)[1:4,1:16],
               matrix(.5, 4, 16, dimnames = list(NULL, NULL)))
  expect_equal(EBImage::imageData(pad12)[1:14, 1:8],
               matrix(.5, 14, 8, dimnames = list(NULL, NULL)))
  expect_equivalent(EBImage::imageData(pad12)[5:14, 9:16],
                    EBImage::imageData(test_img))
  expect_warning(img_pad_to_center(test_img, c(2.5, 6)))
  expect_error(img_pad_to_center(test_img, c(-5, 5)))
  expect_equal(attr(pad12, "operation"),
               list(list(type = "pad", top_bottom = c(8, 0), left_right = c(4, 0),
                         value = 0.5)))
})

test_that("img_resize works", {
  resize_test <- img_resize(test_img, w = ncol(test_img)*2, h = nrow(test_img)*2)
  resize_ebimage <- EBImage::resize(test_img, w = ncol(test_img)*2, h = nrow(test_img)*2)

  expect_equal(EBImage::imageData(resize_test), EBImage::imageData(resize_ebimage))
  expect_equal(attr(resize_test, "operation"),
               list(list(type = "resize", orig_dim = dim(test_img), final_dim = c(16, 20), other_args = list(w = ncol(test_img)*2, h = nrow(test_img)*2))))
})

test_that("img_rotate works", {
  rotate_test <- img_rotate(test_img, angle = 45)
  rotate_ebimage <- EBImage::rotate(test_img, angle = 45)

  expect_equal(EBImage::imageData(rotate_test), EBImage::imageData(rotate_ebimage))
  expect_equal(attr(rotate_test, "operation"),
               list(list(type = "rotate", angle = 45, other_args = list())))
})


test_that("img_translate works", {
  translate_test <- img_translate(test_img, v = c(2, 1))
  translate_ebimage <- EBImage::translate(test_img, v = c(2, 1))

  expect_equal(EBImage::imageData(translate_test), EBImage::imageData(translate_ebimage))
  expect_equal(attr(translate_test, "operation"),
               list(list(type = "translate", vector = c(2, 1), other_args = list())))
})