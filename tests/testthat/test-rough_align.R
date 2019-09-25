# First, set up test objects... Unfortunately, these should be at approximately
# the same resolution as the shoe images (for now) because the parameters don't
# scale well

set.seed(12091038)
img_src <- matrix(1, 1000, 1000)

target_region <- img_src
target_region[450:550,200:800] <- 0

img <- EBImage::Image(target_region)

mask <- 1 - EBImage::erode(img, EBImage::makeBrush(51, "disc"))

img_blotches <- (img_src + rnorm(1000^2, sd = 1.75)) %>%
  EBImage::Image() %>%
  EBImage::gblur(sigma = 15, radius = 101)

noise_brush <- EBImage::makeBrush(5, "disc")
noise_brush <- noise_brush/sum(noise_brush)
img_noise <- abs(0*img_src + rnorm(1000^2, sd = 1)) %>%
  EBImage::Image() %>%
  EBImage::filter2(noise_brush) %>%
  EBImage::thresh(w = 5, h = 5) %>%
  EBImage::opening(EBImage::makeBrush(3, "disc"))

angles <- seq(5, 175, 10)
img_rot <- lapply(angles, function(x) {
  z <- EBImage::rotate(img, x, bg.col = 1, output.dim = dim(img_src))
  z <- EBImage::normalize(z) * img_blotches

  (1 - .6*(1 - z)*img_noise) %>%
    EBImage::gblur(sigma = 3, radius = 5)
})

shift_size <- 25
img_rot_shift <- lapply(angles, function(x) {
  z <- EBImage::rotate(img, x, bg.col = 1, output.dim = dim(img_src))
  z <- EBImage::translate(z, sample(c(shift_size, 0, -shift_size), 2, replace = T), bg.col = 1)
  z <- EBImage::normalize(z) * img_blotches
  (1 - .6*(1 - z)*img_noise) %>%
    EBImage::gblur(sigma = 3, radius = 5)
})

### End object setup ###


exag_imgs <- exaggerate_img_to_mask(img_rot, gaussian_d = 15, threshold_val = 0.14, opening_d = 11, closing_d = 101)
exag_imgs_shift <- exaggerate_img_to_mask(img_rot_shift, gaussian_d = 15, threshold_val = 0.14, opening_d = 11, closing_d = 101)

# par(mfrow = c(6, 6))
# purrr::walk(exag_imgs, plot)
# purrr::walk(exag_imgs_shift, plot)

test_that("exag images have balance of white/black pixels", {
  prop_white <- purrr::map_dbl(exag_imgs, mean)
  expect_true(all(prop_white < .1))
  prop_white_shift <- purrr::map_dbl(exag_imgs_shift, mean)
  expect_true(all(prop_white_shift < .1))
})


est_angles <- align_prcomp(exag_imgs) %>% as.numeric()
est_angles_shift <- align_prcomp(exag_imgs_shift) %>% as.numeric()

test_that("estimated angles are reasonable", {
  errs <- 90 - abs(90 - (angles + est_angles) %% 180) # get abs. distance from 90 deg, then compare to 90 deg. 180 and 0 should both be 90 deg away from 90.
  expect_true(all(errs < 2.5))
  errs_shift <- 90 - abs(90 - (angles + est_angles) %% 180) # get abs. distance from 90 deg, then compare to 90 deg. 180 and 0 should both be 90 deg away from 90.
  expect_true(all(errs_shift < 2.5))
})

est_shifts_0 <- purrr::map_df(exag_imgs, ~calc_shifts(., mask) %>% unlist() %>% rbind() %>% as.data.frame())
est_shifts_1 <- purrr::map_df(exag_imgs_shift, ~calc_shifts(., mask) %>% unlist() %>% rbind() %>% as.data.frame())

test_that("estimated shifts are reasonable", {
  expect_true(all(rowMeans(est_shifts_0)/1000 < 0.02)) # Allow for 2% error in center of image
  expect_true(all(rowMeans(est_shifts_1)/1000 < (0.02 + shift_size/1000))) # Allow for 2% error in center of image
})

test_that("estimated shifts are symmetric", {
  expect_true(all(est_shifts_0$img.top == est_shifts_0$mask.bottom))
  expect_true(all(est_shifts_1$img.top == est_shifts_1$mask.bottom))
  expect_true(all(est_shifts_0$img.left == est_shifts_0$mask.right))
  expect_true(all(est_shifts_1$img.left == est_shifts_1$mask.right))
  expect_true(all(est_shifts_0$mask.top == est_shifts_0$img.bottom))
  expect_true(all(est_shifts_1$mask.top == est_shifts_1$img.bottom))
  expect_true(all(est_shifts_0$mask.left == est_shifts_0$img.right))
  expect_true(all(est_shifts_1$mask.left == est_shifts_1$img.right))
})


rough_alignment <- lapply(img_rot, function(x)
  rough_align(x, mask,
              exaggerate_pars = list(gaussian_d = 15, threshold_val = 0.1, opening_d = 11, closing_d = 101),
              img_fill_value = 1))

rough_alignment_shift <- lapply(img_rot_shift, function(x)
  rough_align(x, mask,
              exaggerate_pars = list(gaussian_d = 15, threshold_val = 0.1, opening_d = 11, closing_d = 101),
              img_fill_value = 1))

# purrr::walk(rough_alignment, ~plot(EBImage::rgbImage(.$img, (1 - .$exag_img), 1 - .$mask)))
# purrr::walk(rough_alignment_shift, ~plot(EBImage::rgbImage(.$img, (1 - .$exag_img), 1 - .$mask)))

test_that("aligned images overlap with mask", {
  expect_true(all(purrr::map_dbl(rough_alignment, ~mean((1 - .$img)*.$mask)/mean(1 - .$img)) > .94))
  expect_true(all(purrr::map_dbl(rough_alignment_shift, ~mean((1 - .$img)*.$mask)/mean(1 - .$img)) > .94))
})

test_that("warning is generated properly", {
  expect_message(res <- rough_align(img_pad(img, c(10, 0, 0, 0), value = 1), mask), "Auto-resizing mask to the size of image")
  expect_equal(dim(res$mask), c(1000, 1016))
})



