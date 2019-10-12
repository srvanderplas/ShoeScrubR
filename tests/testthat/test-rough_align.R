### Set up test objects ###

set.seed(12091538)
img <- EBImage::readImage("test_bkgd.tif", all = F)[,,1]
EBImage::colorMode(img) <- "grayscale"
im_center <- round(dim(img)/2)

w <- 31
h <- 181
x <- rnorm(w*h, mean = .4, sd = .2) %>% matrix(nrow = w, ncol = h)
x <- EBImage::gblur(x, sigma = 1) # Introduce some dependence
x <- x + sample(c(0, 1), size = w*h, replace = T, prob = c(.8, .2)) # Add some censoring
x <- pmax(0, pmin(x, 1))

bar_img <- img*0 + 1
bar_img[im_center[1] + seq(-floor(w/2), floor(w/2), 1),
        im_center[2] + seq(-floor(h/2), floor(h/2), 1)] <- x

angles <- seq(5, 175, 20)
img_rot <- lapply(angles, function(x) {
  z <- EBImage::rotate(bar_img, x, output.dim = dim(img), bg.col = 1)
  zz <- img
  zz[z < 1] <- z[z < 1]
  zz
})
shift_size <- 15
img_rot_shift <- lapply(angles, function(x) {
  z <- EBImage::rotate(bar_img, x, output.dim = dim(img), bg.col = 1) %>%
    EBImage::translate(v = sample(c(-shift_size, 0, shift_size), 2, replace = T), bg.col = 1)
  zz <- img
  zz[z < 1] <- z[z < 1]
  zz
})

mask <- 0*img
mask[im_center[1] + seq(-floor(w/2) - floor(w/10), floor(w/2) + floor(w/10), 1),
     im_center[2] + seq(-floor(h/2) - floor(w/10), floor(h/2) + floor(w/10), 1)] <- 1

# par(mfrow = c(6, 6))
# purrr::walk(img_rot, plot)
# purrr::walk(img_rot_shift, plot)

### End object setup ###



exag_imgs1 <- ShoeScrubR:::exaggerate_img_control(img_rot, gaussian_d = 15, threshold_val = 0.12, opening_d = 5, closing_d = 55)
exag_imgs_shift1 <- ShoeScrubR:::exaggerate_img_control(img_rot_shift, gaussian_d = 15, threshold_val = 0.12, opening_d = 5, closing_d = 55)

# par(mfrow = c(6, 6))
# purrr::walk(exag_imgs1, plot)
# purrr::walk(exag_imgs_shift1, plot)

exag_imgs2 <- ShoeScrubR:::exaggerate_img_auto(img_rot, ngroups = 2, fast = T, epsilon = 1e-04)
exag_imgs_shift2 <- ShoeScrubR:::exaggerate_img_auto(img_rot_shift, ngroups = 2, fast = T)

# par(mfrow = c(6, 6))
# purrr::walk(exag_imgs2, plot)
# purrr::walk(exag_imgs_shift2, plot)

test_that("parameter exag images have balance of white/black pixels", {
  prop_white <- purrr::map_dbl(exag_imgs1, mean)
  prop_white_shift <- purrr::map_dbl(exag_imgs_shift1, mean)

  expect_all_lte(prop_white, .1)
  expect_all_lte(prop_white_shift, .1)
})

test_that("auto exag images have balance of white/black pixels", {
  prop_white <- purrr::map_dbl(exag_imgs2, mean)
  prop_white_shift <- purrr::map_dbl(exag_imgs_shift2, mean)

  expect_all_lte(prop_white, .4)
  expect_all_lte(prop_white_shift, .4)
})

exag_imgs <- c(exag_imgs1, exag_imgs2)
exag_imgs_shift <- c(exag_imgs_shift1, exag_imgs_shift2)

test_that("pca alignment is reasonable", {

  est_angles <- align_prcomp(exag_imgs) %>% as.numeric()
  est_angles_shift <- align_prcomp(exag_imgs_shift) %>% as.numeric()

  errs <- 90 - abs(90 - (rep(angles, times = 2) + est_angles) %% 180) # get abs. distance from 90 deg, then compare to 90 deg. 180 and 0 should both be 90 deg away from 90.
  errs_shift <- 90 - abs(90 - (rep(angles, times = 2) + est_angles) %% 180) # get abs. distance from 90 deg, then compare to 90 deg. 180 and 0 should both be 90 deg away from 90.

  expect_all_lte(errs, 5)
  expect_all_lte(errs_shift, 5)

  expect_lt(abs(-45 - align_prcomp(data.frame(row = 1:50, col = -1*(1:50) + rnorm(50)))), 2)
})

test_that("pca_to_angle works", {
  expect_equal(pca_to_angle(matrix(c(-cos(pi/3), sin(pi/3), -sin(pi/3), -cos(pi/3)), nrow = 2, byrow = T)), 30)
})


test_that("estimated shifts are reasonable and symmetric", {
  est_shifts_0 <- purrr::map_df(exag_imgs, ~calc_shifts(., mask) %>% unlist() %>% rbind() %>% as.data.frame())
  est_shifts_1 <- purrr::map_df(exag_imgs_shift, ~calc_shifts(., mask) %>% unlist() %>% rbind() %>% as.data.frame())

  expect_all_lte(sapply((est_shifts_0)/nrow(img), max), 0.03) # Allow for 3% error in center of image
  expect_all_lte(sapply((est_shifts_1)/nrow(img), max),  (0.03 + 1.2*shift_size/nrow(img)))# Allow for 2% error in center of image

  expect_all_equal(est_shifts_0$img.top, est_shifts_0$mask.bottom)
  expect_all_equal(est_shifts_1$img.top, est_shifts_1$mask.bottom)
  expect_all_equal(est_shifts_0$img.left, est_shifts_0$mask.right)
  expect_all_equal(est_shifts_1$img.left, est_shifts_1$mask.right)
  expect_all_equal(est_shifts_0$mask.top, est_shifts_0$img.bottom)
  expect_all_equal(est_shifts_1$mask.top, est_shifts_1$img.bottom)
  expect_all_equal(est_shifts_0$mask.left, est_shifts_0$img.right)
  expect_all_equal(est_shifts_1$mask.left, est_shifts_1$img.right)
})

test_that("aligned images overlap with mask", {

  rough_alignment <- lapply(img_rot, function(x)
    rough_align(x, mask, img_fill_value = 1, ngroups = 2))

  rough_alignment_shift <- lapply(img_rot_shift, function(x)
    rough_align(x, mask, img_fill_value = 1, ngroups = 2))

  # par(mfrow = c(6, 6))
  # purrr::walk(rough_alignment, ~plot(EBImage::rgbImage(.$img, (1 - .$exag_img), 1 - .$mask)))
  # purrr::walk(rough_alignment_shift, ~plot(EBImage::rgbImage(.$img, (1 - .$exag_img), 1 - .$mask)))

  expect_all_lte((mean(bar_img*(1 - mask))/mean(img) - .01),
                 purrr::map_dbl(rough_alignment, ~mean(.$img*(1-.$mask))/mean(.$img)))
  expect_all_lte((mean(bar_img*(1 - mask))/mean(img) - .01),
                 purrr::map_dbl(rough_alignment_shift, ~mean(.$img*(1-.$mask))/mean(.$img)))
})

test_that("warning is generated properly when mask must be resized during alignment", {
  expect_message(res <- rough_align(img_pad(img + bar_img, c(30, 0, 0, 0), value = 1), mask, ngroups = 2),
                 "Auto-resizing mask to the size of image")
})

test_that("get_mask_arch works", {
  test_img <- cbind(EBImage::makeBrush(25, "disc"), EBImage::makeBrush(25, "disc")) %>%
    EBImage::as.Image() %>%
    img_pad(padding = c(15, 21, 20, 18), value = 0)
  # test_center <- round(dim(test_img)/2)
  mask_min <- get_mask_arch(test_img)
  expect_equal(mask_min, c(33, 41))

  test_img <- EBImage::as.Image(EBImage::makeBrush(25, "line", angle = 45))
  mask_min <- get_mask_arch(test_img)
  expect_equal(mask_min, c(12, 12))
})

test_that("em_thresh works", {
  test_img <- img
  expect_silent(em_thresh(img, ngroups = 2, scale_factor = NULL, N = NULL))
  expect_silent(em_thresh(img, ngroups = 3, scale_factor = NULL, N = 300))
  expect_silent(em_thresh(list(img, img)))

  set.seed(50920803)
  expect_message(em_thresh(img*0 + sample(c(0, 1), length(img), replace = T) * rnorm(length(img), 1, .1)),
                 "EM mixture did not succeed with 3 groups. Retrying with 2 groups.")
  expect_error(
    expect_message(em_thresh(img*0 + sample(c(0, 1), length(img), replace = T, prob = c(.8, .2)) * rnorm(length(img), 1, .1)),
               "Mixture model fitting failed")
  )
})

test_that("est_ppi_film works", {
  expect_equal(est_ppi_film(c(705,1343)), 100)
})

test_that("clean_initial_img threshold search works", {
  expect_message(clean_initial_img(img, threshold_val = 0), "Image cleaned using a threshold of")
})

test_that("binary_center search works", {
  res <- binary_center(list(shoe_mask(brand = "Nike", size = 10, foot = "L"),
                            shoe_mask(brand = "Adidas", size = 10, foot = "L")))
  expect_equal(res[[1]], c(1068, 1890))
  expect_equal(res[[2]], c(958, 1886))

})