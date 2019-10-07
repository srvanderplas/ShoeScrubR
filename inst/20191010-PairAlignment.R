library(EBImage)
library(ShoeScrubR)
library(tidyverse)

set.seed(3142095)

img_output_dir <- "~/Projects/CSAFE/2019-this_is_us/images/shoes/longitudinal/"
lss_dir <- "/lss/research/csafe-shoeprints/ShoeImagingPermanent"

# For a bunch of images...
full_imglist <- list.files("/lss/research/csafe-shoeprints/ShoeImagingPermanent/",
                           pattern = "001\\d{3}[L]_\\d{8}_5_._._.*_.*_.*", full.names = T)
dir <- "/tmp/film-prints"
if (!dir.exists(dir)) dir.create(dir)

file.copy(full_imglist, file.path(dir, basename(full_imglist)))
imglist <- file.path(dir, basename(full_imglist))

shoe_info <- read_csv("~/Projects/CSAFE/2018_Longitudinal_Shoe_Project/Clean_Data/shoe-info.csv") %>%
  filter(ShoeID %in% str_sub(basename(imglist), 1, 3)) %>%
  select(ShoeID, Brand, Size) %>%
  mutate(Size = str_remove(Size, "[ MW]") %>% parse_number()) %>%
  crossing(tibble(Mask_foot = c("R", "L"), Shoe_foot = c("L", "R")), ppi = c(200, 300)) %>%
  mutate(mask = purrr::pmap(list(Brand, Size, Mask_foot, ppi = ppi), shoe_mask))


scan_info <- tibble(
  file = imglist,
  ShoeID = str_extract(basename(file), "^\\d{3}"),
  Shoe_foot = str_extract(basename(file), "\\d{6}[RL]") %>% str_remove_all("\\d"),
  date = str_extract(basename(file), "\\d{8}") %>% parse_date(format = "%Y%m%d"),
  rep = str_extract(basename(file), "5_[12]_1") %>% str_remove("5_|_1")
) %>%
  left_join(unique(select(shoe_info, ShoeID, Brand, Size, Shoe_foot))) %>%
  mutate(
    img = purrr::map(file, EBImage::readImage, all = F),
    img = purrr::map(img, EBImage::channel, "luminance"),
    im_dim = purrr::map(img, dim)
  )

par(mfrow = c(1, 6))
purrr::walk(scan_info$img, plot)

scan_info <- scan_info %>%
  mutate(ppi = purrr::map_dbl(im_dim, est_ppi_film)) %>%
  left_join(shoe_info)

scan_info <- scan_info %>%
  mutate(align = purrr::map2(img, mask, rough_align))


plot_align <- function(df) {
  thresh_intersect <- 1 - thresh((1 - df$img)*df$mask, w = 5, h = 5, offset = 0.02)
  rgbImage(1 - df$mask, df$img, thresh_intersect) %>% plot()
}

par(mfrow = c(1, 6))
purrr::walk(scan_info$align, plot_align)

scan_info <- scan_info %>%
  mutate(aligned_img = purrr::map(align, "img"),
         aligned_exag_img = purrr::map(align, "exag_img"),
         clean_img = purrr::map2(aligned_img, aligned_exag_img, ~.x*(.y != 0) + img_mode(.x)*(.y==0))) %>%
  mutate(aligned_pyr = purrr::map(clean_img, img_pyramid, scale = c(2, 3, 4, 5, 6, 7, 8)))


library(RNiftyReg)
regs <- purrr::map2(scan_info$aligned_pyr[[1]]$img, scan_info$aligned_pyr[[2]]$img, niftyreg, scope = "rigid")

reg1 <- niftyreg(scan_info$aligned_pyr[[1]]$img[[3]], scan_info$aligned_pyr[[2]]$img[[3]], scope = "rigid")
reg2 <- niftyreg(scan_info$aligned_pyr[[1]]$img[[2]], scan_info$aligned_pyr[[2]]$img[[2]], scope = "rigid")
kernel <- EBImage::makeBrush(3, "disc")
gradient1 <- dilate(reg1$image, kernel) - erode(reg1$image, kernel)
gradient2 <- dilate(reg2$image, kernel) - erode(reg2$image, kernel)

transforms <- purrr::map(regs, "forwardTransforms")
purrr::map2(transforms, c(2, 3, 4, 5, 6, 7, 8), ~.x[[1]][1:2, 4]*.y)

align_scaled_matrix <- function(im1, im2, scale, ...) {

  res <- niftyreg(im1, im2, scope = "rigid")

  scale_mat <- matrix(1, nrow = 4, ncol = 4)
  scale_mat[1:2,4] <- scale

  res$forwardTransforms[[1]] <- res$forwardTransforms[[1]]*scale_mat
  res$reverseTransforms[[1]] <- res$reverseTransforms[[1]]*scale_mat

  im1transmat <- matrix(0, nrow = 3, ncol = 2)
  im1transmat[1:2,] <- t(res$reverseTransforms[[1]][1:2, 1:2])
  im1transmat[3, 1:2] <- res$reverseTransforms[[1]][1:2, 4]

  list(matrix = im1transmat, alignment_result = res)
}

align_images_matrix <- function(im1, im2, affine_mat) {

  im1trans <- EBImage::affine(im1, m = affine_mat, output.dim = dim(im1), bg.col = img_mode(im1))

  end_dim <- pmax(dim(im1trans), dim(im2))

  pad1dim <- end_dim - dim(im1trans)
  pad2dim <- end_dim - dim(im2)

  im1trans <- im1trans %>% img_pad(right = pad1dim[1], bottom = pad1dim[2], value = img_mode(.))
  im2trans <- im2 %>% img_pad(right = pad2dim[1], bottom = pad2dim[2], value = img_mode(.))

  list(im1 = im1trans, im2 = im2trans, alignment = res)
}

xx <- align_scaled_matrix(scan_info$aligned_pyr[[1]]$img[[1]],
                          scan_info$aligned_pyr[[2]]$img[[1]],
                          scale = 2)

xximg <- align_images_matrix(scan_info$clean_img[[1]],
                             scan_info$clean_img[[2]],
                             xx$matrix)

plot(rgbImage(xximg[[1]],
              xximg[[2]],
              pmax(xximg[[1]], xximg[[2]])))


im1 <- ShoeScrubR:::em_thresh(scan_info$aligned_pyr[[4]]$img[[1]])

xx <- align_scaled_matrix(
  scan_info$aligned_pyr[[4]]$img[[1]],
  scan_info$aligned_pyr[[3]]$img[[1]],
  scale = 2)

xximg <- align_images_matrix(
  scan_info$clean_img[[4]],
  scan_info$clean_img[[3]],
  xx$matrix)

plot(rgbImage(xximg[[1]],
              xximg[[2]],
              pmax(xximg[[1]], xximg[[2]])))

## Need to figure out how to handle partial images. Can restrict sample space? Align only a portion of the image? Align masks first and then align actual images based on masks? <- that one is plausible with function

# xx <- align_scaled(scan_info$aligned_pyr[[1]]$img[[1]], scan_info$aligned_pyr[[2]]$img[[1]], scale = 2)
# plot(rgbImage(scan_info$img[[1]],
#               EBImage::affine(scan_info$img[[1]],
#                               rbind(xx$reverseTransforms[[1]][1:2,1:2], t(xx$reverseTransforms[[1]][1:2,4]))),
#               scan_info$img[[1]]))
#
# im1trans <- EBImage::affine(scan_info$align[[1]]$img,
#                             rbind(t(xx$reverseTransforms[[1]][1:2,1:2]), xx$reverseTransforms[[1]][1:2,4]),
#                             output.dim = dim(scan_info$align[[1]]$img), bg.col = img_mode(scan_info$img[[1]]))
# end_dim <- pmax(dim(im1trans), dim(scan_info$align[[2]]$img))
# pad1dim <- end_dim - dim(im1trans)
# im1trans <- im1trans %>%
#   img_pad(right = pad1dim[1], bottom = pad1dim[2], value = img_mode(.))
#
# pad2dim <- end_dim - dim(scan_info$align[[2]]$img)
# im2trans <- scan_info$align[[2]]$img %>%
#   img_pad(right = pad2dim[1], bottom = pad2dim[2], value = img_mode(.))
