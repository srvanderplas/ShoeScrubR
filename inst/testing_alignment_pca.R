library(tidyverse)
library(EBImage)
library(ShoeScrubR)

test_img <- matrix(0, nrow = 49, ncol = 49) %>% as.Image()
test_img[24:26,] <- 1
test_img <- pad(test_img, 10, 10, 10, 10)


get_angle <- function(rot) {
  if (hasName(rot, "rotation")) {
    rot <- rot$rotation
  }
  stopifnot(all.equal(dim(rot), c(2,2)))

  idx <- which(abs(rot) > 0.1)
  est1 <- acos(rot[2,2]) %>% `*`(180/pi)
  est2 <- acos(rot[1,2]) %>% `*`(180/pi)
  est3 <- asin(rot[2,1]) %>% `*`(180/pi)
  est4 <- acos(rot[1,1]) %>% `*`(180/pi)
  est5 <- atan2(rot[2,2], rot[2,1]) %>% `*`(180/pi)
  est6 <- atan2(rot[1,2], rot[1,1]) %>% `*`(180/pi)
  tibble(e1 = est1, e2 = est2, e3 = est3, e4 = est4, e5 = est5, e6 = est6)
}

rotations <- tibble(img = list(test_img), angle = list(angle = seq(0, 180-5, 5))) %>%
  unnest(angle, .drop = F) %>%
  mutate(rot_img = purrr::map2(img, angle, rotate, output.dim = dim(test_img), bg.col = 0),
         pca = purrr::map(rot_img, align_prcomp, weighted = F),
         est_angle = purrr::map(pca, get_angle)) %>%
  unnest(est_angle, .drop = F) %>%
  mutate(fix_rot = purrr::map2(rot_img, -e4, rotate, output.dim = dim(test_img), bg.col = 0)) %>%
  nest(angle1 = c(e1:e6)) %>%
  mutate(pca2 = purrr::map(fix_rot, align_prcomp, weighted = F),
         est_angle = purrr::map(pca2, get_angle)) %>%
  unnest(est_angle) %>%
  mutate(fix_rot2 = purrr::map2(fix_rot, -e4, rotate, output.dim = dim(test_img), bg.col = 0)) %>%
  nest(angle2 = c(e1:e6)) %>%
  mutate(pca3 = purrr::map(fix_rot2, align_prcomp, weighted = F),
         est_angle = purrr::map(pca3, get_angle)) %>%
  unnest(est_angle) %>%
  mutate(fix_rot3 = purrr::map2(fix_rot2, -e4, rotate, output.dim = dim(test_img), bg.col = 0))

par(mfrow = c(6, 12))
map(rotations$rot_img, plot)
map(rotations$pca_out, plot)

