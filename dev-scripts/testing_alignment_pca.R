library(tidyverse)
library(EBImage)
library(ShoeScrubR)

test_img <- matrix(0, nrow = 49, ncol = 49) %>% as.Image()
test_img[24:26,] <- 1
test_img <- pad(test_img, 10, 10, 10, 10)



rotations <- tibble(img = list(test_img), angle = list(angle = seq(0, 180 - 5, 5))) %>%
  unnest(angle, .drop = F) %>%
  mutate(rot_img = purrr::map2(img, angle, rotate, output.dim = dim(test_img), bg.col = 0),
         est_angle = purrr::map(rot_img, align_prcomp, weighted = T)) %>%
  mutate(fix_rot = purrr::map2(rot_img, est_angle, rotate, output.dim = dim(test_img), bg.col = 0))

par(mfrow = c(6, 12))
map(rotations$rot_img, plot)
map(rotations$fix_rot, plot)

