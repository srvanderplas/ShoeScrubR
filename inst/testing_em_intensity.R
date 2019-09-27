library(tidyverse)
library(EBImage)
library(ShoeScrubR)

img_output_dir <- "~/Projects/CSAFE/2019-this_is_us/images/shoes/longitudinal/"
lss_dir <- "/lss/research/csafe-shoeprints/ShoeImagingPermanent"

# For a bunch of images...
full_imglist <- list.files("/lss/research/csafe-shoeprints/ShoeImagingPermanent/",
                           pattern = "0[01]\\d{4}[RL]_\\d{8}_5_1_1", full.names = T)
dir <- tempdir()

file.copy(full_imglist, file.path(dir, basename(full_imglist)))
imglist <- file.path(dir, basename(full_imglist))

shoe_info <- read.csv("~/Projects/CSAFE/2018_Longitudinal_Shoe_Project/Clean_Data/shoe-info.csv") %>%
  filter(ShoeID %in% as.numeric(str_sub(basename(imglist), 1, 3))) %>%
  select(ShoeID, Brand, Size) %>%
  mutate(Size = str_remove(Size, "[ MW]") %>% parse_number()) %>%
  crossing(tibble(Mask_foot = c("R", "L"), Shoe_foot = c("L", "R"))) %>%
  mutate(mask = purrr::pmap(list(Brand, Size, Mask_foot, ppi = 300), shoe_mask))

scan_info <- tibble(
  file = imglist,
  ShoeID = str_extract(basename(file), "^\\d{3}") %>% parse_integer(),
  Shoe_foot = str_extract(basename(file), "\\d{6}[RL]") %>% str_remove_all("\\d"),
  date = str_extract(basename(file), "\\d{8}") %>% parse_date(format = "%Y%m%d")
) %>%
  left_join(select(shoe_info, ShoeID, Brand, Size, Shoe_foot)) %>%
  group_by(Shoe_foot, Brand) %>%
  sample_n(5) %>%
  ungroup() %>%
  group_by(ShoeID, Shoe_foot) %>%
  arrange(desc(date)) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(
    img = purrr::map(file, EBImage::readImage, all = F),
    img = purrr::map(img, EBImage::channel, "luminance"),
    im_dim = purrr::map(img, dim)
  )


img <- scan_info$img[[2]]
img10 <- img_resize(img, w = floor(dim(img)[1]/10), h = floor(dim(img)[2]/10))
img20 <- img_resize(img, w = floor(dim(img)[1]/20), h = floor(dim(img)[2]/20))

iem20 <- mixtools::normalmixEM(as.numeric(img20), k = 2)
iem20$posterior
img20_post <- img20
img20_post[1:length(img20)] <- iem20$posterior[,1]
plot(img20_post)

iem10 <- mixtools::normalmixEM(as.numeric(img10), k = 2)
img10_post <- img10
img10_post[1:length(img10)] <- iem10$posterior[,1]
plot(img10_post)


dnorm(seq(0, 1, .01), mean = iem10$mu[1], sd = iem10$sigma[1])
mapply(function(a, b) dnorm(seq(0, 1, .001), mean = a, sd = b), iem10$mu, iem10$sigma) %>%
  apply(1, which.max) %>%
  (function(x) seq(0, 1, .001)[min(which(x == 2))])
