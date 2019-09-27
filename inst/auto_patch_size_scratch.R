library(tidyverse)
library(EBImage)
library(ShoeScrubR)

img_output_dir <- "~/Projects/CSAFE/2019-this_is_us/images/shoes/longitudinal/"
lss_dir <- "/lss/research/csafe-shoeprints/ShoeImagingPermanent"

# Setup and initial cleaning
orig_img <- EBImage::readImage(file.path(lss_dir, "040639L_20180307_5_1_1_boekhoff_pashek_jekruse.tif"))
img <- orig_img %>% channel("luminance")


pyr <- img_pyramid(img, scale = c(1, 3, 5, 7))

pyr$img[[4]]

pyr$data <- purrr::map(pyr$img, as.numeric)

pyr_df <- pyr %>%
  select(scale, data) %>%
  unnest(data)
ggplot(data = pyr_df, aes(x = data, color = factor(scale))) + geom_density()

filter_down <- function(img, scale = 3) {
  offset <- floor(scale/2) + 1 # center of x and y when pulling tiles out
  fimg <- filter2(img, matrix(1/scale^2, nrow = scale, ncol = scale))
  coords <- expand.grid(seq(offset, dim(fimg)[1], by = scale),
                        seq(offset, dim(fimg)[2], by = scale))
  coords$value <- fimg[cbind(coords$Var1, coords$Var2)]
  coords
}


filter_pyr <- tibble(img = list(img),
                     thresh = clean_initial_img(img, threshold_val = 0),
                     filter_d = c(5, 15, 25, 35, 45)) %>%
  mutate(fd = purrr::map2(img, filter_d, filter_down))

ggplot(data = unnest(filter_pyr, fd), aes(x = value, color = factor(filter_d))) +
  geom_density() +
  scale_y_continuous(limits = c(0, 100))


blockwise_stats <- function(img, size) {
  img <- normalize(img)
  dims <- dim(img)
  img_coords <- tibble(idx = 1:length(img), d1 = floor((floor((idx - 1)/dims[2]))/size) + 1, d2 = floor(((idx - 1) %% dims[1])/size) + 1) %>%
    group_by(d1, d2) %>%
    summarize(s = sum(img[idx]), sd = sd(img[idx]), n = n())
}


blockwise <-tibble(img = list(img),
                   thresh = clean_initial_img(img, threshold_val = 0),
                   filter_d = c(5, 15, 25, 35, 45)) %>%
  mutate(bs = purrr::map2(img, filter_d, blockwise_stats))

tmp <- blockwise %>% unnest(bs) %>%
  mutate(p = s/n)

tmp %>%
  ggplot(aes(x = p*(1-p), fill = factor(filter_d))) +
  geom_density(alpha = .5)


tmp %>% group_by(filter_d) %>% summarize(x = mean(p*(1-p) < .025))



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


blockwise <- crossing(scan_info,
                      filter_d = c(25, 35, 45, 75)) %>%
  mutate(bs = purrr::map2(img, filter_d, blockwise_stats))

tmp <- blockwise %>% unnest(bs) %>%
  mutate(p = s/n)

tmp %>%
  ggplot(aes(x = p*(1-p), fill = factor(filter_d), group = interaction(filter_d, file))) +
  geom_density(alpha = .5) +
  facet_wrap(Brand~ShoeID)


tmp %>% group_by(filter_d) %>% summarize(x = mean(p*(1-p) < .025))

tmp %>%
  ggplot(aes(x =sd, fill = factor(filter_d), group = interaction(filter_d, file))) +
  geom_density(alpha = .5) +
  facet_wrap(Brand~ShoeID)
