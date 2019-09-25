library(tidyverse)
library(EBImage)
library(ShoeScrubR)

img_output_dir <- "~/Projects/CSAFE/2019-this_is_us/images/shoes/longitudinal/"
lss_dir <- "/lss/research/csafe-shoeprints/ShoeImagingPermanent"

# Setup and initial cleaning
orig_img <- EBImage::readImage(file.path(lss_dir, "040639L_20180307_5_1_1_boekhoff_pashek_jekruse.tif"))
img <- orig_img %>% channel("luminance")

inv_img <- img %>%
  filter2(makeBrush(25, "gaussian")) %>%
  normalize() %>%
  magrittr::subtract(1, .)
thresh_img <- inv_img > .15

orig_mask <- shoe_mask("Nike", 8, "R", ppi = 300) %>% auto_resize_img(final_dims = dim(img)) %>%
  closing(makeBrush(31, "disc"))

tmp <- rough_align(img, orig_mask)

rgbImage(1 - tmp$img[[1]], tmp$img[[2]], tmp$img[[3]]) %>% plot()



# For a bunch of images...
full_imglist <- list.files("/lss/research/csafe-shoeprints/ShoeImagingPermanent/",
                           pattern = "00\\d{4}[RL]_\\d{8}_5_1_1", full.names = T)
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
  group_by(ShoeID, Shoe_foot) %>%
  arrange(desc(date)) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  mutate(
    img = purrr::map(file, EBImage::readImage, all = F),
    img = purrr::map(img, EBImage::channel, "luminance"),
    im_dim = purrr::map(img, dim)
  ) %>%
  left_join(select(shoe_info, ShoeID, Brand, Size, Shoe_foot)) %>%
  group_by(Shoe_foot, Brand) %>%
  sample_n(2) %>%
  ungroup() %>%
  left_join(select(shoe_info, ShoeID, Shoe_foot, Mask_foot, mask)) %>%
  mutate(mask = purrr::pmap(list(mask, im_dim), ~auto_resize_img(..1, ..2, 0)))

scan_align <- scan_info %>%
  mutate(align = purrr::map2(img, mask, rough_align))

max_dims <- map_df(scan_align$align, ~dim(.$img[[1]]) %>% t() %>% as.data.frame() %>% set_names(c("col", "row"))) %>% summarize_each(max) %>% as.numeric()

scan_align <- scan_align %>%
  mutate(align_resize = purrr::map(align, function(df) {
    mutate(df, im_mode = purrr::map_dbl(img, img_mode),
           img = purrr::map2(img, im_mode, ~img_pad_to_size(.x, value = .y, size = max_dims)))
  }))

plot_align <- function(df) {
  rgbImage(df$img[[1]], 1 - df$img[[2]], (1 - df$img[[4]])) %>% plot()
}

# make into image for show and tell
cols <- 8
plot_dims <- c(pmax(ceiling(nrow(scan_align)/cols), 1), pmin(nrow(scan_align), cols))
png(filename = file.path(img_output_dir, "PCA_Rotate_and_Center_Shift_Before.png"), width = 300*plot_dims[2], height = 300*2*plot_dims[1], units = "px")
par(mfrow = plot_dims)
purrr::pwalk(list(scan_align$img, scan_align$mask), ~rgbImage(..1, (1 - ..2), 1 - (..1)*(1 - ..2)) %>% plot())
dev.off()

# make into image for show and tell
plot_dims <- c(pmax(ceiling(nrow(scan_align)/5), 1), pmin(nrow(scan_align), 5))
png(filename = file.path(img_output_dir, "PCA_Rotate_and_Center_Shift_After.png"), width = 300*plot_dims[2], height = 300*2*plot_dims[1], units = "px")
par(mfrow = plot_dims)
purrr::walk(scan_align$align, plot_align)
dev.off()
