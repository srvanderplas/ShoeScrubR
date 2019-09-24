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
  crossing(Mask_foot = c("R", "L")) %>%
  mutate(Shoe_foot = case_when(Mask_foot == "R" ~"L", Mask_foot == "L" ~ "R")) %>%
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
  sample_n(10) %>%
  mutate(
  img = purrr::map(file, EBImage::readImage, all = F)
)

scan_info <- left_join(scan_info, shoe_info)
scan_align <- scan_info %>%
  mutate(img = purrr::map(img, EBImage::channel, "luminance")) %>%
  mutate(img_dim = purrr::map(img, dim)) %>%
  mutate(mask = purrr::map2(mask, img_dim, auto_resize_img)) %>%
  mutate(align = purrr::map2(img, mask, rough_align))

min_dims <- map_df(scan_align$align, ~dim(.$img[[1]]) %>% t() %>% as.data.frame() %>% set_names(c("col", "row"))) %>% summarize_each(min) %>% as.numeric()

scan_align <- scan_align %>%
  mutate(align_resize = purrr::map(align, function(df) {
    mutate(df, img = img_crop(img, dim = min_dims))
  }))

plot_align <- function(df) {
  rgbImage(1-df$img[[1]], df$img[[2]], df$img[[3]]) %>% plot()
}

# make into image for show and tell
png(filename = file.path(img_output_dir, "PCA_Rotate_and_Center_Shift.png"), width = 300*pmax(nrow(scan_align), 5), height = 300*2*pmax(ceiling(nrow(scan_align)/5), 1), units = "px")
par(mfrow = c(1, 5))
purrr::walk(scan_align$align_resize, plot_align)
dev.off()
