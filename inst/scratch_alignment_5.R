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

orig_mask <- shoe_mask("Nike", 8, "R", ppi = 300) %>% auto_resize_img(final_dims = dim(img))

# Compute center of each image (by # white pixels)
exaggerated_img <- img %>%
  filter2(makeBrush(125, "gaussian")) %>%
  normalize() %>%
  magrittr::subtract(1, .) %>%
  (function(.) . > .125) %>%
  opening(makeBrush(7, "disc")) %>%
  closing(makeBrush(301, "disc"))

img_center <- exaggerated_img %>%
  image_to_df() %>%
  summarize(row = round(mean(row, weight = value)), col = round(mean(col, weight = value))) %>%
  unlist()

mask_center <- image_to_df(orig_mask)  %>%
  summarize(row = round(mean(row)), col = round(mean(col))) %>%
  unlist()

t_dist <- mask_center - img_center
centered_mask <- translate(orig_mask, -t_dist, bg.col = 0)

padded_exag_img <- pad_to_center(exaggerated_img, img_center)
padded_img <- pad_to_center(thresh_img, img_center)
padded_mask <- pad_to_center(centered_mask, img_center)

im_list <- list(orig = img, inv = inv_img, thresh = thresh_img, exaggerated = exaggerated_img, orig_mask = orig_mask, centered_mask = centered_mask)
im_px_mode <- c(1, 0, 0, 0, 0, 0)
im_list_pad <- purrr::map2(im_list, im_px_mode, ~pad_to_center(img = .x, value = .y, center = img_center))

pyramid <- img_pyramid(im_list_pad, scale = c(1, 4, 8, 16, 32)) %>%
  mutate(df = purrr::map(img, image_to_df)) %>%
  mutate(rot_angle_prcomp = purrr::map_dbl(df, align_prcomp))

mask_imgs <- pyramid %>%
  filter(img_name %in% c("centered_mask", "thresh")) %>%
  mutate(rot_img = purrr::map2(img, rot_angle_prcomp, ~EBImage::rotate(.x, angle = .y, output.dim = dim(.x)))) %>%
  mutate(rot_center = purrr::map(rot_img, . %>%
                                   image_to_df() %>%
                                   summarize(row = round(mean(row)), col = round(mean(col))) %>%
                                   unlist()))

masking_img <- mask_imgs %>%
  select(scale, img_name, rot_img) %>%
  tidyr::spread(key = img_name, value = rot_img) %>%
  left_join(
    select(mask_imgs, scale, img_name, rot_center) %>%
      mutate(img_name = paste0(img_name, "_rot_center")) %>%
      tidyr::spread(key = img_name, value = rot_center) %>%
      mutate(mask_shift = purrr::map2(centered_mask_rot_center, thresh_rot_center, ~.x - .y)) %>%
      select(-matches("rot_center"))
  ) %>%
  mutate(recentered_mask = purrr::map2(centered_mask, mask_shift, ~translate(.x, -.y, bg.col = 0)))


masking_img_2 <- masking_img %>%
  select(scale, mask = recentered_mask, img = thresh) %>%
  mutate(by = round(40/(scale)))




get_offsets <- function(dims, by = 2, range = list(width = c(.3, .7), height = c(.3, .7))) {
  width_px_range <- round(dims[1]*range$width)
  height_px_range <- round(dims[2]*range$height)

  width_px <- seq(width_px_range[1], width_px_range[2], by = by) - floor(dims[1]/2)
  height_px <- seq(height_px_range[1], height_px_range[2], by = by) - floor(dims[2]/2)

  tidyr::crossing(row = height_px, col = width_px) %>%
    mutate(id = 1:n()) %>%
    nest(row:col, .key = "offset")
}

tmp <- masking_img_2 %>%
  filter(scale > 5) %>%
  mutate(offsets = map2(img, by, ~get_offsets(dim(.x), by = .y))) %>%
  unnest(offsets) %>%
  mutate(t_mask = map2(mask, offset, ~ translate(.x, -as.numeric(unlist(.y)), bg.col = 0))) %>%
  mutate(overlap = purrr::map2_dbl(img,t_mask, ~(sum(.x*.y)/sum(.x)))) %>%
  group_by(scale) %>%
  arrange(desc(overlap)) %>%
  unnest(offset) %>%
  filter(row_number() <= 1) %>%
  arrange(desc(scale))

tmp <- tmp %>%
  arrange(desc(scale)) %>%
  mutate(color_img = purrr::map2(img, t_mask, ~rgbImage(red = .x*.y, green = .x, blue = .y))) %>%
  mutate(label = purrr::pmap_chr(list(overlap, row, col, scale), ~sprintf("Scale:1/%d\n%.1f%% overlap,\nOS = (%d, %d)",..4,  ..1*100, round(..2)*..4, round(..3)*..4)))

# png(filename = file.path(img_output_dir, "ShiftPyramids.png"), width = 300*nrow(tmp), height = 300*2*1, units = "px")
par(mfrow = c(1, nrow(tmp)))
purrr::map2(tmp$color_img, tmp$label, ~{
  plot(.x)
  text(0, 0, labels = .y, col = "white", adj = c(-.1, 1.2))
})
# dev.off()

get_angles <- function(dims, by = 1, range = c(-5, 5), centered = T) {
  c_rot_div <- ifelse(centered, 2, 1)
  c_rot <- floor(dims/c_rot_div)
  tibble(angle = seq(range[1], range[2], by = by),
         center = list(c_rot)) %>%
    mutate(id = 1:n())
}

tmp2 <- tmp %>%
  ungroup() %>%
  select(scale, img, mask = t_mask, t_overlap = overlap) %>%
  mutate(angles = map(img, ~get_angles(dim(.)))) %>%
  unnest(angles, .drop = F)

tmp <- tmp2 %>%
  mutate(odim = map(img, dim)) %>%
  mutate(r_mask = pmap(list(x = mask, angle = angle, output.dim = odim, output.origin = center) , rotate)) %>%
  mutate(overlap = purrr::map2_dbl(img, r_mask, ~(sum(.x*.y)/sum(.x)))) %>%
  mutate(color_img = purrr::map2(img, r_mask, ~rgbImage(red = .x*.y, green = .x, blue = .y))) %>%
  group_by(scale) %>%
  arrange(desc(overlap)) %>%
  filter(row_number() <= 1) %>%
  ungroup() %>%
  arrange(scale)

par(mfrow = c(3, 11))
map(tmp$color_img, plot)
