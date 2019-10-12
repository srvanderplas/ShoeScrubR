library(tidyverse)
library(EBImage)
library(ShoeScrubR) # SHA: fe17d486

img_output_dir <- "~/Projects/CSAFE/2019-this_is_us/images/shoes/longitudinal/"
lss_dir <- "/lss/research/csafe-shoeprints/ShoeImagingPermanent"

templates <- tibble(file = list.files("inst/templates", full.names = T)) %>%
  mutate(fname = basename(file) %>% str_remove("\\.png$")) %>%
  tidyr::extract(fname, into = c("Shoe", "Model", "Size", "Foot"),
                 regex = "(Adidas|Nike)_(.*)_([\\d\\.]{1,})[MW]_([LR])",
                 remove = T) %>%
  mutate(img = purrr::map(file, ~readImage(.) %>% channel("luminance")))

# png(filename = file.path(img_output_dir, "Templates_Equal_Size.png"), width = 300*8, height = 300*4, units = "px")
par(mfrow = c(2, 8))
templates %>%
  arrange(Shoe, Size, Foot) %>%
  magrittr::extract2("img") %>%
  purrr::map(ShoeScrubR:::auto_resize_img, final_dims = c(1600, 2600)) %>%
  purrr::map(plot)
# dev.off()


# Generate explanation of how this would work

# Setup and initial cleaning
orig_img <- EBImage::readImage(file.path(lss_dir, "040639L_20180307_5_1_1_boekhoff_pashek_jekruse.tif"))
img <- orig_img %>% channel("luminance")

inv_img <- img %>%
  filter2(makeBrush(25, "gaussian")) %>%
  normalize() %>%
  magrittr::subtract(1, .)
thresh_img <- inv_img > .15

orig_mask <- shoe_mask("Nike", 8, "R", ppi = 300) %>% ShoeScrubR:::auto_resize_img(final_dims = dim(img))

# png(filename = file.path(img_output_dir, "Templates_Setup.png"),
    # width = 300*4, height = 300*2, units = "px")
par(mfrow = c(1, 4))
plot(img)
plot(thresh_img)
plot(orig_mask)
plot(rgbImage(green = 1 - thresh_img, blue = 1 - orig_mask, red = 1 - orig_mask))
# dev.off()

# Compute center of each image (by # white pixels)
exaggerated_img <- thresh_img %>%
  opening(makeBrush(5, "disc")) %>%
  closing(makeBrush(301, "disc"))

img_center <- exaggerated_img %>%
  image_to_df() %>%
  summarize(row = round(mean(row)), col = round(mean(col))) %>%
  unlist()

mask_center <- image_to_df(orig_mask)  %>%
  summarize(row = round(mean(row)), col = round(mean(col))) %>%
  unlist()

t_dist <- mask_center - img_center
centered_mask <- translate(orig_mask, -t_dist, bg.col = 0)


# png(filename = file.path(img_output_dir, "Templates_Naive_Centering.png"),
    # width = 300*4, height = 300*2, units = "px")
par(mfrow = c(1, 4))
plot(rgbImage(green = exaggerated_img, blue = thresh_img, red = thresh_img))
points(-img_center[["row"]], img_center[["col"]], col = "blue", cex = 1.5, type = "p", pch = 16)
plot(orig_mask)
points(-mask_center[["row"]], mask_center[["col"]], col = "red", cex = 1.5, type = "p", pch = 16)
plot(rgbImage(green = 1 - orig_mask, blue = 1 - centered_mask, red = 1 - centered_mask))
plot(rgbImage(green = 1 - thresh_img, blue = 1 - centered_mask, red = 1 - centered_mask))
# dev.off()

padded_exag_img <- ShoeScrubR:::pad_to_center(exaggerated_img, img_center)
padded_img <- ShoeScrubR:::pad_to_center(thresh_img, img_center)
padded_mask <- ShoeScrubR:::pad_to_center(centered_mask, img_center)

# Angle alignment

radial_mask_init <- ShoeScrubR:::radial_mask(dim(padded_img), slopes = seq(-90, 90, 10))

radial_masks <- purrr::map(
  c(-6, -2, 0, 2, 6),
  ~ShoeScrubR:::radial_mask(dim(padded_img),
                            slopes = seq(-90, 90, 10) + .)
)


# png(filename = file.path(img_output_dir, "Templates_Radial_Mask_Single.png"),
    # width = 300*1, height = 300*2, units = "px")
purrr::map(radial_masks[3], ~plot(rgbImage(red = ., green = padded_img, blue = padded_mask)))
# dev.off()


# png(filename = file.path(img_output_dir, "Templates_Radial_Mask_Multiple.png"),
    # width = 300*5, height = 300*2, units = "px")
par(mfrow = c(1, 5))
purrr::map(radial_masks, ~plot(rgbImage(red = ., green = padded_img, blue = padded_mask)))
# dev.off()

odim <- floor(dim(padded_img)/10)
padded_img_small <- EBImage::resize(padded_img, h = odim[2], w = odim[1])
exag_img_small <- EBImage::resize(padded_exag_img, h = odim[2], w = odim[1])
padded_mask_small <- EBImage::resize(padded_mask, h = odim[2], w = odim[1])

img_df <- image_to_df(padded_img_small) %>%
  mutate(row = -row)
exag_img_df <- image_to_df(exag_img_small) %>%
  mutate(row = -row)
mask_df <- image_to_df(padded_mask_small) %>%
  mutate(row = -row)

angle_options <- seq(-15, 15, 1)
radial_masks <- purrr::map(
  angle_options,
  ~ShoeScrubR:::radial_mask_df(odim,
                            slopes = seq(-90, 90, 5) + .)
)

radial_mask_img <- ShoeScrubR:::radial_mask_df(dim(exag_img_small), slopes = seq(-90, 90, 5))

img_data <- img_df %>% inner_join(radial_mask_img) %>%
  rename(img_val = val)
exag_img_data <- exag_img_df %>% inner_join(radial_mask_img) %>% rename(img_val = val)
mask_data <- purrr::map(radial_masks, ~inner_join(., mask_df) %>% rename(mask_val = val))

dists <- purrr::map_df(mask_data, ~mask_line_distance(., exag_img_data))
dists <- dists %>%
  mutate(angle = angle_options) %>%
  arrange(desc(sum_white_in_mask))

rgbImage(padded_mask %>% EBImage::rotate(-7, output.dim = dim(padded_img), output.origin = floor(dim(padded_mask)/2)), padded_exag_img, 1-padded_img) %>% plot()
rgbImage(padded_mask %>% EBImage::rotate(15, output.dim = dim(padded_img), output.origin = floor(dim(padded_mask)/2)), padded_img, 0*padded_img) %>% plot()

dists %>%
  gather(key = "measure", value = "value", -angle) %>%
  ggplot(aes(x = angle, y = value, color = measure)) + geom_line() + facet_grid(measure~., space = "free", scales = "free")

rgbImage(padded_mask, padded_img, 0*padded_img) %>% plot()
