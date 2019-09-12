library(ShoeScrubR)
library(EBImage)
full_imglist <- list.files("/lss/research/csafe-shoeprints/ShoeImagingPermanent/",
                           pattern = "0\\d{5}[RL]_\\d{8}_5_1_1", full.names = T)

dir <- tempdir()

file.copy(full_imglist, file.path(dir, basename(full_imglist)))
imglist <- list.files(dir, "\\..*", full.names = T)

img <- EBImage::readImage(file.path(dir, "001351R_20171211_5_1_1_csafe_jekruse.tif"))
orig_img <- EBImage::readImage(file.path(dir, "040639L_20180307_5_1_1_boekhoff_pashek_jekruse.tif"))
# img <- readImage(imglist[grepl("089", imglist)][1])
orig_mask <- shoe_mask("Nike", 8, "R", ppi = 300)

actual_img <- img %>% channel("luminance")

search_img <- actual_img %>%
  filter2(makeBrush(5, "gaussian")) %>%
  magrittr::subtract(1, .) %>%
  normalize()
# search_img <- search_img > .25 # quick and dirty threshold...


max_dims <- pmax(dim(orig_mask), dim(actual_img))

mask <- auto_resize_img(orig_mask, max_dims)
actual_img <- auto_resize_img(actual_img, max_dims)

theta_res <- 1
img_downsample_factor <- 10
theta_range <- c(-5, 5)

# Compute radial samples
mdim <- dim(mask)
mcenter <- round(mdim/2)
max_radius <- round(sqrt(sum((mdim - mcenter)^2)))
slopes <- unique(seq(-90, 90, theta_res) %% 180)

test_slopes <- seq(theta_range[1], theta_range[2], by = theta_res) %% 180

img_slopes <- slopes[slopes %% img_downsample_factor == 0]

points <- crossing(slope = slopes,
                   r = seq(-max_radius, max_radius, by = 1)) %>%
  mutate(theta = slope/180*pi) %>%
  mutate(row = mcenter[1] + round(r*sin(theta)),
         col = mcenter[2] + round(r*cos(theta)))

matrix_points <- points %>%
  filter(col <= mdim[2], row <= mdim[1], row > 0, col > 0)

img_points <- points %>%
  filter(slope %in% img_slopes) %>%
  mutate(col = col - mcenter[2] + round(ncol(img)/2),
         row = row - mcenter[1] + round(nrow(img)/2)) %>%
  filter(col <= ncol(img), row <= nrow(img), row > 0, col > 0)

masked_matrix <- matrix_points
masked_matrix$mask_val <- mask[as.matrix(matrix_points[,4:5])]

masked_img <- img_points
masked_img$img_val <- search_img[as.matrix(img_points[,4:5])]

# Create single line objects for each theta
mask_lines <- masked_matrix %>%
  select(mask_theta = slope, r, mask_val) %>%
  nest(-mask_theta, .key = "mask.data")

img_lines <- masked_img %>%
  select(img_theta = slope, r, img_val) %>%
  nest(-img_theta, .key = "img.data")

# Make mask w/ white at (points$x, points$y)
angle_sample_mask <- mask * 0
angle_sample_mask[as.matrix(matrix_points[,4:5])] <- 1

img_mask <- img*0
img_mask[as.matrix(img_points[,c(4, 5)])] <- 1



# Set up combinations of thetas to test

res <- crossing(dtheta = test_slopes, img_theta = img_slopes) %>%
  mutate(mask_theta = (img_theta + dtheta) %% 180) %>%
  left_join(mask_lines, by = "mask_theta") %>%
  left_join(img_lines, by = "img_theta") %>%
  mutate(dist = purrr::map2(mask.data, img.data, mask_line_distance)) %>%
  select(dtheta, img_theta, mask_theta, dist) %>%
  unnest(dist) %>%
  group_by(dtheta) %>%
  summarize(total_white_in_mask_region = sum(sum_white_in_mask),
            total_white_outside_mask = sum(sum_white_out_mask),
            capture_loss_ratio = total_white_in_mask_region/total_white_outside_mask)


arrange(res, total_white_outside_mask)
arrange(res, total_white_in_mask_region)

rot_angle <- arrange(res, desc(capture_loss_ratio))$dtheta[1]

if (rot_angle > 90) rot_angle <- rot_angle - 180


# par(mfrow = c(2, 3))
# plot(search_img)
# plot(mask)
# plot(dilate(angle_sample_mask, makeBrush(5, "box")))
#
# plot(normalize(search_img + dilate(angle_sample_mask, makeBrush(5, "box"))))
# plot(normalize(mask + dilate(angle_sample_mask, makeBrush(5, "box"))))


rotated_mask <- EBImage::rotate(mask, angle = rot_angle, output.origin = mcenter)

mask_sized_img <- resize_img(1 - search_img, dim(rotated_mask))

plot(rgbImage(red = mask_sized_img, green = 1 - rotated_mask, blue = 1 - mask))


