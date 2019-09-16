#' Get a specific shoe mask
#'
#' @param brand Adidas or Nike
#' @param size Numerical size: 7, 7.5, 10, 10.5 for Adidas, 8, 8.5, 10, 10.5
#'        for Nike
#' @param foot R or L, determined from the image (e.g. if the image is a film
#'        print, this would be the opposite foot from the image name because
#'        film prints are reversed)
#' @param ppi pixels per inch (linearly) for the returned mask image
#' @export
#' @examples
#' shoe_mask("Nike", 10.5, "R") %>% plot()
#' plot(normalize(shoe_mask("Nike", 10, "R") + shoe_mask("Nike", 10.5, "R")))
shoe_mask <- function(brand, size, foot, ppi = 300) {
  stopifnot(brand %in% c("Adidas", "Nike"))
  if (brand == "Adidas") {
    stopifnot(size %in% c(7, 7.5, 10, 10.5))
  } else {
    stopifnot(size %in% c(8, 8.5, 10, 10.5))
  }
  stopifnot(foot %in% c("R", "L"))

  brand_model <- gsub("Adidas", "Adidas_Seeley", brand) %>%
    gsub("Nike", "Nike_Winflow", .)

  size_x <- ifelse(size < 9, paste0(size, "W"), paste0(size, "M"))

  filename <- paste0(paste(brand_model, size_x, foot, sep = "_"), ".png")
  full_filename <- system.file("templates", filename, package = "ShoeScrubR")

  stopifnot(file.exists(full_filename))
  im <- EBImage::readImage(full_filename)
  scale <- ppi/200
  output_dim <- floor(dim(im)*scale)
  EBImage::resize(im, w = output_dim[1], h = output_dim[2])
}


#' Optimize the mask position on the image by rotation
#'
#' @param img image
#' @param mask mask
#' @param theta_parms list of parameters determining angles used for radial sampling of image and mask
#' @param offset_parms list of parameters determining mask offsets used for radial sampling
#' @param debug get summary results (FALSE) or full results (TRUE). Full results
#'        are more useful when optimizing over mask_offset or another parameter
#'        in addition to the angle of the mask image
#' @export
#' @examples
#' \dontrun{
#' orig_img <- EBImage::readImage(file.path("/lss/research/csafe-shoeprints/ShoeImagingPermanent/",
#' "040639L_20180307_5_1_1_boekhoff_pashek_jekruse.tif"))
#' orig_mask <- shoe_mask("Nike", 8, "R", ppi = 300)
#' img <- orig_img %>% channel("luminance") %>%
#'   filter2(makeBrush(5, "gaussian")) %>%
#'   magrittr::subtract(1, .) %>%
#'   normalize()
#' max_dims <- pmax(dim(orig_mask), dim(orig_img))
#' mask <- auto_resize_img(orig_mask, img)
#' img <- auto_resize_img(img, max_dims)
#'
#' res <- mask_image_align(img, mask, theta_parms = list(theta_res = 1, theta_range = c(-10, 10), img_res = 10), offset_parms = list(offset_radius = c(500, 500), offset_res = c(50, 50)), debug = F)
#' }
mask_image_align <- function(img, mask,
                             theta_parms = list(),
                             offset_parms = list(), debug = F) {

  if (length(theta_parms) == 0) {
    theta_res <- 2.5
  } else {
    theta_res <- theta_parms[["theta_res"]]
  }

  offsets <- do.call("expand_offsets", offset_parms)

  angles <- do.call("expand_angles", theta_parms)

  img_thetas <- unique(angles$img_theta)
  img_lines <- radial_mask_img(img, img_thetas) %>%
    dplyr::rename(img_theta = theta, img.data = data)

  mask_thetas <- angles %>% select(-img_theta) %>% nest(-dtheta) %>%
    dplyr::mutate(slopes = purrr::map(data, ~as.numeric(unlist(.)))) %>%
    dplyr::select(-data)

  max_radius <- pmin(get_max_radius(img),
                     get_max_radius(mask))

  mask_combos_init <- tidyr::crossing(cent_offset = offsets, mask_thetas) %>%
    dplyr::mutate(combo_row = 1:n())

  if (nrow(mask_combos_init) > 50000) warning("More than 50000 combinations of factors; this may take a while")

  mask_combos_lines <- purrr::map2(mask_combos_init$cent_offset,
                                   mask_combos_init$slopes,
                                   ~radial_mask_img(mask, slopes = .y,
                                                    cent_offset = .x) %>%
                                     dplyr::rename(mask_theta = theta, mask.data = data))

  mask_combos_lines <- purrr::map(mask_combos_lines, ~bind_cols(., img_lines))

  mask_combos <- tibble::tibble(combo_row = 1:length(mask_combos_lines), data = mask_combos_lines) %>%
    tidyr::unnest(data) %>%
    left_join(select(mask_combos_init, combo_row, cent_offset, dtheta))


  full_res <- mask_combos %>%
    dplyr::mutate(dist = purrr::map2(mask.data, img.data, mask_line_distance))

  full_res <- full_res %>%
    dplyr::select(combo_row, cent_offset, dtheta, img_theta, mask_theta, dist) %>%
    tidyr::unnest(dist) %>%
    dplyr::group_by(dtheta, combo_row) %>%
    dplyr::summarize(
      total_white_in_mask_region = sum(sum_white_in_mask),
      total_white_outside_mask = sum(sum_white_out_mask),
      capture_rate = total_white_in_mask_region/(total_white_in_mask_region + total_white_outside_mask)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(capture_rate)) %>%
    dplyr::left_join(dplyr::select(mask_combos_init, combo_row, cent_offset))

  if (debug) return(full_res)

  full_res %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::select(dtheta, cent_offset) %>%
    dplyr::mutate(row_offset = purrr::map_dbl(cent_offset, 1),
                  col_offset = purrr::map_dbl(cent_offset, 2)) %>%
    dplyr::select(-cent_offset)
}


#' @param mask_offset number of pixels (vector of length 2) to offset the
#'        center of the mask by
expand_offsets <- function(offset_radius = c(400, 400), offset_res = c(50, 50), offset_quadrant = 1:4) {
  if (length(offset_radius == 1)) {
    offset_radius <- rep(offset_radius, 2)
  } else if (length(offset_radius != 2)) {
    stop("Cannot understand offset_radius argument: should be length 2 or 1")
  }

  offset_rows <- seq(-offset_radius[1], offset_radius[1], by = offset_res[1])
  offset_cols <- seq(-offset_radius[2], offset_radius[2], by = offset_res[2])

  offsets <- tidyr::crossing(row_offset = offset_rows, col_offset = offset_cols)

  if (length(offset_quadrant) < 4) {
    if (!1 %in% offset_quadrant) {
      offsets <- filter(offsets, !(row_offset > 0 & col_offset > 0))
    }
    if (!2 %in% offset_quadrant) {
      offsets <- filter(offsets, !(row_offset > 0 & col_offset < 0))
    }
    if (!3 %in% offset_quadrant) {
      offsets <- filter(offsets, !(row_offset < 0 & col_offset < 0))
    }
    if (!4 %in% offset_quadrant) {
      offsets <- filter(offsets, !(row_offset < 0 & col_offset > 0))
    }
  }

  offsets %>%
    dplyr::mutate(vec = purrr::map2(row_offset, col_offset, ~c(.x, .y))) %>%
    `[[`("vec")
}

#' Get all angle combinations for mask, image, and offset
#'
#' @param theta_res resolution of optimization angle, in degrees
#' @param theta_range range of potential rotation angles to test (vector of
#'        length 2)
#' @param img_res resolution of theta to be used in the image radial sampling,
#'        should be divisible by theta_res.
#' @export
expand_angles <- function(theta_res = 2.5, theta_range = c(-5, 5), img_res = 10*theta_res) {

  # Need to test that values make sense, e.g.
  # theta_res << 90,
  # theta_res < diff(theta_range),
  # img_res is reasonable

  test_slopes <- seq(theta_range[1], theta_range[2], by = theta_res)
  img_slopes <- seq(-90, 90 - img_res, by = img_res)

  tidyr::crossing(dtheta = test_slopes, img_theta = img_slopes) %>%
    dplyr::mutate(mask_theta = (img_theta + dtheta))
}


get_max_radius <- function(x, ...) UseMethod("get_max_radius", x)

get_max_radius.numeric <- function(x, cent_offset = NULL, ...) {
  center <- floor(x/2)
  if (!is.null(cent_offset)) center <- center + cent_offset

  return(round(pmax(sqrt(sum((x - center)^2)), sqrt(sum((1 - center)^2)))))
}

get_max_radius.integer <- function(x, cent_offset = NULL, ...) {
  center <- floor(x/2)
  if (!is.null(cent_offset)) center <- center + cent_offset

  return(round(pmax(sqrt(sum((x - center)^2)), sqrt(sum((1 - center)^2)))))
}

get_max_radius.Image <- function(x, ...) {
  dd <- dim(x)
  return(get_max_radius(dd, ...))
}

#' Apply mask data frame to an image
radial_mask_img <- function(im, nest = T, ...) {
  xx <- radial_mask_df(dims = dim(im), ...)
  xx <- dplyr::mutate(xx, value = imageData(im)[cbind(xx$row, xx$col)])
  if (nest) {
    dplyr::select(xx, theta = slope, r, val = value) %>%
      tidyr::nest(-theta, .key = "data")
  } else {
    xx
  }

}

#' Make mask image
radial_mask <- function(dims, ...) {
  tmp <- radial_mask_df(dims = dims, ...)
  im <- matrix(0, nrow = dims[1], ncol = dims[2])
  im[cbind(tmp$row, tmp$col)] <- 1
  im <- Image(im)
  im
}

#' Helper function - create data frame of mask coordinates given center, dims, slopes
radial_mask_df <- function(dims, slopes = seq(-90, 90, 5),
                           expandBrush = EBImage::makeBrush(5, "disc"),
                           cent_offset = NULL) {

  # Center arg allows for some "jittering" of what the image center really is --
  # will make optimization faster?
  center <- floor(dims/2)
  if (!is.null(cent_offset)) center <- center + cent_offset

  max_radius <- get_max_radius(dims, cent_offset = cent_offset)

  tmp <- tidyr::crossing(slope = slopes,
                  r = seq(-max_radius, max_radius, by = 1)) %>%
    dplyr::mutate(theta = slope/180*pi) %>%
    dplyr::mutate(row = center[1] + round(r*sin(theta)),
                  col = center[2] + round(r*cos(theta))) %>%
    unique() %>%
    dplyr::filter(col <= dims[2], row <= dims[1], row > 0, col > 0)

  expandBrush_df <- image_to_df(expandBrush) %>%
    dplyr::mutate(row = row - median(row),
                  col = col - median(col))

  tmp2 <- tmp %>%
    dplyr::mutate(expanded = purrr::map2(row, col,
                                  function(x, y) dplyr::select(expandBrush_df, row, col) %>%
                                    dplyr::mutate(r.index = 1:n(),
                                                  row = row + x,
                                                  col = col + y))) %>%
    dplyr::select(-row, -col) %>%
    tidyr::unnest(expanded) %>%
    dplyr::filter(col <= dims[2], row <= dims[1], row > 0, col > 0)

  tmp2
}

#' Function to get distance from line datasets
mask_line_distance <- function(mask_data, img_data) {
  inner_join(mask_data, img_data, by = c("r", "r.index")) %>%
    summarize(sum_white_in_mask = sum(img_val*mask_val, na.rm = T),
              sum_white_out_mask = sum(img_val*(1-mask_val), na.rm = T))
}

plot_rotated_mask_img <- function(img, mask, align_parms) {
  angle <- align_parms$dtheta
  offset <- c(align_parms$row_offset, align_parms$col_offset)

  rotated_mask <- EBImage::rotate(mask, angle = angle, output.origin = floor(dim(mask)/2))
  translated_mask <- EBImage::translate(rotated_mask, -offset, bg.col = 0)

  mask_sized_img <- auto_resize_img(1 - img, dim(translated_mask))
  plot(rgbImage(red = mask_sized_img, green = 1 - translated_mask, blue = 1 - mask))
}

auto_resize_img <- function(x, final_dims) {
  img_dims <- dim(x)
  diffs <- img_dims - final_dims

  new_img <- x

  # Rows first
  if (diffs[1] > 0) {
    # Need to crop
    profile_row <- rowMeans(new_img)
    # Get average of the middle 2/3 of the image
    middle <- mean(profile_row[round(.16*length(profile_row)):round(.84*length(profile_row))])
    # Get maximum offset size
    max_offset <- abs(diffs[1])
    # Calculate mean value of the cropped image given an offset of x
    offset_vals <- sapply(1:max_offset, function(x) mean(profile_row[x:(x + final_dims[1] - 1)]))
    # Find which offset is closest to the middle 2/3 of the image
    best_offset <- which.min((middle - offset_vals)^2)
    # Crop using the best offset
    new_img <- new_img[best_offset:(best_offset + final_dims[1] - 1), ]
  } else if (diffs[1] < 0) {
    # Need to pad
    # Get pad value - median pixel from 10 rows on each edge of the image
    value <- median(new_img[c(1:10,(img_dims[1] - 10):img_dims[1]),])
    # Create new image
    temp_img <- matrix(value, nrow = final_dims[1], ncol = ncol(new_img)) %>% as.Image()
    # Determine offset:
    offset <- floor(abs(diffs[1]/2))
    # Add in old values
    temp_img[offset:(img_dims[1] + offset - 1), ] <- new_img
    # Save bigger image as the new_img object
    new_img <- temp_img
  }

  # Cols next
  if (diffs[2] > 0) {
    # Need to crop
    profile_col <- colMeans(new_img)
    # Get average of the middle 2/3 of the image
    middle <- mean(profile_col[round(.16*length(profile_col)):round(.84*length(profile_col))])
    # Get maximum offset size
    max_offset <- abs(diffs[2])
    # Calculate mean value of the cropped image given an offset of x
    offset_vals <- sapply(1:max_offset, function(x) mean(profile_col[x:(x + final_dims[2] - 1)]))
    # Find which offset is closest to the middle 2/3 of the image
    best_offset <- which.min((middle - offset_vals)^2)
    # Crop using the best offset
    new_img <- new_img[, best_offset:(best_offset + final_dims[2] - 1)]
  } else if (diffs[2] < 0) {
    # Need to pad
    # Get pad value - median pixel from 10 rows on each edge of the image
    value <- median(new_img[, c(1:10,(img_dims[2] - 10):img_dims[2])])
    # Create new image
    temp_img <- matrix(value, nrow = nrow(new_img), ncol = final_dims[2]) %>% as.Image()
    # Determine offset:
    offset <- floor(abs(diffs[2]/2))
    # Add in old values
    temp_img[, offset:(img_dims[2] + offset - 1)] <- new_img
    # Save bigger image as the new_img object
    new_img <- temp_img
  }

  new_img
}

