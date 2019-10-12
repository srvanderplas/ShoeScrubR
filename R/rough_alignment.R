#' Align an image and a mask based on principal components
#'
#' One of img or img_df must be supplied. If img is supplied, additional
#' arguments to img_to_df may also be supplied using ....
#'
#' @param img image or data frame of the locations of nonzero pixels in an image
#'               (columns row, col, value*). If img is a data frame and has
#'               column value, value will be used to weight the results. If img
#'               is a list, each image in the list will be handled separately
#' @param weighted should weighted calculation be used?
#' @param ... additional arguments to image_to_df
#' @return angle to rotate the image by (in degrees)
#' @export
align_prcomp <- function(img = NULL, weighted = T, ...) {
  if (is.list(img) & !is.data.frame(img)) {
    return(lapply(img, align_prcomp, weighted = weighted, ...))
  }

  if (EBImage::is.Image(img)) {
    img_df <- image_to_df(img, ...)
  } else {
    img_df <- img
  }
  stopifnot(!is.null(img_df))

  if (hasName(img_df, "value") & weighted) {
    rowmean <- mean(img_df$row, weight = img_df$value, na.rm = T)
    colmean <- mean(img_df$col, weight = img_df$value, na.rm = T)
    img_df$row <- img_df$row - rowmean
    img_df$col <- img_df$col - colmean
    weight <- img_df$value/sum(img_df$value)
    center_vals <- c(rowmean, colmean)
    center <- F
  } else {
    weight <- rep(1, nrow(img_df))
    center <- T
  }

  df <- img_df[,c("row", "col")] %>% as.matrix()

  pca <- prcomp(df*sqrt(weight), center = center, scale = F)

  if (sum(sign(pca$rotation)) < -1) pca$rotation <- pca$rotation * -1 # Make most PC vals positive

  # if (!center) pca$center <- center_vals

  # pca

  pca_to_angle(pca)
}


pca_to_angle <- function(rot) {
  if (hasName(rot, "rotation")) {
    rot <- rot$rotation
  }
  stopifnot(all.equal(dim(rot), c(2,2)))

  if (sum(sign(rot)) < -1) rot <- rot * -1 # Make most PC vals positive

  if (abs(rot[1,1]) > abs(rot[2,1])) mag <- asin(abs(rot[1,1])) else mag <- acos(abs(rot[2,1]))
  if (sign(rot[1,1]) == sign(rot[2,1])) coef <- 1 else coef <- -1

  return(mag*coef * 180/pi)
}

#' Initial alignment of mask and image
#'
#' @param img Image - greyscale, white(ish) as background, black as signal
#' @param mask Mask; white (signal) and black (bkgd)
#' @param img_fill_value value to fill any image padding with; defaults to the
#'        mode value (as calculated by \code{\link{img_mode}}).
#' @param ...
#' @return list containing image and mask which have been aligned to center
#'         value and rotated via principal components so that PC1 is the
#'         positive y-axis.
#' @export
rough_align <- function(img, mask, img_fill_value = img_mode(img), ...) {

  if (!all(dim(img) == dim(mask))) {
    message("Auto-resizing mask to the size of image, preserving mask scaling")
    mask <- auto_resize_img(mask, final_dims = dim(img))
  }

  exag_img <- exaggerate_img_auto(img, ...)

  # plot(rgbImage(img, exag_img, mask))

  img_angle <- align_prcomp(exag_img)
  mask_angle <- align_prcomp(mask)

  img_align <- img_rotate(img, img_angle, bg.col = img_fill_value, output.dim = dim(img))
  exag_align <- img_rotate(exag_img, img_angle, bg.col = 0, output.dim = dim(img))
  mask_align <- img_rotate(mask, mask_angle, bg.col = 0, output.dim = dim(img))

  exag_align_center <- get_mask_arch(exag_align)
  mask_align_center <- get_mask_arch(mask_align)

  # plot(rgbImage(img_align, exag_align, mask_align))

  # plot(rgbImage(img, exag_img, mask))
  shifts <- calc_shifts(exag_align, mask_align, img_center = exag_align_center, mask_center = mask_align_center)

  padded_img <- img_pad(img_align, padding = shifts$img, value = img_fill_value)
  padded_exag_img <- img_pad(exag_align, padding = shifts$img, value = 0)
  padded_mask <- img_pad(mask_align, padding = shifts$mask, value = 0)

  # plot(rgbImage(padded_img, padded_exag_img, padded_mask))

  list(img = padded_img, exag_img = padded_exag_img, mask = padded_mask)
}

#' Get arch location (approx) in a shoe mask
#'
#' This function assumes that the shoe is oriented such that the toe-heel axis
#' is vertical. The arch would then be the local minimum of mask width near the
#' center of the image.
#'
#' @param img Image
#' @return center (row, col) of the location of the arch in the mask
#' @examples
#' library(ShoeScrubR)
#' get_mask("Nike", 10, "L") %>%
#' get_mask_arch()
get_mask_arch <- function(img) {
  mask_width <- img %>% colSums() %>% smooth.spline(df = 10)

  d1 <- mask_width %>% predict(deriv = 1) %>% tibble::as_tibble() %>%
    dplyr::rename(idx = x, d1 = y) %>%
    dplyr::mutate(
      d1 = round(d1, 8),
      fd_sign = sign(d1) - sign(dplyr::lag(d1, 1)))
  d2 <- mask_width %>% predict(deriv = 2) %>% tibble::as_tibble() %>%
    dplyr::rename(idx = x, d2 = y) %>%
    dplyr::mutate(d2 = round(d2, 8))
  ds <- dplyr::left_join(d1, d2, by = "idx") %>%
    dplyr::filter(d2 > 0, fd_sign > 0) %>%
    dplyr::filter(dplyr::row_number() == which.min((idx - dim(img)[2]/2)^2))

  if (nrow(ds) > 0) {
    row <- ds$idx
    col <- round(mean(which(img[,row] > 0)))
  } else {
    row <- col <- NaN
  }

  if (nrow(ds) == 0 | is.nan(row) | is.nan(col)) {
    center <- binary_center(img, trim = T)
    row <- ifelse(is.nan(row), center[2], row)
    col <- ifelse(is.nan(col), center[1], col)
  }

  c(col, row)
}

calc_shifts <- function(img, mask,
                        img_center = binary_center(img),
                        mask_center = binary_center(mask, trim = trim_mask),
                        trim_mask = F) {

  mask_pad_left_top <- pmax(img_center, mask_center) - mask_center
  mask_pad_right_bottom <- mask_center - pmin(img_center, mask_center)

  img_pad_left_top <- pmax(img_center, mask_center) - img_center
  img_pad_right_bottom <- img_center - pmin(img_center, mask_center)

  # if (!all.equal(img_center + img_pad_left_top,
  #                mask_center + mask_pad_left_top)) {
  #   warning("Image center alignment appears to have failed")
  # }
  return(list(img = c(top = img_pad_left_top[2],
                      bottom = img_pad_right_bottom[2],
                      left = img_pad_left_top[1],
                      right = img_pad_right_bottom[1]),
              mask = c(top = mask_pad_left_top[2],
                       bottom = mask_pad_right_bottom[2],
                       left = mask_pad_left_top[1],
                       right = mask_pad_right_bottom[1])))
}
