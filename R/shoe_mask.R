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
  img_resize(im, w = output_dim[1], h = output_dim[2])
}

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
  if (is.list(img)) {
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
#' @param exaggerate_pars list of parameters to use to exaggerate the original
#'        image to a mask-like appearance
#' @param img_fill_value value to fill any image padding with; defaults to the
#'        mode value (as calculated by \code{\link{img_mode}}).
#' @return list containing image and mask which have been aligned to center
#'         value and rotated via principal components so that PC1 is the
#'         positive y-axis.
#' @export
rough_align <- function(img, mask, exaggerate_pars = list(),
                        img_fill_value = img_mode(img)) {

  if (!all(dim(img) == dim(mask))) {
    message("Auto-resizing mask to the size of image, preserving mask scaling")
    mask <- auto_resize_img(mask, final_dims = dim(img))
  }

  exaggerate_pars$img <- img
  exag_img <- do.call("exaggerate_img_to_mask", exaggerate_pars)

  # plot(rgbImage(img, exag_img, mask))

  img_angle <- align_prcomp(exag_img)
  mask_angle <- align_prcomp(mask)

  img_align <- img_rotate(img, img_angle, bg.col = img_fill_value, output.dim = dim(img))
  exag_align <- img_rotate(exag_img, img_angle, bg.col = 0, output.dim = dim(img))
  mask_align <- img_rotate(mask, mask_angle, bg.col = 0, output.dim = dim(img))

  # plot(rgbImage(img_align, exag_align, mask_align))

  # plot(rgbImage(img, exag_img, mask))
  shifts <- calc_shifts(exag_align, mask_align)

  padded_img <- img_pad(img_align, padding = shifts$img, value = img_fill_value)
  padded_exag_img <- img_pad(exag_align, padding = shifts$img, value = 0)
  padded_mask <- img_pad(mask_align, padding = shifts$mask, value = 0)

  # plot(rgbImage(padded_img, padded_exag_img, padded_mask))

  list(img = padded_img, exag_img = padded_exag_img, mask = padded_mask)
}

calc_shifts <- function(img, mask, trim_mask = F) {

  img_center <- binary_center(img)
  mask_center <- binary_center(mask, trim = trim_mask)

  mask_pad_left_top <- pmax(img_center, mask_center) - mask_center
  mask_pad_right_bottom <- mask_center - pmin(img_center, mask_center)

  img_pad_left_top <- pmax(img_center, mask_center) - img_center
  img_pad_right_bottom <- img_center - pmin(img_center, mask_center)

  if (!all.equal(img_center + img_pad_left_top,
                 mask_center + mask_pad_left_top)) {
    warning("Image center alignment appears to have failed")
  }
  return(list(img = c(top = img_pad_left_top[2],
                      bottom = img_pad_right_bottom[2],
                      left = img_pad_left_top[1],
                      right = img_pad_right_bottom[1]),
         mask = c(top = mask_pad_left_top[2],
                  bottom = mask_pad_right_bottom[2],
                  left = mask_pad_left_top[1],
                  right = mask_pad_right_bottom[1])))
}


#' Invert, smooth, and threshold an image
#'
#' @param img Image
#' @param gaussian_d diameter of brush to use for gaussian blur
#' @param threshold_val threshold value to use on normalized, inverted, blurred
#'        image
#' @export
clean_initial_img <- function(img, gaussian_d = 25, threshold_val = .15) {
  if (is.list(img)) {
    return(lapply(img, clean_initial_img, gaussian_d = gaussian_d, threshold_val = threshold_val))
  }

  tmp <- img %>%
    EBImage::filter2(EBImage::makeBrush(gaussian_d, "gaussian")) %>%
    EBImage::normalize() %>%
    magrittr::subtract(1, .)

  if (is.null(threshold_val) | threshold_val == 0) {
    thresholds <- seq(.05, .5, .01)
    t_mean <- sapply(thresholds, function(x) mean(tmp > x))
    threshold_val <- thresholds[which.min(abs(t_mean - 0.055))]
    message("Image cleaned using a threshold of ", threshold_val)
  }
  tmp <- tmp %>%
    (function(.) . > threshold_val)

  attr(tmp, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "clean",
                                             gaussian_d = gaussian_d,
                                             threshold_val = threshold_val)))

  tmp
}

#' Exaggerate an image to a mask-like appearance
#'
#' @param img Image
#' @param gaussian_d diameter of brush to use for gaussian blur
#' @param threshold_val threshold value to use on normalized, inverted, blurred
#'        image
#' @param opening_d diameter to use for image opening (despeckling)
#' @param closing_d diameter to use for image closing (exaggeration)
#' @export
exaggerate_img_to_mask <- function(img, gaussian_d = 125, threshold_val = .125,
                                   opening_d = 7, closing_d = 301) {
  if (is.list(img)) {
    return(lapply(img, exaggerate_img_to_mask, gaussian_d = gaussian_d,
                  threshold_val = threshold_val, opening_d = opening_d,
                  closing_d = closing_d))
  }

  tmp <- clean_initial_img(img, gaussian_d = gaussian_d, threshold_val = threshold_val) %>%
    EBImage::opening(EBImage::makeBrush(opening_d, "disc")) %>%
    EBImage::closing(EBImage::makeBrush(closing_d, "disc"))

  attr(tmp, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "exaggerate",
                                             opening_d = opening_d,
                                             closing_d = closing_d)))

  tmp
}

#' Binary image center
#'
#' @param img image/matrix
#' @param trim Trim 5% from each side of the image? (Useful for removing page boundary issues)
#' @export
binary_center <- function(img, trim = T) {
  if (is.list(img)) {
    return(lapply(img, binary_center, trim = trim))
  }

  d1sum <- apply(img != 0, 1, sum)
  if (trim) {
    d1sum_trim <- rep(0, length(d1sum))
    d1sum_trim[ceiling(.05*length(d1sum)):floor(.95*length(d1sum))] <- 1
  } else {
    d1sum_trim <- 1
  }

  d1sum <- d1sum*d1sum_trim
  d1sum <- d1sum/sum(d1sum)
  d1sum <- sum((1:(dim(img)[1])) * d1sum)

  d2sum <- apply(img != 0, 2, sum)
  if (trim) {
    d2sum_trim <- rep(0, length(d2sum))
    d2sum_trim[ceiling(.05*length(d2sum)):floor(.95*length(d2sum))] <- 1
  } else {
    d2sum_trim <- 1
  }

  d2sum <- d2sum*d2sum_trim
  d2sum <- d2sum/sum(d2sum)
  d2sum <- sum((1:(dim(img)[2])) * d2sum)

  round(c(d1sum, d2sum))
}

#' Get most common pixel value (approx)
#'
#' @param img Image
#' @param digits precision of numerical values
#' @param size size of sample - if img has more than size pixels, it will be
#'        reduced to a random sample of size before tabulation
#' @export
img_mode <- function(img, digits = 2, size = 50000) {
  if (is.list(img)) {
    return(lapply(img, img_mode, digits = digits, size = size))
  }

  v <- img %>%
    EBImage::normalize() %>%
    EBImage::imageData() %>%
    as.vector()

  if (length(v) > size) v <- sample(v, size = size, replace = F)
  v %>%
    round(., digits = digits) %>%
    table() %>%
    sort(decreasing = T) %>%
    names() %>% `[`(1) %>%
    as.numeric()
}
