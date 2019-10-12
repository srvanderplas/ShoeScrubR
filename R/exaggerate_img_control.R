#' Exaggerate an image to a mask-like appearance with control of parameters
#'
#' @param img Image
#' @param gaussian_d diameter of brush to use for gaussian blur
#' @param threshold_val threshold value to use on normalized, inverted, blurred
#'        image
#' @param opening_d diameter to use for image opening (despeckling)
#' @param closing_d diameter to use for image closing (exaggeration)
#' @export
exaggerate_img_control <- function(img, gaussian_d = 125, threshold_val = .125,
                                   opening_d = 7, closing_d = 301) {

  tmp <- clean_initial_img(img, gaussian_d = gaussian_d,
                           threshold_val = threshold_val)

  tmp <- img_open_close(tmp, opening_d = opening_d, closing_d = closing_d)

  tmp
}

img_open_close <- function(img, opening_d, closing_d, opening_shape = "disc", closing_shape = "disc") {
  if (is.list(img)) {
    return(lapply(img, img_open_close, opening_d = opening_d, closing_d = closing_d,
                  opening_shape = opening_shape, closing_shape = closing_shape))
  }

  tmp <- img %>%
    EBImage::opening(EBImage::makeBrush(opening_d, shape = opening_shape)) %>%
    EBImage::closing(EBImage::makeBrush(closing_d, shape = closing_shape))


  attr(tmp, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "exaggerate",
                                             opening_d = opening_d,
                                             closing_d = closing_d)))
  tmp
}

#' Invert, smooth, and threshold an image
#'
#' @param img Image
#' @param gaussian_d diameter of brush to use for gaussian blur
#' @param threshold_val threshold value to use on normalized, inverted, blurred
#'        image. If threshold is 0, it will be automatically determined using
#'        a heuristic approach.
#' @export
clean_initial_img <- function(img, gaussian_d = 25, threshold_val = .15) {
  . <- NULL
  if (is.list(img)) {
    return(lapply(img, clean_initial_img, gaussian_d = gaussian_d, threshold_val = threshold_val))
  }

  tmp <- img %>%
    EBImage::filter2(EBImage::makeBrush(gaussian_d, "gaussian")) %>%
    EBImage::normalize() %>%
    magrittr::subtract(1, .)

  if (threshold_val == 0) {
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

#' Get the center of a binary image
#'
#' @description Calculate the center of mass of a binary image, with or without trimming
#' @param img image/matrix
#' @param trim Trim 5\% from each side of the image? (Useful for removing page boundary issues)
#'
#' @export
binary_center <- function(img, trim = T) {
  if (is.list(img)) {
    return(lapply(img, binary_center, trim = trim))
  }

  # stopifnot(all(unique(img) %in% c(0, 1)))

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
