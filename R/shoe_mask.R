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
#'               column value, value will be used to weight the results.
#' @param weighted should weighted calculation be used?
#' @param ... additional arguments to image_to_df
#' @return angle to rotate the image by (in degrees)
#' @export
align_prcomp <- function(img = NULL, weighted = T, ...) {
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
#' @param img Image - greyscale, white as background, black as signal
#' @param mask Mask; white (signal) and black (bkgd)
#' @param thresh_pars list of parameters to use to clean the original image
#' @param exaggerate_pars list of parameters to use to exaggerate the original
#'        image to a mask-like appearance
#' @return list containing image and mask which have been aligned to center
#'         value and rotated via principal components so that PC1 is the
#'         positive y-axis.
#' @export
rough_align <- function(img, mask, thresh_pars = list(), exaggerate_pars = list()) {

  if (!all.equal(dim(img), dim(mask))) {
    message("Auto-resizing mask to the size of image, preserving mask scaling")
    mask <- auto_resize_img(mask, final_dims = dim(img))
  }


  exaggerate_pars$img <- img
  exag_img <- do.call("exaggerate_img_to_mask", exaggerate_pars)


  im_mode <- img_mode(img) # Get image fill color

  img_angle <- align_prcomp(exag_img)
  mask_angle <- align_prcomp(mask)

  img_align <- img_rotate(img, img_angle, bg.col = im_mode, output.dim = dim(img))
  exag_align <- img_rotate(exag_img, img_angle, bg.col = 0, output.dim = dim(img))
  mask_align <- img_rotate(mask, mask_angle, bg.col = 0, output.dim = dim(img))

  img_center <- binary_center(exag_align)
  mask_center <- binary_center(mask_align)

  trans_dist <- mask_center - img_center
  centered_mask <- img_translate(mask_align, -trans_dist, bg.col = 0, output.dim = dim(img))
  mask_center <- binary_center(centered_mask)

  padded_img <- img_pad_to_center(img_align, img_center, value = im_mode)
  padded_exag_img <- img_pad_to_center(exag_align, img_center)
  padded_mask <- img_pad_to_center(centered_mask, img_center)

  thresh_pars$img <- padded_img
  thresh_img <- do.call("clean_initial_img", thresh_pars)

  tibble(name = c("img", "thresh", "mask"),
         img = list(padded_img, thresh_img, padded_mask))
}

#' Invert, smooth, and threshold an image
#'
#' @param img Image
#' @param gaussian_d diameter of brush to use for gaussian blur
#' @param threshold_val threshold value to use on normalized, inverted, blurred
#'        image
#' @export
clean_initial_img <- function(img, gaussian_d = 25, threshold_val = .15) {
  img %>%
    EBImage::filter2(EBImage::makeBrush(gaussian_d, "gaussian")) %>%
    EBImage::normalize() %>%
    magrittr::subtract(1, .) %>%
    (function(.) . > threshold_val)
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
  clean_initial_img(img, gaussian_d = gaussian_d, threshold_val = threshold_val) %>%
    EBImage::opening(EBImage::makeBrush(opening_d, "disc")) %>%
    EBImage::closing(EBImage::makeBrush(closing_d, "disc"))
}

#' Binary image center
#'
#' @param x image/matrix
#' @export
binary_center <- function(x) {
  d1sum <- apply(x != 0, 1, sum)
  d1sum <- d1sum/sum(d1sum)
  d1sum <- sum((1:(dim(x)[1])) * d1sum)
  d2sum <- apply(x != 0, 2, sum)
  d2sum <- d2sum/sum(d2sum)
  d2sum <- sum((1:(dim(x)[2])) * d2sum)

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
