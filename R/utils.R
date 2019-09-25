#' Image to data frame
#'
#' @param img image
#' @param filter_val value(s) of pixels to remove
#' @param row_neg negate the rows so that the image is "right side up"?
#' @export
image_to_df <- function(img, filter_val = 0, row_neg = F) {
  imdim <- dim(img)
  px_idx <- 0:(length(img) - 1)
  df <- cbind(
    row = px_idx %% imdim[1] + 1,
    col = c(1, -1)[row_neg + 1] * (floor(px_idx/imdim[1]) + 1),
    frame = floor(px_idx/(imdim[1]*imdim[2])) + 1,
    value = as.vector(img))

  if (!is.null(filter_val)) {
    df <- df[df[,4] != filter_val,]
  }

  res <- as.data.frame(df)
  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "convert to df",
                                        filter_val = 0)))
  res
}

#' Pad image so that specified coordinates are in the center of the image
#'
#' @param img image
#' @param center coordinates describing the desired "center" of the image (row, col)
#' @param value fill value to use for padding the image.
#'              If value has the same length as the number of frames in the
#'              image, it will be applied frame-wise.
#' @export
#' @examples
#' par(mfrow = c(1, 2))
#' im <- EBImage::readImage(system.file('images', 'nuclei.tif', package='EBImage'))[,,1:3]
#' EBImage::colorMode(im) <- "Color"
#' dim(im)
#' plot(im)
#'
#' im_pad <- pad_to_center(im, center = c(200, 200), value = c(0, 1, 0))
#' dim(im_pad)
#' plot(im_pad, all = T)
img_pad_to_center <- function(img, center = round(dim(img)/2), value = 0) {
  if (is.list(img)) {
    lapply(img, img_pad_to_size, size = size, value = value)
  }
  stopifnot(EBImage::is.Image(img))

  if (sum(center %% 1 > 0) > 0) {
    warning("non-integer center coordinates will be rounded to the nearest integer")
    center <- round(center)
  }

  # center <- abs(center)
  if (sum(center <= 0) != 0) stop("center coordinates must be nonzero.")

  dist_bottom_right <- dim(img)[1:2] - center
  dist_top_left <- center
  padding <- dist_bottom_right - dist_top_left
  pad_top <- pmax(0, padding[1])
  pad_bottom <- pmax(0, -padding[1])
  pad_left <- pmax(0, padding[2])
  pad_right <- pmax(0, -padding[2])

  img_pad(img, top = pad_top, bottom = pad_bottom, left = pad_left, right = pad_right, value = value)
}

#' Pad an image
#'
#' @param img image
#' @param top number of pixels to add to the top
#' @param bottom number of pixels to add to the bottom
#' @param left number of pixels to add to the left
#' @param right number of pixels to add to the right
#' @param value fill value to use for padding the image.
#'              If value has the same length as the number of frames in the
#'              image, it will be applied frame-wise.
#' @export
#' @importFrom EBImage getFrames colorMode Image
#' @importFrom abind abind
#' @examples
#' par(mfrow = c(1, 2))
#' im <- EBImage::readImage(system.file('images', 'nuclei.tif', package='EBImage'))[,,1:3]
#' EBImage::colorMode(im) <- "Color"
#' dim(im)
#' plot(im)
#'
#' im_pad <- img_pad(im, top = 15, bottom = 10, left = 5, right = 0, value = c(0, 1, 0))
#' dim(im_pad)
#' plot(im_pad)
img_pad <- function(img, top = 0, bottom = 0, left = 0, right = 0, value = 0) {
  if (is.list(img)) {
    return(lapply(img, img_pad, top = top, bottom = bottom, left = left,
                  right = right, value = value))
  }
  stopifnot(EBImage::is.Image(img))

  y <- EBImage::getFrames(img)

  y <- mapply(img_pad_frame, y, value,
              MoreArgs = list(top = top, bottom = bottom,
                              left = left, right = right), SIMPLIFY = F)

  res <- abind::abind(y, along = length(dim(img)), new.names = dimnames(img)) %>%
    EBImage::Image(colormode = EBImage::colorMode(img))
  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "pad",
                                        top_bottom = c(top, bottom),
                                        left_right = c(left, right),
                                        value = value)))

  res
}

img_pad_frame <- function(x, top = 0, bottom = 0, left = 0, right = 0, value = 0) {
  rbind(
    matrix(value, nrow = left, ncol = ncol(x) + top + bottom),
    cbind(matrix(value, ncol = top, nrow = nrow(x)),
          x,
          matrix(value, ncol = bottom, nrow = nrow(x))),
    matrix(value, nrow = right, ncol = ncol(x) + top + bottom)
  ) %>%
    EBImage::Image()
}

#' Resize an image (and record metadata)
#'
#' @param img Image
#' @param ... additional arguments to EBImage::resize: width w, height h,
#'        output.dim = c(w, h), output.origin = c(0, 0), antialias = F, ...
#'        (other arguments to affine, including bg.col, antialias, filter)
#' @return resized image with attribute "operation" that records the original
#'         and final dimensions of the image
#' @export
img_resize <- function(img, ...) {
  if (is.list(img)) {
    return(lapply(img, img_resize, ...))
  }
  stopifnot(EBImage::is.Image(img))

  res <- EBImage::resize(img, ...)
  args <- list(...)
  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "resize",
                                        orig_dim = dim(img),
                                        final_dim = dim(res),
                                        other_args = args)))

  res
}

#' Translate an image (and record metadata)
#'
#' @param img Image
#' @param v translation vector (2 numbers)
#' @param ... additional arguments to EBImage::translate: vector v, filter, ...
#'        (other arguments to affine, including bg.col, antialias, filter)
#' @return resized image with attribute "operation" that records the original
#'         and final dimensions of the image
#' @export
img_translate <- function(img, v, ...) {
  if (is.list(img)) {
    return(lapply(img, img_translate, ...))
  }
  stopifnot(EBImage::is.Image(img))

  res <- EBImage::translate(img, v = v, ...)
  args <- list(...)
  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "translate",
                                        vector = v, other_args = args)))

  res
}

#' Rotate an image (and record metadata)
#'
#' @param img Image
#' @param angle rotation angle in degrees
#' @param ... additional arguments to EBImage::rotate: output.dim,
#'        output.origin, ... (other arguments to affine, including bg.col,
#'        antialias, filter)
#' @return resized image with attribute "operation" that records the original
#'         and final dimensions of the image
#' @export
img_rotate <- function(img, angle, ...) {
  if (is.list(img)) {
    return(lapply(img, img_rotate, angle = angle, ...))
  }
  stopifnot(EBImage::is.Image(img))

  res <- EBImage::rotate(img, angle = angle, ...)
  args <- list(...)
  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "rotate",
                                        angle = angle,
                                        other_args = args)))

  res
}

#' Crop an image (and record metadata)
#'
#' @param img Image or list of images
#' @param dim New dimension of the image
#' @param center point to use around which cropping is symmetrical
#' @export
img_crop <- function(img, dim, center = NULL) {
  if (is.list(img)) {
    return(lapply(img, img_crop, dim = dim))
  }

  stopifnot(EBImage::is.Image(img))
  if (is.null(center)) center <- floor(dim(img)/2)

  current_dim <- dim(img)[1:2] # handle extra frames
  center <- center[1:2]

  to_crop <- current_dim - dim[1:2]
  if (all(to_crop == c(0, 0))) return(img)

  left_top_prop <- center/current_dim
  left_top <- floor(to_crop*left_top_prop)
  right_bottom <- ceiling(to_crop*(1 - left_top_prop))

  stopifnot(all.equal(to_crop, left_top + right_bottom)) # just to be sure

  left_top_coord <- pmin(pmax(left_top + 1, c(1, 1)), current_dim)
  right_bottom_coord <- pmin(pmax(current_dim - right_bottom, c(1, 1)), current_dim)

  img_frames <- getFrames(img)

  img_frame_crop <- mapply(crop_frame, img_frames,
                           MoreArgs = list(left_top = left_top_coord,
                                           right_bottom = right_bottom_coord),
                           SIMPLIFY = F)

  res <- abind::abind(img_frame_crop, along = length(dim(img)), new.names = dimnames(img)) %>%
    EBImage::Image(colormode = EBImage::colorMode(img))
  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "crop",
                                             old_dim = current_dim,
                                             new_dim = dim,
                                             center = center,
                                             top_corner = left_top_coord,
                                             bottom_corner = right_bottom_coord)))

  res
}

crop_frame <- function(x, left_top, right_bottom) {
  x[left_top[1]:right_bottom[1],left_top[2]:right_bottom[2]]
}



#' Image pyramid
#'
#' @param img image (or list of images). If image list is named, resulting
#'            tibble will have an extra column, img_name.
#' @param scale vector of numeric scaling factors
#' @return tibble containing columns img, scale, dim, and (if original image
#'         list is named) img_name. The img column will contain the scaled image
#' @export
img_pyramid <- function(img, scale, ...) {

  imgdf <- tibble::tibble(img = img)
  if (!is.null(names(img))) {
    imgdf$img_name <- names(img) %>% make.unique()
  }

  imgdf <- tidyr::crossing(imgdf, scale = scale) %>%
    dplyr::mutate(
      dim = purrr::map2(img, scale, ~floor(dim(.x)/.y)),
      img = purrr::map2(img, dim, ~img_resize(.x, w = .y[1], h = .y[2]))
    )

  stopifnot(c("img", "scale", "dim") %in% names(imgdf))

  return(imgdf)
}

#' Automatically resize image to a specific size - crop or pad as necessary
#'
#' This function resizes the image in rows first, and then columns. If the
#' final dimension is smaller than the current dimension, the offset which
#' minimizes the distance between the average pixel value of the middle 2/3 of
#' the original image and the average pixel value of the new image will be used.
#' If the final dimension is larger than the current dimension, the image will
#' be symmetrically padded on that dimension with pixels of the value specified
#' or (if not specified) the median pixel value of the 10 outmost rows/columns
#' of the image. This function operates with the assumption that the center of
#' the image contains the majority of the useful information.
#'
#' @param img image. If img has more than one frame (color or otherwise), the
#'          operation will be performed on each frame of the image separately.
#' @param final_dims numerical vector of length two giving the width and height
#'                   of the output image, respectively.
#' @param value fill value to use for padding the image if necessary (if NULL,
#'              will be automatically set to the median value of the 10 pixels
#'              on the left and right edge of the image). If value has the same
#'              length as the number of frames in img, value will be applied
#'              frame-wise.
#' @export
#' @importFrom EBImage getFrames colorMode Image
#' @importFrom abind abind
#' @examples
#'
#' par(mfrow = c(1, 3))
#' im <- EBImage::readImage(system.file('images', 'nuclei.tif', package='EBImage'))
#' dim(im)
#' plot(im, all = T)
#'
#' im_sm <- auto_resize_img(im, c(510, 490))
#' dim(im_sm)
#' plot(im_sm, all = T)
#'
#' im_big <- auto_resize_img(im, c(510, 550), value = c(0, .25, .5, .75))
#' dim(im_big)
#' plot(im_big, all = T)
auto_resize_img <- function(img, final_dims, value = NULL) {
  y <- EBImage::getFrames(img)

  if (is.null(value)) {
    y <- lapply(y, auto_resize_frame, final_dims = final_dims, value = value)
  } else {
    y <- mapply(auto_resize_frame, y, value,
                MoreArgs = list(final_dims = final_dims), SIMPLIFY = F)
  }

  res <- abind::abind(y, along = length(dim(img)), new.names = dimnames(img)) %>%
    EBImage::Image(colormode = EBImage::colorMode(img))

  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "resize",
                                        orig_dim = dim(img),
                                        final_dim = dim(res),
                                        value = value)))
  res
}

auto_resize_frame <- function(x, final_dims, value = NULL) {
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
    if (is.null(value)) {
      value <- median(new_img[c(1:10,(img_dims[1] - 10):img_dims[1]),])
    }
    # Create new image
    temp_img <- matrix(value, nrow = final_dims[1], ncol = ncol(new_img)) %>% EBImage::as.Image()
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
    # Get pad value - median pixel from 10 cols on each edge of the image
    if (is.null(value)) {
      value <- median(new_img[, c(1:10,(img_dims[2] - 10):img_dims[2])])
    }
    # Create new image
    temp_img <- matrix(value, nrow = nrow(new_img), ncol = final_dims[2]) %>% EBImage::as.Image()
    # Determine offset:
    offset <- floor(abs(diffs[2]/2))
    # Add in old values
    temp_img[, offset:(img_dims[2] + offset - 1)] <- new_img
    # Save bigger image as the new_img object
    new_img <- temp_img
  }

  new_img
}
