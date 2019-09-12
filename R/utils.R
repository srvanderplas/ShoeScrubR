#' Image to data frame
#'
#' @param img image
#' @param filter_val value(s) of pixels to remove
#' @export
image_to_df <- function(img, filter_val = 0) {
  tmp <- suppressMessages({
    img %>%
    EBImage::imageData() %>%
    tibble::as_tibble(., .name_repair = "unique") %>%
    dplyr::mutate(row = 1:dplyr::n()) %>%
    magrittr::set_names(gsub(x = names(.), pattern = "\\.{1,}", replacement = "")) %>%
    tidyr::gather(key = col, val = val, -row) %>%
    dplyr::mutate(col = as.numeric(col)) %>%
    dplyr::mutate(row = -row)
  })

  if (is.null(filter_val)) tmp

  tmp %>%
    dplyr::filter(!val %in% filter_val)
}


#' Pad image so that specified coordinates are in the center of the image
#'
#' @param img image
#' @param center coordinates describing the desired "center" of the image (row, col)
#' @param value value to use for padding
#' @export
pad_to_center <- function(img, center = round(dim(img)/2), value = 0) {
  if (sum(center %% 1 > 0) > 0) {
    warning("non-integer center coordinates will be rounded to the nearest integer")
    center <- round(center)
  }

  center <- abs(center)
  if (sum(center <= 0) != 0) stop("center coordinates must be nonzero.")

  dist_bottom_right <- dim(img) - center
  dist_top_left <- center
  padding <- dist_bottom_right - dist_top_left
  pad_top <- pmax(0, padding[1])
  pad_bottom <- pmax(0, -padding[1])
  pad_left <- pmax(0, padding[2])
  pad_right <- pmax(0, -padding[2])

  pad(img, top = pad_top, bottom = pad_bottom, left = pad_left, right = pad_right, value = value)
}

#' Pad an image
#'
#' @param img image
#' @param top number of pixels to add to the top
#' @param bottom number of pixels to add to the bottom
#' @param left number of pixels to add to the left
#' @param right number of pixels to add to the right
#' @param value value to pad with
#' @export
pad <- function(img, top = 0, bottom = 0, left = 0, right = 0, value = 0) {
  rbind(
    matrix(value, nrow = left, ncol = ncol(img) + top + bottom),
    cbind(matrix(value, ncol = top, nrow = nrow(img)),
          img,
          matrix(value, ncol = bottom, nrow = nrow(img))),
    matrix(value, nrow = right, ncol = ncol(img) + top + bottom)
  ) %>%
    Image()
}
