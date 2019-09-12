#' Create initial mask for film and powder print
#'
#' Creates initial mask by blurring image, then inverting, then thresholding.
#' @param img image
#' @param blur width of blur brush
#' @param w width of threshold region
#' @param h height of threshold region
#' @param offset threshold offset
#' @export
#' @return mask of film and powder print image
#' TODO: Add example
film_mask <- function(img, blur = 45, w = 100, h = 100, offset = 0.0125) {
  blur_brush <- EBImage::makeBrush(blur, shape = "disc", step = T)^2
  blur_brush <- blur_brush/sum(blur_brush)

  img %>%
    EBImage::filter2(filter = blur_brush) %>%
    (function(x) (1 - x)) %>%
    EBImage::thresh(w = w, h = h, offset = offset)
}

#' Prunes small, disjoint regions from the mask
#'
#' @param mask image mask with disjoint regions labeled using EBImage::bwlabel()
#'        or other similar function
#' @param prop_limit areas which are smaller than this proportion of the image
#'        will be pruned if they are disjoint
#'
#' @export
film_prune <- function(mask, prop_limit = 0.02) {
  counts <- table(mask)
  categories <- sort(unique(as.numeric(mask)))
  prop <- counts/sum(counts)
  cat_replace <- categories[prop < prop_limit]
  mask[mask %in% cat_replace] <- 0
  mask
} # Get rid of really small areas - assumes large areas have merged...

#' Cleans up initial mask using dilation, erosion, and pruning of small regions
#'
#' @param mask binary image (usually from film_mask)
#' @param d1 diameter to use for mask erosion
#' @param d2 diameter to use for mask dilation
#' @param prop maximum proportion of the image which can be pruned
#' @export
film_mask_clean <- function(mask, d1 = 5, d2 = 91, prop = 1.5*pi*d2^2/length(mask)) {
  f1 <- makeBrush(d1, shape = "disc")
  f2 <- makeBrush(d2, shape = "disc")
  mask %>%
    EBImage::erode(f1) %>%
    EBImage::dilate(f2) %>%
    EBImage::bwlabel() %>%
    film_prune(prop_limit = prop)
}

#' Expands image mask
#'
#' This function "fills in" any pixels which are not part of the mask but are
#' between selected pixels in the same row or column. The resulting image is
#' then dilated using a disc-shaped brush of expand pixels
#' @param mask image mask
#' @param expand diameter of expansion brush. Should be odd if nonzero
#' @examples
#' mat <- matrix(0, nrow = 30, ncol = 30)
#' mat[13:17,13:17] <- matrix(c(0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0), ncol = 5, byrow = T)
#' mat <- mat %>% EBImage::as.Image()
#' par(mfrow = c(2, 2))
#' plot(mat)
#' plot(film_expand_mask(mat, expand = 0))
#' plot(film_expand_mask(mat, expand = 3))
#' @export
film_expand_mask <- function(mask, expand = 51) {

  useful <- mask != 0

  if (sum(useful) > 0) {
    useful_df <- tibble::as_tibble(useful, rownames = NA) %>%
      dplyr::mutate(row = 1:n()) %>%
      tidyr::gather(key = column, value = value, -row) %>%
      dplyr::mutate(column = column %>% str_remove("V") %>%
                      readr::parse_integer()) %>%
      dplyr::filter(value > 0)

    useful_cols <- useful_df %>%
      dplyr::group_by(row) %>%
      dplyr::summarize(mincol = min(column), maxcol = max(column)) %>%
      dplyr::mutate(column = purrr::map2(
        mincol, maxcol, ~tibble::tibble(column = .x:.y))) %>%
      dplyr::select(-mincol, -maxcol) %>%
      tidyr::unnest(column) %>%
      dplyr::select(row, column)

    useful_rows <- useful_df %>%
      dplyr::group_by(column) %>%
      dplyr::summarize(minrow = min(row), maxrow = max(row)) %>%
      dplyr::mutate(row = purrr::map2(
        minrow, maxrow, ~tibble::tibble(row = .x:.y))) %>%
      dplyr::select(-minrow, -maxrow) %>%
      tidyr::unnest(row) %>%
      dplyr::select(row, column)

    tmp <- bind_rows(useful_rows, useful_cols) %>%
      unique() %>%
      dplyr::mutate(value = 1) %>%
      tidyr::complete(row = 1:nrow(mask), column = 1:ncol(mask), fill = list(value = 0)) %>%
      tidyr::spread(key = column, value = value, fill = NA) %>%
      dplyr::arrange(row) %>%
      dplyr::select(-row) %>%
      as.matrix() %>%
      EBImage::as.Image()

  } else {
    tmp <- mask
  }

  if (expand > 0) {
    EBImage::dilate(tmp, kern = makeBrush(expand, "disc"))
  } else {
    tmp
  }

} # Shorthand for erode, dilate, bwlabel, then "convex hull" calculation

#' Mask borders of an image
#'
#' @param img image
#' @param d distance from border to mask
#' @param fill color to mask with
#' @export
mask_borders <- function(img, d = 30, fill = median(image)) {
  img[1:d,] <- img[(nrow(img) - d):nrow(img),] <- fill
  img[,1:d] <- img[,(ncol(img) - d):ncol(img)] <- fill

  img
}


film_expand_mask_alt <- function(mask, expand = 51) {

  test_mat <- 0*mask

  nr <- nrow(mask)
  nc <- ncol(mask)

  for (i in 1:nrow(mask)) {
    ne0 <- which(mask[i,] != 0)

    test_mat[i, min(ne0):max(ne0)] <- test_mat[i, min(ne0):max(ne0)] + 1
  }

  for (j in 1:ncol(mask)) {
    ne0 <- which(mask[,j] != 0)

    test_mat[min(ne0):max(ne0), j] <- test_mat[min(ne0):max(ne0), j]  + 1
  }

  test_mat > 0
}
