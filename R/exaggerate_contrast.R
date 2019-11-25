#' Exaggerate darkness of light images
#'
#' This function finds the pixel value above which background is more likely,
#' and below which the print is more likely. The image values are shifted so
#' that this value is at 0.5, with maximum value 1 and minimum value 0, and then
#' the image is renormalized to fully utilize the 0-1 range.
#'
#' @param img Image
#' @param mask Mask, where object of interest has value 1 and everything else has value 0
#' @param eps threshold below which values are considered to be zero
#'
#' @export
#' @importFrom dplyr mutate filter lead lag
#' @importFrom EBImage normalize
#' @importFrom stats density
#'
#'
exaggerate_contrast <- function(img, mask, eps = 1e-4) {
  print_dens <- bkgd_dens <- x <- NULL
  likely_print <- stats::density(img[mask == 1], bw = .001, n = 1001, from =  0, to = 1)
  likely_bkgd <- stats::density(img[mask == 0], bw = .001, n = 1001, from =  0, to = 1)

  kdens_diff <- data.frame(x = likely_print$x,
                           print_dens = likely_print$y,
                           bkgd_dens = likely_bkgd$y) %>%
    dplyr::mutate(diff = print_dens - bkgd_dens) %>%
    dplyr::mutate(diff = ifelse(abs(diff) <= eps, 0, diff))

  icept1 <- kdens_diff %>%
    dplyr::filter(round(print_dens) > 0, round(bkgd_dens) > 0) %>%
    dplyr::mutate(lead = (dplyr::lead(x, 1) - x), lag = x - dplyr::lag(x, 1)) %>%
    dplyr::filter(lead == lag) %>%
    dplyr::filter(c(0, diff(sign(diff))) == -2)

  EBImage::normalize(pmax(pmin(img - icept1$x + .5, 1), 0))
}


