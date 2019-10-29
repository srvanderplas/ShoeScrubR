#' FFT-based image registration
#'
#' @param im1 Image 1
#' @param im2 Image 2
#' @export
fft_register <- function(im1, im2) {
  img_size <- pmax(dim(im1), dim(im2))
  im1 <- img_pad_to_size(im1, img_size, value = 1)
  im2 <- img_pad_to_size(im2, img_size, value = 1)

  # Apply Hanning window to reduce edge discontinuities
  img1 <- im1 * (1 - hanning_2d(dim(im1)))
  img2 <- im2 * (1 - hanning_2d(dim(im2)))

  # FFT image
  i1 <- fft(img1)
  i2 <- fft(img2)
  cross_power <- i1*Conj(i2)
  normalized_cross_power <- cross_power/sqrt(cross_power * Conj(cross_power))



  ifft <- Re(fft(ecp/sqrt(ecp * Conj(ecp)), inv = T))



}

#' Convert image to polar coordinates
#'
#' @param img Image
img_to_polar <- function(img) {

}

#' 2d Hanning Window
#'
#' Gives the coefficients for a hanning window for an image of size imdim.
#' @param imdim vector of length 2 giving the dimensions of the image
#' @param widen_root vector of length 2 giving the root to widen the image by
#'                   (in each dimension). It may be useful to have this be in
#'                   proportion to the image aspect ratio to ensure that the
#'                   filter is proportionate to the image.
#' @export
hanning_2d <- function(imdim, widen_root = c(2, 1)) {
  hanning_row <- 1 - cos(2*pi*(1:imdim[1])/imdim[1])
  hanning_col <- 1 - cos(2*pi*(1:imdim[2])/imdim[2])

  mat <- (.5)^(1/mean(widen_root)) *
    (hanning_row)^(1/widen_root[1]) %*%
    (t(hanning_col))^(1/widen_root[2])

  EBImage::as.Image(mat)
}

#' FFT shift
#'
#' Shift zero-frequency component to center of the spectrum
#' @param mat Matrix from fft
#' @export
fftshift <- function(mat) {
  newmat <- mat * 0
  mid <- round(dim(mat)/2)
  imdim <- dim(mat)
  idx1 <- 1:mid[1]
  ridx1 <- rev(imdim[1] - idx1 + 1)
  idx2 <- 1:mid[2]
  ridx2 <- rev(imdim[2] - idx2 + 1)
  newmat[ridx1, ridx2] <- mat[idx1, idx2]
  newmat[idx1, idx2] <- mat[ridx1, ridx2]
  newmat[ridx1, idx2] <- mat[idx1, ridx2]
  newmat[idx1, ridx2] <- mat[ridx1, idx2]

  newmat
}