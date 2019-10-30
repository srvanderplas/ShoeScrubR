#' FFT-based image registration
#'
#' @param im1 Image 1
#' @param im2 Image 2
#' @param ... extra arguments passed to e.g. hanning_2d providing the widening
#'            exponents.
#' @export
#' @examples
#' img1 <- matrix(1, nrow = 150, ncol = 150) %>% EBImage::as.Image()
#' img1[40:120, 70:90] <- 0
#' img1[30:35, 40:120] <- 0
#' img1[83:85, 90:105] <- .5
#' img2 <- img_rotate(img1, 15, output.origin = c(75,75), output.dim = c(150, 150), bg.col = 1) %>%
#'           img_translate(v = c(5, 5), output.dim = c(150, 150), bg.col = 1)
#'
#' img1 <- img1 + rnorm(length(img1), sd = .15)
#' img2 <- img2 + rnorm(length(img1), sd = .15)
#'
#' plot(EBImage::rgbImage(img2, img1, pmin(img1, img2)))
#' res <- fft_align(img1, img2)
#' plot(EBImage::rgbImage(res[[1]], res[[2]], pmin(res[[1]], res[[2]])))
fft_align <- function(im1, im2, ...) {
  img_size <- pmax(dim(im1), dim(im2))
  im1 <- img_pad_to_size(im1, img_size, value = 1)
  im2 <- img_pad_to_size(im2, img_size, value = 1)

  # Apply Hanning window to reduce edge discontinuities
  img1 <- (1 - im1) * (hanning_2d(dim(im1), ...))
  img2 <- (1 - im2) * (hanning_2d(dim(im2), ...))

  # FFT image
  i1 <- fft(img1) %>% fftshift()
  i2 <- fft(img2) %>% fftshift()

  # Use the FT magnitude as a new image to get rotation
  m1 <- Re(sqrt(i1 * Conj(i1)))
  m2 <- Re(sqrt(i2 * Conj(i2)))

  nt <- 720
  i1_rot_polar <- img_to_polar(m1, ntheta = nt)
  i2_rot_polar <- img_to_polar(m2, ntheta = nt)

  thetas <- seq(0, 2*pi, length.out = nt)

  i1_rot_polar_hanning <- normalize(i1_rot_polar^.2) * (hanning_2d(dim(i1_rot_polar)))
  i2_rot_polar_hanning <- normalize(i2_rot_polar^.2) * (hanning_2d(dim(i2_rot_polar)))

  angle_reg <- fft_register(i1_rot_polar_hanning, i2_rot_polar_hanning)

  theta <- thetas[which(rowSums(angle_reg == max(angle_reg)) == 1)]*180/pi

  # affine_mat <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta), 0, 0),
                       # nrow = 3, ncol = 2, byrow = T)

  # if (theta %% 180 != 0) {
    im1_fix <- img_rotate(im1, angle = -theta,
                          bg.col = 1, output.dim = dim(im1),
                          output.origin = floor(dim(im1)/2))
  # } else {
  #   im1_fix <- im1
  # }

  img1_fix <- (1 - im1_fix) * (hanning_2d(dim(im1), ...))
  i1_fix <- fft(img1_fix)

  cross_power_rot <- i1_fix*Conj(i2)
  cross_power_mag <- Re(sqrt(cross_power_rot * Conj(cross_power_rot)))
  normalized_cross_power <- cross_power_rot/cross_power_mag


  ifft <- Re(fft(normalized_cross_power, inv = T))

  fix_coord <- function(val, dim) {
    ifelse(val > dim/2, val - dim, val)
  }

  solutions <- tibble::tibble(ifft_val = sort(ifft, decreasing = T)[1:10]) %>%
    dplyr::mutate(shift_row = purrr::map_dbl(ifft_val, ~which.max(rowSums(ifft == .))),
                  shift_col = purrr::map_dbl(ifft_val, ~which.max(colSums(ifft == .)))) %>%
    dplyr::mutate(shift_row = fix_coord(shift_row, dim(ifft)[1]),
                  shift_col = fix_coord(shift_col, dim(ifft)[2]))

  # affine_mat[3,] <- as.numeric(solutions[1,2:3])
  return(list(im1 = img_translate(im1_fix, v = -as.numeric(solutions[1, 2:3]),
                                  bg.col = 1, output.dim = dim(im1)),
              im2 = im2,
              theta = theta,
              solutions = solutions))
}

fft_register <- function(im1, im2) {
  # assumes hanning window already applied
  # FFT image
  i1 <- fft(im1)
  i2 <- fft(im2)

  cross_power <- i1*Conj(i2)
  cross_power_mag <- Re(sqrt(cross_power * Conj(cross_power)))
  normalized_cross_power <- cross_power/cross_power_mag


  ifft <- Re(fft(normalized_cross_power, inv = T))

  return(ifft)
}

#' Convert image to polar coordinates
#'
#' @param img Image
#' @param ntheta number of theta values to use
img_to_polar <- function(img, ntheta = 360) {

  center <- round(dim(img)/2)
  out_grid <- expand.grid(theta = seq(0, 2*pi, length.out = ntheta),
                          r = 1:ceiling(sqrt(sum(pmax(center, rev(dim(img)))^2))))
  cart_grid <- data.frame(x = with(out_grid, center[1] + r*cos(theta)),
                          y = with(out_grid, center[2] + r*sin(theta))) %>%
    dplyr::mutate(x = pmin(pmax(x, 1), dim(img)[1]), y = pmin(pmax(y, 1), dim(img)[2]))

  r <- matrix(0, nrow = length(unique(out_grid$theta)), ncol = length(unique(out_grid$r))) %>% EBImage::as.Image()
  r[] <- img[as.matrix(cart_grid)]

  return(EBImage::as.Image(r))
}

#' 2d Hanning Window
#'
#' Gives the coefficients for a hanning window for an image of size imdim.
#' @param imdim vector of length 2 giving the dimensions of the image
#' @param widen_root vector of length 2 giving the root to widen the image by
#'                   (in each dimension). It may be useful to have this be in
#'                   proportion to the image aspect ratio to ensure that the
#'                   filter is proportionate to the image.
#' @param ... extra parameters which will be dropped without warning
#' @export
hanning_2d <- function(imdim, widen_root = c(pmax(imdim[1]/imdim[2], 1), pmax(imdim[2]/imdim[1], 1)), ...) {
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
