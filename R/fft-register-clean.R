#' FFT-based image registration
#'
#' @param im1 Image 1
#' @param im2 Image 2
#' @param signal What intensity is the signal in the image?
#'        1 = white (or light pixels), 0 = black (or dark pixels)
#' @param angle Recover angle alignment first?
#' @param theta_res Resolution of theta alignment in degrees
#' @return a list containing the aligned images, the angle value (if applicable), the shift value (if applicable), and the angle and shift FFT matrices
#' @export
#' @examples
#' img1 <- matrix(1, nrow = 1500, ncol = 1500) %>% EBImage::as.Image()
#' img1[400:1200, 700:900] <- 0
#' img1[300:350, 400:1200] <- 0
#' img1[830:850, 900:1050] <- .5
#' img2 <- img1 %>%
#'           img_translate(v = c(50, 50), output.dim = c(1500, 1500), bg.col = 1) %>%
#'           img_rotate(15, output.origin = c(750,750), output.dim = c(1500, 1500), bg.col = 1)
#'
#' noise1 <- sample(c(0, 1), length(img1), prob = c(.9, .1), replace = T) * rnorm(length(img1), sd = .15)
#' noise2 <- sample(c(0, 1), length(img1), prob = c(.9, .1), replace = T) * rnorm(length(img1), sd = .15)
#' img1 <- pmin(pmax(img1 + noise1, 0), 1)
#' img2 <- pmin(pmax(img2 + noise2, 0), 1)
#'
#' res <- fft_align(img1, img2)
#'
#' par(mfrow = c(2, 2))
#' # Initial images
#' plot(EBImage::rgbImage(img2, img1, pmin(img1, img2)))
#' # First, show what the post-angle-recovery images look like
#' plot(EBImage::rgbImage(res$angle[[1]], res$angle[[2]], pmin(res$angle[[1]], res$angle[[2]])))
#' # Fully transformed imgs
#' plot(EBImage::rgbImage(res$res[[1]], res$res[[2]], pmin(res$res[[1]], res$res[[2]])))
fft_align <- function(im1, im2, signal = 1, angle = T, theta_res = 1/3) {
  # Recover angle
  if (angle) {
    # Preprocess
    imlst_han <- fft_preprocess(list(im1, im2), signal = signal)
    # FFT
    imlst_fft_shift <- fft_and_shift(imlst_han)
    # Calculate magnitude and transform to polar
    imlst_polar <- mag_to_polar(imlst_fft_shift, ntheta = round(360/theta_res))

    # From here, could (in theory) call fft_align recursively
    res_angle <- fft_align(imlst_polar[[1]], imlst_polar[[2]], angle = F, signal = signal)
    # Hanning and such
    imlst_polar_pre <- fft_preprocess(imlst_polar, signal = 1)
    # FFT and shift
    imlst_polar_fft <- fft_and_shift(imlst_polar_pre)
    # Get cross-power, then invert
    angle_cross <- fft_crosspower(imlst_polar_fft)
    angle_ifft <- ifft(angle_cross, shift = F)
    # Get angle back out
    idxs <- get_value_idxs(angle_ifft)
    thetas <- seq(0, 2*pi, length.out = round(360/theta_res))*180/pi
    angle_val <- thetas[idxs[1]]
    # Rotate im2
    im2_adj <- img_rotate(im2, angle_val, bg.col = 0, output.dim = dim(im1))
  } else {
    im2_adj <- im2
  }

  # Hanning filter, adjust images to same size
  imlst_han <- fft_preprocess(list(im1, im2_adj), signal = signal)

  # FFT
  imlst_fft_shift <- fft_and_shift(imlst_han)

  # get cross-power and invert
  imlst_cross <- fft_crosspower(imlst_fft_shift)
  # Handle 0 values in cross-power spectrum. Not sure why this is necessary...
  imlst_cross[is.na(imlst_cross)] <- imlst_cross[which.min(Mod(imlst_cross))]
  imlst_ifft <- ifft(imlst_cross, shift = F)

  idxs <- get_value_idxs(imlst_ifft)
  translate_im2_by <- fix_coord(val = idxs, dim = dim(imlst_cross))

  res <- list(aligned = list(im1, img_translate(im2_adj, translate_im2_by, bg.col = 0, output.dim = dim(im1))))

  if (angle) {
    res$angle <- angle_val
    res$fft_angle <- angle_ifft
  }

  res$translate <- translate_im2_by
  res$fft_shift <- imlst_ifft

  res
}

fft_preprocess <- function(imlst, signal = 1) {
  # Force to use black-background encoding
  if (signal == 0) {
    imlst <- lapply(imlst, function(x) 1 - x)
  }

  # Pad so that images are the same size
  imlst_pad <- pad_img_match(imlst[[1]], imlst[[2]], value = 0)

  # Use hann filter
  imlst_han <- purrr::map(imlst_pad, im_hann_filter)

  imlst_han
}

fft_and_shift <- function(imlst) {
  purrr::map(imlst, stats::fft) %>%
    purrr::map(fftshift, inv = F)
}

fft_crosspower <- function(imlst) {
  x <- imlst[[1]] * Conj(imlst[[2]])
  y <- x/Mod(x)
  # y[is.na(y)] <- 0
  y
}

ifft <- function(x, shift = T) {
  # Use normalized inverse fft instead of R default unnormalized
  y <- Re(fft(x, inverse = T)/length(x))

  # Shift back
  if (shift) fftshift(y, inv = T) else y
}

mag_to_polar <- function(imlst, ntheta = 720) {
  # Get FT Mag:
  imlst_fft_mag <- purrr::map(imlst, Mod)
  imlst_polar <- purrr::map(imlst_fft_mag, img_to_polar, ntheta = ntheta)
}


im_hann_filter <- function(im, exaggerate = 1) {
  asp <- dim(im)[1]/dim(im)[2]
  widen_root <- c(pmax(asp, exaggerate), pmax(1/asp, exaggerate))

  row_frac <- seq(0, 2*pi, length.out = dim(im)[1])
  col_frac <- seq(0, 2*pi, length.out = dim(im)[2])

  hanning_row <- (.5*(1 - cos(row_frac)))^(1/widen_root[1])
  hanning_col <- (.5*(1 - cos(col_frac)))^(1/widen_root[2])

  mat <- hanning_row %*% t(hanning_col)
  mat_img <- im * 0
  mat_img[] <- mat

  mat_img * im
}

#' Pad each image so that the two are the same size
#'
#' @param im1 image 1
#' @param im2 image 2
#' @param value pixel value to use to pad (generally 0 if 1 is signal)
#' @export
pad_img_match <- function(im1, im2, value = 0) {
  img_size <- pmax(dim(im1), dim(im2))
  im1 <- img_pad_to_size(im1, img_size, value = value)
  im2 <- img_pad_to_size(im2, img_size, value = value)

  stopifnot(all.equal(dim(im1), dim(im2)))

  return(list(im1, im2))
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


#' FFT shift
#'
#' Shift zero-frequency component to center of the spectrum
#' @param mat Matrix from fft
#' @param inv Is this un-doing a previous shift? (if so, use ceiling instead of floor so everything is put back correctly)
#' @export
fftshift <- function(mat, inv = F) {
  if (inv) {
    f1 <- ceiling
  } else {
    f1 <- floor
  }

  newmat <- mat * 0
  imdim <- dim(mat)
  mid <- f1(imdim/2)

  # Define indexes for each quadrant of the image
  idx1 <- 1:mid[1]
  ridx1 <- sort(imdim[1] - idx1 + 1)
  idx2 <- 1:mid[2]
  ridx2 <- sort(imdim[2] - idx2 + 1)

  # Swap quadrants
  newmat[ridx1, ridx2] <- mat[idx1, idx2]
  newmat[idx1, idx2] <- mat[ridx1, ridx2]
  newmat[ridx1, idx2] <- mat[idx1, ridx2]
  newmat[idx1, ridx2] <- mat[ridx1, idx2]

  newmat
}

# find row and column of maximum value in the matrix
get_value_idxs <- function(mat, val = max(mat, na.rm = T)) {
  sapply(1:2, function(x) {
    apply(mat == val, x, function(xx) sum(xx)) %>%
      which.max()
  })
}

# Helper function
fix_coord <- function(val, dim) {
  ifelse(val > dim/2, val - dim + 1, val)
}


#' Convenience for plotting
#'
#' @param imlst List of images
#'
#' @description Plot a list of 2-3 images using rgbImage. If 2 images, the B channel will
#'   show the intersection between the two images. If there are 3 images, each
#'   image will be shown in a separate channel. Invisibly returns the RGB image.
#' @export
plot_imlist <- function(imlst) {
  stopifnot(length(imlst) >= 2 & length(imlst) < 4)
  if (length(imlst) == 2) imlst[[3]] <- pmin(imlst[[1]], imlst[[2]])
  plot(do.call(EBImage::rgbImage, imlst))
  invisible(do.call(EBImage::rgbImage, imlst))
}