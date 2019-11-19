#' FFT-based image registration
#'
#' @param im1 Image 1
#' @param im2 Image 2
#' @param signal What intensity is the signal in the image? 1 = white, 0 = black
#' @param show_bg Should rotations and translations of images use a visible background?
#' @param angle.lim Restrict search to angles between (-angle.lim, angle.lim)
#'                  (in radians). Defaults to pi/18, or ~10 degrees.
#' @param trans.lim Restrict search to coordinates between
#'                  (-trans.lim[1], trans.lim[1]) x (-trans.lim[2], trans.lim[2]).
#'                  Defaults to (200, 200)
#' @param ... extra arguments passed to e.g. hanning_2d providing the widening
#'            exponents.
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
#' res <- fft_align(1 - img1, 1 - img2, angle.lim = pi)
#'
#' par(mfrow = c(2, 2))
#' # Initial images
#' plot(EBImage::rgbImage(img2, img1, pmin(img1, img2)))
#' # First, show what the post-angle-recovery images look like
#' plot(EBImage::rgbImage(res$angle[[1]], res$angle[[2]], pmin(res$angle[[1]], res$angle[[2]])))
#' # Fully transformed imgs
#' plot(EBImage::rgbImage(res$res[[1]], res$res[[2]], pmin(res$res[[1]], res$res[[2]])))
fft_align <- function(im1, im2, signal = 1, show_bg = F,
                      angle.lim = pi/18, trans.lim = c(200, 200), ...) {

  # Ensure same size, origin
  ims <- pad_img_match(im1, im2)
  img_size <- dim(ims[[1]])
  img_origin <- round(img_size/2)

  # Keep bkgd?
  bg_col <- ifelse(show_bg, 0.5, 1 - signal)

  # Apply Hanning window to reduce edge discontinuities
  ims_han <- hanning_list(ims, ...)

  # FFT image
  ffts <- fft_list(ims_han)

  # Recover angle via polar transform of fft magnitude
  nt <- 720
  thetas <- seq(0, 2*pi, length.out = nt)
  angle_reg <- angle_recover(ffts, nt = nt)
  angle_reg2 <- angle_reg
  if (!is.null(angle.lim)) {
    # Enforce angle constraint
    invalid_thetas <- cos(thetas) < cos(angle.lim)
    angle_reg2[invalid_thetas,] <- -Inf
  }
  theta <- thetas[get_value_idxs(angle_reg2)[1]]*180/pi

  # Update with rotated image 2
  ims_fix <- ims
  ims_fix[[2]] <- img_rotate(ims[[2]], angle = theta, bg.col = bg_col,
                             output.dim = img_size,
                             output.origin = img_origin)

  ims_han_fix <- hanning_list(ims_fix, invert = signal != 1, ...)

  # Recover shift
  ifft <- do.call(fft_register, ims_han_fix) # 1 on 2
  if (!is.null(trans.lim)) {
    # Enforce translation constraint
    idx_mat <- ifft * 0 - Inf
    valid_shifts_x <- 1:trans.lim[1]
    valid_shifts_y <- 1:trans.lim[2]
    idx_mat[valid_shifts_x, valid_shifts_y] <- 0
    idx_mat[img_size[1] - valid_shifts_x + 1, valid_shifts_y] <- 0
    idx_mat[valid_shifts_x, valid_shifts_y] <- 0
    idx_mat[img_size[1] - valid_shifts_x + 1, img_size[2] - valid_shifts_y + 1] <- 0
  } else {
    idx_mat <- ifft * 0
  }

  transform_2_by <- mapply(fix_coord, get_value_idxs(ifft + idx_mat), img_size)

  ims_fix_trans <- ims_fix
  ims_fix_trans[[1]] <- img_translate(ims_fix[[1]], v = -pmin(0, transform_2_by),
                                      bg.col = bg_col, output.dim = img_size)
  ims_fix_trans[[2]] <- img_translate(ims_fix[[2]], v = pmax(0, transform_2_by),
                                      bg.col = bg_col, output.dim = img_size)

  return(list(res = ims_fix_trans,
              angle = ims_fix, # debugging purposes
              orig = ims, # debugging purposes
              angle_fft_res = angle_reg, # debugging purposes
              rot_img_2_by = theta,
              shift_fft_res = ifft, # debugging purposes
              translate_img_2_by = transform_2_by))
}

# find row and column of maximum value in the matrix
get_value_idxs <- function(mat, val = max(mat)) {
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
#' Plot a list of 2-3 images using rgbImage. If 2 images, the B channel will
#' show the intersection between the two images. If there are 3 images, each
#' image will be shown in a separate channel. Invisibly returns the RGB image.
#' @param imlst List of images
#' @export
plot_imlist <- function(imlst) {
  stopifnot(length(imlst) >= 2 & length(imlst) < 4)
  if (length(imlst) == 2) imlst[[3]] <- pmin(imlst[[1]], imlst[[2]])
  plot(do.call(EBImage::rgbImage, imlst))
  invisible(do.call(EBImage::rgbImage, imlst))
}

#' Pad each image so that the two are the same size
#'
#' @param im1 image 1
#' @param im2 image 2
#' @export
pad_img_match <- function(im1, im2) {
  img_size <- pmax(dim(im1), dim(im2))
  im1 <- img_pad_to_size(im1, img_size, value = 1)
  im2 <- img_pad_to_size(im2, img_size, value = 1)

  stopifnot(all.equal(dim(im1), dim(im2)))

  return(list(im1, im2))
}

# Apply hanning filter to a list of 2 images.
# If necessary, invert the image first (so white = signal)
hanning_list <- function(lst, invert = T, ...) {
  if (invert) lst <- lapply(lst, function(x) 1 - x)
  lapply(lst, function(x) x * (hanning_2d(dim(x), ...)))
}

# Apply fft to a list
fft_list <- function(lst, shift = T, inverse = F) {
  # Shift first if inverse fft
  if (inverse & shift) lst <- lapply(lst, fftshift)

  fftlst <- lapply(lst, function(x) fft(x, inverse = inverse))
  # Shift last if forward fft
  if (!inverse & shift) fftlst <- lapply(fftlst, fftshift, inv = T)
  fftlst
}

# Get the alignment angle from a set of fft'd images
angle_recover <- function(imlst, nt = 720, ...) {
  # Use the FT magnitude as a new image to get rotation
  im_fft_mag <- lapply(imlst, function(x) Re(sqrt(x * Conj(x))))
  # Convert to polar coords
  im_polar <- lapply(im_fft_mag, function(x) img_to_polar(x, ntheta = nt))
  # Apply hanning filter
  im_polar_hanning <- hanning_list(im_polar, invert = F, ...) # signal is already white
  # register ffts
  do.call(fft_register, im_polar_hanning)
}

# Perform fft-based image registration (simple translation)
fft_register <- function(im1, im2) {
  # assumes hanning window already applied
  # FFT image
  i1 <- fft(im1)
  i2 <- fft(im2)

  cross_power <- i1*Conj(i2)
  normalized_cross_power <- cross_power / Mod(cross_power)


  ifft <- Re(fft(cross_power, inv = T))/length(cross_power)
  # Force same format
  ifft <- 0*im1 + ifft

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
