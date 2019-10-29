#' Automatically crop blank space from the edges of an image
#'
#' @param img Image
#' @param bkgd background pixel value (defaults to 1 or 0, whichever is in the top corner)
#' @param pad how much padding to leave?
#' @export
img_autocrop <- function(img, bkgd = NULL, pad = 0) {
  multilayer <- length(dim(img))
  if (multilayer > 2) {

    dim_var <- apply(img, 3, function(x) mean(x %in% c(0,1)))
    dim_bin <- dim_var >= .9
    if (sum(dim_bin) > 0) {
      idx <- which(dim_bin)[1]
      newimg <- img[,,idx]
      colorMode(newimg) <- 0

    } else {
      warning("Color image has no binary layers")
      newimg <- img
    }
  } else {
    newimg <- img
  }
  if (is.null(bkgd)) bkgd <- newimg[1]

  img_dim1 <- apply(newimg != bkgd, 1, sum)
  crop_dim1 <- c(min(which(img_dim1 != 0)), max(which(img_dim1 != 0))) + c(-pad, pad)
  crop_dim1 <- pmin(pmax(1, crop_dim1), dim(img)[1])

  img_dim2 <- apply(newimg != bkgd, 2, sum)
  crop_dim2 <- c(min(which(img_dim2 != 0)), max(which(img_dim2 != 0))) + c(-pad, pad)
  crop_dim2 <- pmin(pmax(1, crop_dim2), dim(img)[2])

  if (length(dim(img)) > 2) {
    res <- EBImage::as.Image(img[crop_dim1[1]:crop_dim1[2], crop_dim2[1]:crop_dim2[2],])
  } else {
    res <- EBImage::as.Image(img[crop_dim1[1]:crop_dim1[2], crop_dim2[1]:crop_dim2[2]])
  }

  attr(res, "operation") <- append(attr(img, "operation"),
                                   list(list(type = "autocrop",
                                             top_bottom = crop_dim2,
                                             left_right = crop_dim1,
                                             value = bkgd)))

  res
}