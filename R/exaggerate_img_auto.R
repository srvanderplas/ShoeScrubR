#' Exaggerate an image to a mask-like appearance automatically
#'
#' @param img Image
#' @param ... extra arguments, primarily to em_thresh
exaggerate_img_auto <- function(img, ...) {

  img_masked <- img_em_clean(img, ...)

  img_blur_mask <- img_clean_blur(img_masked)

  img_mask_clean(img_blur_mask)
}

img_em_clean <- function(img, ...) {
  if (is.list(img)) {
    return(lapply(img, img_em_clean, ...))
  }

  emt <- em_thresh(img, ...)

  stopifnot(length(emt$img_ratios) >= 1)

  labeled_img <- emt$img_ratios[[1]] %>% # EM algorithm to threshold
    (function(.) { . > 10 }) %>% # binarization
    EBImage::as.Image() %>%
    img_auto_clean() %>%
    EBImage::bwlabel() %>%
    clean_img_corners() # This takes ~9 seconds for a regular-size film scan
  # and ~22 seconds for a larger film scan
  # The time difference is due to the auto-cleaning

  img * (labeled_img > 0) + (labeled_img == 0)

}

#' Compute likely pixel category based on EM algorithm clustering of intensity
#'
#' EM algorithm is for normally distributed groups; this assumption is likely
#' not accurate, but it is fast and effective.
#' If both N and scale_factor are NULL, the full image will be used, which will be slow.
#' @param img Image
#' @param N Number of points to sample from the image (speeds up computational time)
#' @param scale_factor Alternative to N - scales image by a factor of scale_factor
#' @param ngroups Number of clusters
#' @param quiet suppress output from normalmixEM using sink()?
#' @param ... additional arguments to mixtools::normalmixEM
#' @importFrom mixtools normalmixEM
#' @importFrom stats dnorm
#' @importFrom abind abind
#' @importFrom EBImage as.Image
em_thresh <- function(img,
                      scale_factor = 10,
                      N = ifelse(is.null(scale_factor), pmin(length(img), 30000), NULL),
                      ngroups = 3, quiet = T, ...) {
  if (is.list(img)) {
    return(lapply(img, em_thresh, ...))
  }

  # imsmall <- sample(as.numeric(img), size = N, replace = F)
  if (!is.null(scale_factor)) {
    imsmall <- img_resize(img, w = floor(dim(img)[1]/scale_factor),
                          h = floor(dim(img)[2]/scale_factor))
  } else if (!is.null(N)) {
    new_dim <- round(dim(img) * sqrt(N/length(img)))
    imsmall <- img_resize(img, w = new_dim[1], h = new_dim[2])
  } else {
    imsmall <- img
  }

  if (quiet) {
    sink("/dev/null")
  }
  em <- try(mixtools::normalmixEM(imsmall, k = ngroups, ...), silent = T)
  if (quiet) {
    sink()
  }

  while (ngroups > 2 & "try-error" %in% class(em)) {
    ngroups <- ngroups - 1
    message("EM mixture did not succeed with ", ngroups + 1,
            " groups. Retrying with ", ngroups, " groups.")
    if (quiet) {
      sink("/dev/null")
    }
    em <- try(mixtools::normalmixEM(imsmall, k = ngroups, ...), silent = T)
    if (quiet) {
      sink()
    }
  }

  if ("try-error" %in% class(em)) stop("Mixture model fitting failed\n", em)

  values <- cbind(em$x, em$posterior)

  mean_idx <- order(em$mu, decreasing = F)

  # get probs for full img
  img_dens_functions <- lapply(mean_idx, function(x) {
    0*img + dnorm(img, mean = em$mu[x], sd = em$sigma[x])
  })

  img_dens_total <- abind::abind(img_dens_functions, along = 3) %>%
    apply(c(1, 2), sum)

  img_dens_ratio <- lapply(img_dens_functions, function(x) {
    EBImage::as.Image(x/(img_dens_total - x))
  })

  return(list(img_ratios = img_dens_ratio, em = em))
}

img_auto_clean <- function(img) {
  ## automatic in the sense that it doesn't require tuning parameters to
  ## image size/resolution
  img %>%
    EBImage::erode(EBImage::makeBrush(size = 1)) %>%
    EBImage::dilate(EBImage::makeBrush(size = 3, shape = "line")) %>%
    EBImage::erode(EBImage::makeBrush(size = 5, shape = "disc")) %>%
    EBImage::dilate(EBImage::makeBrush(size = 9, shape = "disc"))
}


clean_img_corners <- function(labeled_img, len = NULL) {
  # TODO: Add checks for image being labeled...

  if (is.null(len)) len <- pmax(5, floor(sqrt(length(labeled_img))/40))

  idx_rows <- list(
    unique(pmin(pmax(1:len, 1), nrow(labeled_img))),
    unique(pmin(pmax((nrow(labeled_img) - len):nrow(labeled_img), 1),
                nrow(labeled_img))))
  idx_cols <- list(
    unique(pmin(pmax(1:len, 1), ncol(labeled_img))),
    unique(pmin(pmax((ncol(labeled_img) - len):ncol(labeled_img), 1),
                ncol(labeled_img))))
  idxs <- expand.grid(row = idx_rows, col = idx_cols)

  imgtbls <- lapply(1:4, FUN = function(i) {
    labeled_img[idxs$row[[i]], idxs$col[[i]]]
  }) %>%
    unlist() %>%
    unique()

  pct_label <- data.frame(
    label = imgtbls,
    pct = sapply(imgtbls, function(x) mean(labeled_img == x)))

  pct_label <- pct_label[pct_label$pct < .1, ]

  if (length(pct_label) > 0) labeled_img[labeled_img %in% pct_label$label] <- 0

  labeled_img
}

img_clean_blur <- function(img) {
  if (is.list(img)) {
    return(lapply(img, img_clean_blur))
  }

  img_blur_labels <- img %>%
    EBImage::gblur(sigma = sqrt(length(img))/50, boundary = "replicate") %>%
    (function(.) . < median(.)) %>%
    EBImage::bwlabel() %>%
    EBImage::fillHull()

  big_label <- table(img_blur_labels[img_blur_labels != 0]) %>%
    sort(decreasing = T) %>%
    names() %>% as.numeric()

  img_blur_labels == big_label[1]
}


img_mask_clean <- function(img) {
  if (is.list(img)) {
    return(lapply(img, img_mask_clean))
  }

  brushsize <- floor(sqrt(length(img))/5)
  brushsize <- brushsize + as.numeric(brushsize %% 2 == 0) # ensure brush size is odd

  EBImage::fillHull(img) %>%
    EBImage::opening(EBImage::makeBrush(brushsize, "disc"))
}
