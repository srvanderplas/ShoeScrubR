#' Get a specific shoe mask
#'
#' @param brand Adidas or Nike
#' @param size Numerical size: 7, 7.5, 10, 10.5 for Adidas, 8, 8.5, 10, 10.5
#'        for Nike
#' @param foot R or L, determined from the image (e.g. if the image is a film
#'        print, this would be the opposite foot from the image name because
#'        film prints are reversed)
#' @param ppi pixels per inch (linearly) for the returned mask image
#' @export
#' @examples
#' shoe_mask("Nike", 10.5, "R") %>% plot()
#' plot(normalize(shoe_mask("Nike", 10, "R") + shoe_mask("Nike", 10.5, "R")))
shoe_mask <- function(brand, size, foot, ppi = 300) {
  stopifnot(brand %in% c("Adidas", "Nike"))
  if (brand == "Adidas") {
    stopifnot(size %in% c(7, 7.5, 10, 10.5))
  } else {
    stopifnot(size %in% c(8, 8.5, 10, 10.5))
  }
  stopifnot(foot %in% c("R", "L"))

  brand_model <- gsub("Adidas", "Adidas_Seeley", brand) %>%
    gsub("Nike", "Nike_Winflow", .)

  size_x <- ifelse(size < 9, paste0(size, "W"), paste0(size, "M"))

  filename <- paste0(paste(brand_model, size_x, foot, sep = "_"), ".png")
  full_filename <- system.file("templates", filename, package = "ShoeScrubR")

  stopifnot(file.exists(full_filename))
  im <- EBImage::readImage(full_filename)
  scale <- ppi/200
  output_dim <- floor(dim(im)*scale)
  EBImage::resize(im, w = output_dim[1], h = output_dim[2])
}

#' Align an image and a mask based on principal components
#'
#' One of img or img_df must be supplied. If img is supplied, additional
#' arguments to img_to_df may also be supplied using ....
#'
#' Adapted from the dudi.pca function in the ade4 package, which performs
#' weighted principal components analysis.
#'
#' @param img_df data frame of the locations of nonzero pixels in an image
#'               (columns row, col, value*). Value is optional, and if supplied,
#'               will be used to weight the results.
#' @param img image
#' @param weighted should weighted calculation be used?
#' @param ... additional arguments to image_to_df
#'
align_prcomp <- function(img = NULL, img_df = NULL, weighted = F, ...) {
  if (!is.null(img) & is.null(img_df)) {
    img_df <- image_to_df(img, ...)
  }
  stopifnot(!is.null(img_df))

  if (hasName(img_df, "value") & weighted) {
    rowmean <- mean(img_df$row, weight = img_df$value, na.rm = T)
    colmean <- mean(img_df$col, weight = img_df$value, na.rm = T)
    img_df$row <- img_df$row - rowmean
    img_df$col <- img_df$col - colmean
    weight <- img_df$value/sum(img_df$value)
  } else {
    weight <- 1/nrow(img_df)
  }

  df <- img_df[,c("row", "col")] %>% as.matrix()

  pca <- prcomp(df*sqrt(weight), center = F, scale = F)
  # angle <- pca$rotation[2,1] %>% acos() %>% `*`(180/pi)
  #
  # angle
}

#' Weighted principal components analysis
#'
#' @param
wpca <- function(df, row.w = rep(1, nrow(df))/nrow(df))
{
  dfx <- df
  row.w <- if (hasName(df, "value")) df$value else rep(1, nrow(df))
  df <- df[,c("row", "col")] %>%
    na.omit()
  center <- apply(df, 2, function(v) sum(v * row.w)/sum(row.w))
  df <- sweep(df, 2, center)

  df <- as.matrix(df)
  df.ori <- df
  df <- df * sqrt(row.w)
  df <- crossprod(df, df)

  eig1 <- eigen(df, symmetric = TRUE)

  acos(1 -eig1$vectors[1,which.max(eig1$values)]) * 180/pi
}

