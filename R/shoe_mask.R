#' Get a specific shoe mask
#'
#' @param brand Adidas or Nike
#' @param size Numerical size: 7, 7.5, 10, 10.5 for Adidas, 8, 8.5, 10, 10.5
#'        for Nike
#' @param foot R or L, determined from the image (e.g. if the image is a film
#'        print, this would be the opposite foot from the image name because
#'        film prints are reversed)
#' @export
#' @example
#' shoe_mask("Nike", 10.5, "R") %>% plot()
#' plot(normalize(shoe_mask("Nike", 10, "R") + shoe_mask("Nike", 10.5, "R")))
shoe_mask <- function(brand, size, foot) {
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
  EBImage::readImage(full_filename)
}