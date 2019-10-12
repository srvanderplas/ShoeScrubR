library(EBImage)

# Make all templates black-and-white
template_list <- list.files("inst/templates", "(Adidas|Nike).*.png", full.names = T)

bw_template <- function(file) {
  img <- readImage(file, all = F)
  if (length(dim(img)) > 2) {
    img <- img[,,1]
  }


  # colorMode(img) <- 0
  # channel(img, "luminance") %>%
  #   round() %>%
  #   # (function(.) . > mean(.[1:10, 1:10])) %>%
  img %>%
    writeImage(file)
  # writeImage(newim, file)
}

purrr::map(template_list, bw_template)

img <- readImage("inst/templates/Adidas_Seeley_10.5M_L.png", )


# Orient all templates the same

straighten_image <- function(img) {
  s1 <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), nrow = 3, byrow = T)
  s2 <- matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), nrow = 3, byrow = T)

  gx <- round(filter2(img, s1), 4)
  gy <- round(filter2(img, s2), 4)

  edges <- sqrt(gx^2 + gy^2)
  edge_locations <- edges %>%
    image_to_df()


  # im_center <- edge_locations %>% select(-val) %>%
  #   summarize_all(mean)

  slope <- coef(lm(row ~ col, data = edge_locations))[2]

  slope <- round(atan(as.numeric(slope)))*180/pi
  if (abs(slope) > 0) {
    message("Rotating by %0.2f degrees", 90 + slope)
    img_r <- EBImage::rotate(img, 90 + slope, bg.col = 0)
  } else {
    img_r <- img
  }

  img_r
}

imgs <- purrr::map(template_list, readImage)

rotated_imgs <- purrr::map(imgs, straighten_image)
