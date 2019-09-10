library(EBImage)

# Make all templates black-and-white
template_list <- list.files("inst/templates", "(Adidas|Nike).*.png", full.names = T)

bw_template <- function(file) {
  readImage(file) %>%
    channel("Red") %>%
    writeImage(file)
}

purrr::map(template_list, bw_template)
template <- readImage("inst/Adidas_Seeley_10.5M_L.png") %>%
  channel('Red')


# Orient all templates the same

straighten_image <- function(img) {
  s1 <- matrix(c(-1, 0, 1, -2, 0, 2, -1, 0, 1), nrow = 3, byrow = T)
  s2 <- matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), nrow = 3, byrow = T)

  gx <- round(filter2(img, s1), 4)
  gy <- round(filter2(img, s2), 4)

  edges <- sqrt(gx^2 + gy^2)
  edge_locations <- edges


  slope <- coef(lm(row ~ col, data = edge_locations))[2]

  if (round(atan(slope)) > 0)
  im_center <- edge_locations %>% select(-val) %>%
    summarize_all(mean)
}