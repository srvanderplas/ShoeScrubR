library(ShoeScrubR)

full_imglist <- list.files("/lss/research/csafe-shoeprints/ShoeImagingPermanent/",
                      pattern = "0\\d{5}[RL]_\\d{8}_5_1_1", full.names = T)

dir <- tempdir()

file.copy(full_imglist, file.path(dir, basename(full_imglist)))
imglist <- list.files(dir, "\\..*", full.names = T)

shoe_info <- read.csv("~/Projects/CSAFE/2018_Longitudinal_Shoe_Project/Clean_Data/shoe-info.csv")

shoe_images <- shoe_info %>%
  group_by(Brand, Size) %>%
  filter(ShoeID == min(ShoeID))

template_imgs <- tibble(
  file = full_imglist,
  shoe_id = str_extract(file, "\\d{6}[RL]"),
  date = str_extract(file, "\\d{8}") %>% lubridate::ymd()
) %>%
  mutate(shoe_id_num = str_extract(shoe_id, "^\\d{3}")) %>%
  filter(shoe_id_num %in% sprintf("%03d", shoe_images$ShoeID)) %>%
  mutate(shoe = str_extract(shoe_id, "[RL]")) %>%
  group_by(shoe_id_num, shoe) %>%
  filter(date == min(date)) %>%
  filter(row_number() == 1)

file.copy(template_imgs$file, file.path("extra", "template_source_images", basename(template_imgs$file)))




imgs <- tibble(
  file = imglist,
  shoe_id = str_extract(file, "\\d{6}[RL]"),
  date = str_extract(file, "\\d{8}") %>% lubridate::ymd()
) %>%
  group_by(shoe_id, date) %>%
  filter(row_number() == 1) %>%
  filter(str_detect(shoe_id, "L")) %>%
  group_by(shoe_id) %>%
  arrange(date) %>%
  mutate(orig = purrr::map(file, readImage, all = F))

imgs <- imgs %>%
  mutate(norm = purrr::map(orig, ~normalize(.) %>% channel('luminance'))) %>%
  mutate(thresh_init = purrr::map(norm, film_mask)) %>%
  mutate(maskborders = purrr::map(thresh_init, mask_borders, d = 80, fill = 0)) %>%
  mutate(mask = purrr::map(maskborders, film_mask_clean, d2 = 91, prop = .05))

imgs <- imgs %>%
  mutate(maskborders =  purrr::map(mask, mask_borders, d = 80, fill = 0) %>%
           purrr::map(function(x) dilate(x, makeBrush(80, "box"))))

imgs <- imgs %>%
  mutate(hull = purrr::map(mask, fillHull)) %>%
  mutate(full_mask = purrr::map(hull, film_expand_mask))

imgs <- imgs %>%
  mutate(cleaned = purrr::map2(norm, full_mask, function(x, y) as.Image(x*y + (1 - y)*median(x)))) %>%
  mutate(cleaned_thresh = purrr::map(cleaned, ~1 - thresh(1 - ., w = 150, h = 150, offset = .05)))
