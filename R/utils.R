#' Image to data frame
#'
#' @param img image
#' @param filter_val value(s) of pixels to remove
#' @export
image_to_df <- function(img, filter_val = 0) {
  tmp <- img %>%
    EBImage::imageData() %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::mutate(row = 1:dplyr::n()) %>%
    magrittr::set_names(str_replace(names(.), "\\.{1,}", "")) %>%
    tidyr::gather(key = col, val = val, -row) %>%
    dplyr::mutate(col = as.numeric(col)) %>%
    dplyr::mutate(row = -row)

  if (is.null(filter_val)) tmp

  tmp %>%
    dplyr::filter(!val %in% filter_val)
}
