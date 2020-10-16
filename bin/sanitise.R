#' Title Sanitise LaTeX formatted url table entries
#'
#' @param table
#' @param col_names
#' @param pattern
#'
#' @return
#' @export
#'
#' @examples
sanitise <- function(table     = NULL,
                     col_names = FALSE,
                     pattern   = "href") {

if (col_names == FALSE) {
  col_names <- colnames(table)
  }

for (col in col_names) {
  col_as_string <- paste(table[[col]], collapse = "")
  if (grepl(pattern, col_as_string, fixed = TRUE)) {

  url <- paste0(col, "_url")
  table[[url]] <- gsub("\\href\\{", "", table[[col]])
  table[[col]] <- sub(".*\\{(.*)\\}.*", "\\1", table[[col]])
  table[[url]] <- gsub("\\\\", "", table[[url]])
  table[[url]] <- gsub("\\}.*", "", table[[url]])
  table[[col]] <- gsub("\\href\\{", "", table[[col]])

  # prepare url for DT::datatable
  table[[col]] <- paste0("<a href='",
                      table[[url]],
                      "'",
                      "target='blank",
                      "'>",
                      table[[col]],
                      "</a>") 
  as.data.frame(table) %>% 
    dplyr::mutate_all(stringr::str_replace_all,"<a href=''target='blank'></a>", "") -> table
    }
  }
return(table)
}