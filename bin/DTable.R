#' Title
#'
#' @param df 
#' @param nDigits_after_decimal 
#' @param table_caption 
#' @param escape 
#' @param font_family 
#'
#' @return
#' @export
#'
#' @examples
DTable <- function(df = df,
                   nDigits_after_decimal = 1,
                   table_caption = ""  ,
                   escape = FALSE,
                   font_family = "sans-serif"){ # absolut path to file
  
  
  # > READING FILE INTO DF; DF TRANSFORMATIONS:
  # Read file
  
  # > Customize interactive DT::datatable
  DT::datatable(df, 
                rownames   = FALSE,
                escape     = FALSE,
                fillContainer = FALSE,
                filter     = "bottom",
                caption    = table_caption,
                extensions = c('FixedColumns','Scroller', 'Buttons'),
                # OPTIONS:
                options = list(
                  
                  # Does not allow columnful dataframes go rogue and tucks them in to fit page width
                  scrollX = TRUE,
                  
                  # Defines all capabilities
                  dom        = 'PBRMDCT<"clear">lfrtip',
                  
                  autoWidth  = TRUE,
                  ColReorder = TRUE,
  
                  #   columnDefs = list(list(targets = length(colnames(df)), visible = TRUE)))),
                  
                  lengthMenu = list(c(10, 50, -1), c('10', '50', 'All')),

                  buttons    = list('copy','print', 
                                    list(extend  = 'collection',
                                         buttons = c('csv', 'excel', 'pdf'),
                                         text    = 'Export')),
                  # Black header container for colnames
                  initComplete = JS(
                    "function(settings, json) {",
                    "$(this.api().table().header()).css({'color': '#fff','font-family': 'sans-serif', 'background-color': '#4e4b4c'});",
                    "$(this.api().table().body()).css({'color': '#4e4b4c','font-family': 'sans-serif',   'text-align' : 'center'});",
                    "$(this.api().table().footer()).css({'color': '#fff','font-family': 'sans-serif'});",
                    "$(this.api().table().container()).css({'color': '#fff','font-family': 'sans-serif', 'outline-color' : '#4e4b4c' });",
                    "$(this.api().table().node()).css({'color': '#fff','font-family': 'sans-serif'});",
                    "}") )) %>% 
    
    # Change fontsize of cell values
    formatStyle(columns    = seq_along(colnames(df)), 
                fontSize   = "85%",
                fontFamily = "sans-serif")   %>%

    formatSignif(
      columns = unlist(lapply(df, is.numeric)),
      digits = 2,
      interval = 3,
      mark = ",",
      dec.mark = getOption("OutDec"))                 -> fancyDatatable

  
  return(fancyDatatable)  
}
