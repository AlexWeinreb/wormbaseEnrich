#' A plot function for TEA.
#'
#'
#' @param df dataframe as output by `enrichment_analysis()`.
#' @param labels_inside_columns write the labels on top of the bars, instead of on the side. Leave NULL for automatic selection.
#' @param max_bars maximal number of bars to be shown.
#'
#'
#' @returns a ggplot object
#'
#' @import ggplot2
#' @export
plot_enrichment_results <- function(df, labels_inside_columns = NULL,
                                    max_bars = 15L) {

  if( ! all(c("FDR", "term_name") %in% colnames(df))){
    stop("`df` must be the result of `enrichment_analysis()` and contain columns named 'FDR' and 'term_name'")
  }


  if(nrow(df) > max_bars) df <- df[1:max_bars,]

  df$term_name <- factor(df$term_name, levels = rev(df$term_name))

  if(is.null(labels_inside_columns)){
    labels_inside_columns <- (max(nchar(as.character(df$term_name))) > 29)
  }


  ggplot(df) +
    theme_classic() +
    geom_col(aes(x = .data$term_name,
                 y = -log10( .data$FDR ) )) +
    (if(labels_inside_columns){
      geom_text(aes(x = .data$term_name, y = .1, label = .data$term_name),
                hjust = 0,
                color = "white")
    }) +
    (if(labels_inside_columns){
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    }) +
    coord_flip() +
    ylab(expression(-log[10](q))) +
    xlab(NULL)

}

