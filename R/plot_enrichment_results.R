#' A plot function for TEA.
#'
#'
#' @param df dataframe as output by `enrichment_analysis()`.
#' @param labels_on_columns write the labels on top of the bars, instead of on the side. Leave NULL for automatic selection.
#' @param n_bars maximal number of bars to be shown.
#' @param filename if not NULL, file name to save the plot as png.
#'
#'
#' @returns a ggplot object
#'
#' @import ggplot2
#' @export
plot_enrichment_results <- function(df, labels_on_columns = NULL,
                                    n_bars = 15L, filename = NULL) {

  if( ! all(c("Q value", "Term") %in% colnames(df))){
    stop("`df` must be the result of `enrichment_analysis()` and contain columns named 'Q value' and 'Term'")
  }


  df <- df[1:n_bar,]

  df$Term2 <- gsub(" WBbt:\\d{7}$", "", df$Term)
  df$Term2 <- gsub(" WBPhenotype:\\d{7}$", "", df$Term2)
  df$Term2 <- gsub(" GO:\\d{7}$", "", df$Term2)

  df$Term2 <- factor(df$Term2, levels = rev(df$Term2))

  if(is.null(labels_on_columns)){
    labels_on_columns <- (max(nchar(as.character(df$Term2))) > 29)
  }


  gg <- ggplot(df) +
    theme_classic() +
    geom_col(aes(x = Term2,
                 y = -log10(`Q value`))) +
    (if(labels_on_columns){
      geom_text(aes(x = Term2, y = .1, label = Term2),
                hjust = 0,
                color = "white")
    }) +
    (if(labels_on_columns){
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    }) +
    coord_flip() +
    ylab(expression(-log[10](q))) +
    xlab(NULL)

  if(!is.null(filename)){
    if(!(is.character(filename) && nchar(filename) > 0)){
      stop("If filename is not NULL, it should be a valid file name.")
    }
    if(! endsWith(filename, ".png")){
      filename <- paste0(filename, ".png")
    }
    ggsave("filename", gg, device = "png", width = 10, height = 10, units = "cm")
  }

  gg
}

