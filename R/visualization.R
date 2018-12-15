#' @title Plots the output from \code{calculate_ACE()}
#'
#' @description Plots the output from \code{calculate_ACE()} which contains the A.C.E. model of heritabily (with confidence intervals).
#'
#' @param df A \code{data.frame} that is the output from \code{\link[mherit]{calculate_ACE}}.
#'
#' @param cutoff The minimum value of \code{A} to display. Defaults to the cutoff at which the top 25\% of taxonomic features are displayed.
#'
#' @param color The taxonomic level to color the plot by. Defaults to the first \code{\link[phyloseq]{rank_names}}.
#'
#' @param facet Logical. If TRUE, the output plot will be faceted by taxonomic rank.
#'
#' @param verbose Logical. If FALSE, messages will be supressed.
#'
#' @return This function returns a \code{ggplot} plot of the taxnomic features sorted by their heritability and optionally faceted.
#'
#' @examples
#' p <- plot_heritability(df, cutoff = 0.3, color = "Class", facet = TRUE, verbose = TRUE)
#' p
#'
plot_heritability <- function(df, cutoff, color, facet = TRUE, verbose = TRUE){

  # Check Inputs ----------

  if(is.null(df) | class(df)[1] != "data.frame"){
    message("The data.frame not found or is not an object of class 'data.frame'.")
  }

  a_cols <- c("A_Estimate", "A_2.5%", "A_97.5%")
  if(sum(!colnames(df) %in% a_cols) == 3){
    message("The data.frame does not contain the expected colnames(). Try running `calculate_ACE()` again. ")
  }

  # Sort and Trim Input data.frame ----------

  # remove any instances were the calculations failed
  df$`A_2.5%`[is.nan(df$`A_2.5%`)] <- NA
  df$`A_97.5%`[is.nan(df$`A_97.5%`)] <- NA
  # remove taxa were CI overlaps with 0
  df <- df[df$`A_2.5%` >= 0,]
  # remove rows containing NA
  df <- df[complete.cases(df[,c(a_cols)]),]
  # reorder data.frame
  df <- df[with(df, order(-A_Estimate)), ]
  # replace underscores names
  df$names <- gsub("_", " ", rownames(df))
  # relevel names so they display decreasing top to bottom
  df$names  <- factor(df$names , levels =rev(df$names))

  # color
  if(is.null(color)){
    color <- rank_names(ps)[1]
    if(verbose == TRUE){message(paste0("`color` variable not found. Defaulting to ", color))}
  }

  # cutoff
  if(is.null(cutoff)){
    cutoff <- round(df[ceiling(nrow(df) * .25), "A_Estimate"], 2)
    if(verbose == TRUE){message(paste0("`cutoff` variable not found. Defaulting to ", cutoff, ", which will display the top 25% more heritible taxonomic features."))}
  }
  df <- df[df$A_Estimate >= cutoff,]

  # Report Plotted Details ----------

  if(verbose == TRUE){
    ranks_plural <- c("ASVs", "Genera", "Families", "Orders", "Classes", "Phyla")
    names(ranks_plural) <- c("ASV", "Genus", "Family", "Order", "Class", "Phylum")
    rank_table <- table(df$rank)
    rank_table <- rank_table[rank_table != 0]
    names(rank_table) <- ranks_plural[names(rank_table)]
    message(paste0("Plot output contains ", paste(rank_table,names(rank_table), collapse = ", "), "."))
  }

  # Plot! ----------

  p <- ggplot(df)
  p <- p + geom_hline(yintercept = cutoff, linetype ="dashed", size = .5)
  p <- p + geom_point(aes_string(x = "names", y = "A_Estimate", color = color), size = 2)
  p <- p + geom_errorbar(aes_string(x = "names", ymin="`A_2.5%`", ymax="`A_97.5%`", color = color), width=.5, size=1)
  p <- p + ylab("Heritability (A)")
  p <- p + xlab("Taxonomic Feature")
  p <- p + coord_flip()
  p <- p + guides(color=guide_legend(ncol=1))
  p <- p + theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    legend.background = element_rect(fill = "transparent",colour = NA),
    legend.key = element_rect(fill = "transparent",colour = NA))

  if(facet == TRUE){
    p <- p + facet_grid(~rank)
  }

  #p

  return(p)
}

#

