#' Given two sets of nodes, determine if one shows enrichment.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#' @param hit.functor A functor selecting the "white" elements.
#' @param draw.functor A functor selecting the "drawn" elements.
#'
#' @return A vector with enrichment properties.
#' @export
node_enrichment <- function(chia.obj, hit.functor, draw.functor) {
    # How many total nodes?
    total = number_of_nodes(chia.obj)
    
    # How many are drawn?
    selected.subset = draw.functor(chia.obj)
    selected.number = number_of_nodes(selected.subset)
    
    # How many have the desired property in the complete set?
    hit.all = number_of_nodes(hit.functor(chia.obj))
    
    # How many have the desired property in the selected set?
    hit.selected = number_of_nodes(hit.functor(selected.subset))
    
    # How many hits would we expect by chance?
    hit.expected = (hit.all/total) * selected.number
    
    # What is the enrichment?
    hit.enrichment = hit.selected / hit.expected
    
    # Calculate enrichment
    enrich = 1 - phyper(hit.selected, hit.all, total - hit.all, selected.number)
    
    return(c(Total        = total, 
             Selected     = selected.number, 
             Hit          = hit.all, 
             Hit.selected = hit.selected, 
             Expected     = hit.expected, 
             Fold.change  = log2(hit.enrichment),
             p.value      = enrich))
}

#' Given a data-frame combining the results of multiple calls to node_enrichment, plot those enrichments.
#'
#' @param enrichment.df A dataframe combining the results of multiple node_enrichment calls.
#' @param filename The name of the file where the plot should be saved.
#' @param label A label for the categories represented by enrichment.df's rows.
#' @export
chia_plot_enrichment <- function(enrichment.df, filename, label="Category") {
  # Reorder enrichment by fold-change.
  categories = rownames(enrichment.df)
  enrichment.df$Category = factor(categories, categories[order(enrichment.df$Fold.change)])

  # Determine which enrichments are significative.
  enrichment.df$Significative = factor(ifelse(enrichment.df$p.value < 0.05, "Significant", "Not significant"), levels=c("Not significant", "Significant"))
  
  # Plot.
  ggplot(enrichment.df) +
    geom_bar(aes(x=Category, y=Fold.change, fill = Significative), colour = "black", stat = "identity") +
    ylab("log2(Hits/Expected hits)") + xlab(label) +
    ggtitle("Transcription factor enrichment") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size = unit(1, "cm")) +
    coord_flip()
    
  ggsave(filename, width=7, height=14)
}

#' Calculate enrichment of transcription factors on a selection of nodes.
#'
#' @param chia.obj The network on which TF enrichment needs to be performed.
#' @param draw.function A function which, given a ChIA object, returns a subset of the "drawn" nodes.
#' @return A data-frame giving the enrichments of all transcription factors at the specified nodes.
#' @export
tf_enrichment <- function(chia.obj, draw.function) {
  results = list()
  
  # Loop over all transcription factors, calculating enrichment for each.
  for(tf in get_tf_names(chia.obj)) {
      # Define a function for selecting the nodes bearing that factor.
      tf.select = function(x) {return(select_by_tf(x, tf))}
      
      # Perform enrichment of the TF.
      results[[tf]] = node_enrichment(chia.obj, tf.select, draw.function)
  }
  
  # Combine the results into a data-frame, and rename it appropriately.
  retval = do.call(rbind.data.frame, results)
  colnames(retval) <- names(results[[1]])

  return(retval)
}
