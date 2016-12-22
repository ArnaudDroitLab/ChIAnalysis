#' Determines if certain attributes are coherent within certain categories of nodes.
#'
#' @param chia.obj The base ChIA object on which to assess intra-category coherence.
#' @param coherence.function A function which takes in a ChIA object and returns a 
#'   vector of a 2-level factor of the same length as the number fo regions in the 
#'   ChIA object.
#' @param node.categories A data-frame with a number of rows equal to the number of regions
#'   in chia.obj, where each column represents a set of indices indicating which nodes belong
#'   in a given category.
#' @param m.hyper The m parameter for the hypergeometric test, namely the total number of coherent genes.
#' @param output.file The name of the file where the resulting graph should be saved.
#'
#' @return A list containing the computed data and the generated plot object.
#' 
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' @export
coherence.test <- function(chia.obj, coherence.function, node.categories, m.hyper, output.file=NULL, ...) {
  # Apply the coherence function
  coherence.list = category.apply(chia.obj, coherence.function, node.categories, ...)
  
  # Count the success/failures from each category.
  coherence.counts = lapply(coherence.list, function(x) {
    tmp = c(sum(x==levels(x)[1], na.rm=TRUE), sum(x==levels(x)[2], na.rm=TRUE))
    names(tmp) = levels(x)
    return(tmp)
  })
  
  # Turn the count list into a data-frame.
  coherence.df = data.frame(matrix(unlist(coherence.counts), ncol=2, byrow=TRUE))
  colnames(coherence.df) <- names(coherence.counts[[1]])
  rownames(coherence.df) <- names(node.categories)
  
  # Calculate additional statistics
  coherence.df$Ratio = coherence.df[,2] / coherence.df[,1]
  coherence.df$Diff = coherence.df[,2] - coherence.df[,1]
  
  # Calculate a pValue with a hypergeometric test for each network
  # q = the number of coherent genes in the network
  # m = the number of coherent genes in total
  # n = the number of other genes in total
  # k = the number of genes in the network
  q <- coherence.df[,1] + coherence.df[,2]
  m <- m.hyper
  n <- length(genes(TxDb.Hsapiens.UCSC.hg19.knownGene)) - m
  k <- aggregate(Gene.Representative ~ Component.Id, chia.obj$Regions, sum)$Gene.Representative
  coherence.df$PValue <- phyper(q, m, n, k, lower.tail = FALSE)
  
  # Add the names of genes in a list.
  gene.names = lapply(coherence.list, function(x) {
    level1 <- names(x)[x == levels(x)[1]]
    level1.genes <- paste(level1[!(is.na(level1))], collapse = ", ")
    level2 <- names(x)[x == levels(x)[2]]
    level2.genes <- paste(level2[!(is.na(level2))], collapse = ", ")
    others <- names(x)[is.na(x)]
    no.change <- paste(others[!(is.na(others))], collapse = ", ")
    genes.list <- c(level1.genes, level2.genes, no.change)
    names(genes.list) <- c(paste0(levels(x)[1], ".genes"), paste0(levels(x)[2], ".genes"), "No.changes")
    return(genes.list)
  })
  
  # Add the names of genes in the data frame.
  coherence.df[,(ncol(coherence.df) + 1):(ncol(coherence.df) + 3)] <- data.frame(t(sapply(gene.names, function(x){return(x)})))
  
  # Remove categories which do not have at least two elements.
  coherence.df.subset = coherence.df[coherence.df[,1] + coherence.df[,2] > 1,]
  
  scores = (pmax(coherence.df.subset[,2], coherence.df.subset[,1]))/(coherence.df.subset[,2]+coherence.df.subset[,1])
  scores = ifelse(coherence.df.subset[,2] > coherence.df.subset[,1], scores, -scores)
  coherence.df.subset$Score = scores
  ordering = order(scores, -coherence.df.subset[,1], coherence.df.subset[,2])
  coherence.df.subset$Category = factor(rownames(coherence.df.subset), rownames(coherence.df.subset)[ordering])
  coherence.df.subset$Threshold = ifelse(scores <= -0.75, "Down", ifelse(scores >= 0.75, "Up", "Non-coherent"))
  
  # Generate a plot representing the data.
  ymin = -max(coherence.df.subset[,1])
  ymax = max(coherence.df.subset[,2])
  
  ribbon.df = coherence.df.subset
  ribbon.df$Category = as.numeric(ribbon.df$Category)
  
  plot.obj = ggplot(coherence.df.subset, aes(x=Category)) +
    geom_bar(mapping=aes(y=eval(parse(text=colnames(coherence.df)[2]))), stat="identity", fill="white", color=rgb(0, 0, 255, maxColorValue=255)) +
    geom_bar(mapping=aes(y=eval(parse(text=paste0("-",colnames(coherence.df)[1])))), stat="identity", fill="white", color=rgb(255, 0, 0, maxColorValue=255)) +
    geom_ribbon(data=ribbon.df, mapping=aes(x=Category, fill=Threshold, alpha=Threshold), ymin=ymin, ymax=ymax) +
    scale_fill_manual(values=c(Down=rgb(39, 170, 225, maxColorValue=255), Up=rgb(39, 170, 225, maxColorValue=255), "Non-coherent"=rgb(255, 255, 225, maxColorValue=255))) +
    scale_alpha_manual(values=c(Down=0.3, Up=0.3, "Non-coherent"=0)) +
    labs(y="Gene attribute") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          axis.text = element_text(color="black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  if(!is.null(output.file)) {
    ggsave(output.file, width=14, height=7, plot=plot.obj)
  }
  
  return(list(Plot=plot.obj, Data=coherence.df, Genes=gene.names))
}

#' Given a ChIA object, determine which genes are active and which are not.
#'
#' @param chia.obj The base ChIA object on which to assess gene activity.
#'
#' @return A 2-level vector with values "Inactive", "Active" or NA for each node.
#'
#' @importFrom plyr revalue
#' @export
active.gene.coherence <- function(chia.obj) {
    tmp = ifelse(chia.obj$Regions$Gene.Representative, chia.obj$Regions$Is.Gene.Active, NA)
    tmp = revalue(factor(tmp, levels=c(FALSE, TRUE)), c("FALSE"="Inactive", "TRUE"="Active"))
    names(tmp)[chia.obj$Regions$Gene.Representative] <- chia.obj$Regions$SYMBOL[chia.obj$Regions$Gene.Representative]
    return(tmp)
}

#' Generates a function which takes a chia object, and determine which of its nodes
#' have an positive or negative fold-change based on the given column name.
#'
#' @param fc.column The column the new function should check for fold-change values.
#'
#' @return A function taking a ChIA object, and returning a 2-level vector with 
#'   "Upregulated", "Downregulated" or NA depending on the value of the given 
#'   fold-change column.
#'
#' @export
fold.change.coherence <- function(fc.column) {
    force(fc.column)
    function(chia.obj) {
        values = chia.obj$Regions[,fc.column]
        direction = factor(ifelse(values > 0, "Upregulated", "Downregulated"), levels=c("Downregulated", "Upregulated"))
        
        names(direction)[chia.obj$Regions$Gene.Representative] <- chia.obj$Regions$SYMBOL[chia.obj$Regions$Gene.Representative]
        
        return(direction)
    }
}

#' Categorizes regions according to their expression level.
#'
#' Given expression threshold, this function classifies all regions as 
#' belonging to a gene "below" the bottom expression level or "above" the top
#' expression level.
#'
#' @param chia.obj The chia object whose regions must be classified.
#' @param top The threshold above which a region is classified as "Top".
#' @param bottom The threshold under which a region is classified as "Bottom".
#'
#' @return A vector with a 2-level factor, where regions are classified either
#'   as "Top", "Bottom" or NA according to their expression level.
#' @export
expression.coherence <- function(chia.obj, top, bottom) {
  
  # Create the vector
  expression <- rep(NA, times = nrow(chia.obj$Regions))
  expression[(chia.obj$Regions$Expr.mean <= bottom) & chia.obj$Regions$Gene.Representative] <- "Bottom"
  expression[(chia.obj$Regions$Expr.mean >= top) & chia.obj$Regions$Gene.Representative] <- "Top"
  # Factorize the vector
  expression <- factor(expression, levels = c("Bottom", "Top"))
  # Name the vector with the genes symbols
  #if (names){
  names(expression)[chia.obj$Regions$Gene.Representative] <- chia.obj$Regions$SYMBOL[chia.obj$Regions$Gene.Representative]
  #}
  return(expression)
}