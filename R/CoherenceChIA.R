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
  coherence.df$Category = factor(rownames(coherence.df), rownames(coherence.df)[order(coherence.df$Diff)])

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
  
  # Generate a plot representing the data.
  plot.obj = ggplot(coherence.df.subset, aes(x=Category)) +
    geom_bar(mapping=aes(y=eval(parse(text=colnames(coherence.df)[2]))), stat="identity", fill="Blue", color="black") +
    geom_bar(mapping=aes(y=eval(parse(text=paste0("-",colnames(coherence.df)[1])))), stat="identity", fill="Red", color="black") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
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

#' Generates a fucntion which takes a chia object, and determine which of its nodes
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

#' Generates a fucntion which takes a chia object, and determine which of its nodes
#' are in the top 25% and in the bottom 25% of expression.
#'
#' @param chia.obj The base ChIA object on which to assess gene activity.
#' @param top The expression value limit of the top 25% expressed genes.
#' @param bottom The expression value limit of the bottom 25% expressed genes.
#'
#' @return A 2-level vector with 
#'   "Top", "Bottom" or NA depending on the value of the given 
#'   fold-change column.
#'
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