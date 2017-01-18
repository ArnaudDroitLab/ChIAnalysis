#' Return the left part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "left side".
#' of the original data.
#' @export
chia_left <- function(chia.obj) {
    return(chia.obj$Regions[as_edgelist(chia.obj$Graph)[,1],])
}

#' Return the right part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "right side".
#' of the original data.
#' @export
chia_right <- function(chia.obj) {
    return(chia.obj$Regions[as_edgelist(chia.obj$Graph)[,2],])
}

#' Return a data-frame representing all edges and the annotation of the 
#' vertices it connects.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by 
#'   \code{\link{load_chia}}.
#'
#' @return A data-frame representing all edges of the graph and their annotation.
#' @export
get_annotated_edges <- function(chia.obj) {
    left.df = as.data.frame(chia_left(chia.obj))
    colnames(left.df) <- paste("Left", colnames(left.df), sep=".")
    right.df = as.data.frame(chia_right(chia.obj))
    colnames(right.df) <- paste("Right", colnames(right.df), sep=".")
    return(data.frame(left.df, right.df))
}

#' Return the right part of the ChIA-PET data.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A \linkS4class{GRanges} object with the \code{Regions} from the "\code{chia.obj}" parameter corresponding the the "right side".
#' of the original data.
#' @export
get_granges <- function(chia.obj) {
    return(GRanges(chia.obj$Regions))
}

#' Determines if the given chia.obj has associated chromatin states.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated chromatin states.
#' @export
has_chrom_state <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Chrom.State))
}

#' Determines if the given chia.obj has associated components.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated components.
#' @export
has_components <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Component.Id) && !is.null(chia.obj$Regions$Component.size))
}

#' Determines if the given chia.obj has associated gene specificities.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated gene specificities.
#' @export
has_gene_specificity <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative) &&
           !is.null(chia.obj$Regions$Expression.Tau) &&
           !is.null(chia.obj$Regions$Expression.Category))
}

#' Determines if the given chia.obj has associated node degrees.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated node degrees.
#' @export
has_degree <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Degree))
}

#' Determines if the given chia.obj has associated gene representatives.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated gene representatives.
#' @export
has_gene_representative <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative))
}

#' Determines if the given chia.obj has associated expression levels.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated expression levels.
#' @export
has_expression_levels <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Gene.Representative) &&
           !is.null(chia.obj$Regions$Expr.mean))
}

#' Determines if the given chia.obj has associated gene annotations.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated gene annotations.
#' @export
has_gene_annotation <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Simple.annotation))
}

#' Determines if the given chia.obj has associated node centralities.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has associated node centralities.
#' @export
has_centrality <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Is.central))
}

#' Determines if the given chia.obj has TF binding information.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has TF binding information.
#' @export
has_transcription_factors <- function(chia.obj) {
    return(sum(grepl("^TF", colnames(chia.obj$Regions))) > 0)
}

#' Determines if the given chia.obj has polymerase binding information.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has polymerase binding information.
#' @export
has_polymerases <- function(chia.obj) {
    return(sum(grepl("^POL", colnames(chia.obj$Regions))) > 0)
}

#' Determines if the given chia.obj has_histone_marks binding information.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has_histone_marks binding information.
#' @export
has_histone_marks <- function(chia.obj) {
    return(sum(grepl("^HIST", colnames(chia.obj$Regions))) > 0)
}

#' Determines if the given chia.obj has TF binding information.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return True if the object has TF binding information.
#' @export
has_fitness <- function(chia.obj) {
    return(!is.null(chia.obj$Regions$Fitness))
}

#' Return the number_of_nodes in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The number_of_nodes in the chia object.
#' @export
number_of_nodes <- function(chia.obj) {
    node.count = vcount(chia.obj$Graph)
    stopifnot(nrow(chia.obj$Regions) == node.count)

    return(vcount(chia.obj$Graph))
}

#' Return the number_of_contacts in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The number_of_contacts in the chia object.
#' @export
number_of_contacts <- function(chia.obj) {
    return(ecount(chia.obj$Graph))
}

#' Return the number_of_components in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The number_of_components in the chia object.
#' @export
number_of_components <- function(chia.obj) {
    return(components(chia.obj$Graph)$no)
}

#' Return the mean component size of a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The mean component size of the chia object.
#' @export
average_component_size <- function(chia.obj) {
    return(mean(components(chia.obj$Graph)$csize))
}

#' Return the number_of_genes represented in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The number_of_genes represented in the chia object.
#' @export
number_of_genes <- function(chia.obj) {
    return(sum(chia.obj$Regions$Gene.Representative))
}

#' Return the number of active genes represented in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The number of active genes represented in the chia object.
#' @export
number_active_genes <- function(chia.obj) {
    return(sum(chia.obj$Regions$Gene.Representative & chia.obj$Regions$Is.Gene.Active))
}

#' Return the number genes per component in the CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The number of active genes represented in the chia object.
#' @export
genes_by_component <- function(chia.obj) {
    stopifnot(has_gene_representative(chia.obj) && has_components(chia.obj))
    return(aggregate(Gene.Representative~Component.Id, as.data.frame(chia.obj$Regions), sum))
}

#' Return the proportion of nodes representing genes in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The number of active genes represented in the chia object.
#' @export
proportion_genes <- function(chia.obj) {
    return(number_of_genes(chia.obj) / number_of_nodes(chia.obj))
}

#' Return the proportion of active genes among all genes in a CHIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return The proportion of active genes among all genes in the chia object.
#' @export
proportion_active_genes <- function(chia.obj) {
    return(number_active_genes(chia.obj) / number_of_genes(chia.obj))
}

#' Returns a list of all statistics for a given ChIA object.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A named list of calculated properties for the given ChIA object.
#' @export
get_all_statistics <- function(chia.obj) {
    results = lapply(list("Number of nodes"            = number_of_nodes,
                          "Number of contacts"         = number_of_contacts,
                          "Number of components"       = number_of_components,
                          "Average component size"     = average_component_size,
                          "Number of genes"            = number_of_genes,
                          "Number of active genes"     = number_active_genes,
                          "Proportion of genes"        = proportion_genes,
                          "Proportion of active genes" = proportion_active_genes),
                    function(f) f(chia.obj))

    return(results)
}

#' Returns the shortest_distances between a set of nodes and another one.
#'
#' @param chia.obj The ChIA object on which distances must be computed.
#' @param from A vector of indices describing the starting nodes in the distance calculations.
#' @param to A vector of indices describing the destination nodes in the distance calculations.
#'
#' @return A vector containing the distance from all "from" nodes to all "to" nodes.
#' @export
shortest_distance <- function(chia.obj, from, to) {
  dist.mat = distances(chia.obj$Graph, v=which(from), to=which(to))
  
  # Get the shortest_distance for all 'from' nodes to any 'to' nodes.
  shortest.dist = apply(dist.mat, 1, min)
  
  return(shortest.dist)
}

#' Returns a summary of the distances between two sets of nodes.
#'
#' @param chia.obj The ChIA object on which distances must be computed.
#' @param from A vector of indices describing the starting nodes in the distance calculations.
#' @param to A vector of indices describing the destination nodes in the distance calculations.
#' @param max The maximum distance to be reported.
#' @param values Report absolute values instead of proportions.
#'
#' @return A vector describing the number of/proportion of nodes from 'from'
#'    with the corresponding distance to nodes in 'to'.
#' @export
#' @import igraph
summarize_distances <- function(chia.obj, from, to, max=3, values=FALSE) {
  # Get the matrix of pairwise distances.
  dist.mat = distances(chia.obj$Graph, v=which(from), to=which(to))
  
  # Get the shortest_distance for all 'from' nodes to any 'to' nodes.
  shortest.dist = shortest_distance(chia.obj, from, to)
  
  # Create table of distancs.
  dist.table = table(shortest.dist)
  
  # Make sure all entries up to max are present.
  missing.entries = setdiff(0:max, names(dist.table))
  missing.values = rep(0, length(missing.entries))
  names(missing.values) = missing.entries
  dist.table = c(dist.table, missing.values)
  
  # Sort back entries.
  dist.table = dist.table[order(as.integer(names(dist.table)))]
  
  # Collapse all entries greater than or equal to max.
  collapsed = sum(dist.table[-(1:max)])
  dist.table = dist.table[1:max]
  dist.table = c(dist.table, collapsed)
  names(dist.table)[length(dist.table)] = paste0(max, '+')
  
  if(values) {
    return(dist.table)
  } else {
    return(dist.table / sum(dist.table))
  }
}

#' Determines which regions of a ChIA object are part of gene networks.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#' @param min.gene The minimum number of gene in a component for it to be considered
#'   a gene network.
#'
#' @return A logical vector indicating whcih regiosn are part of a gene network.
#' @export
regions_in_gene_network <- function(chia.obj, min.gene=2) {
    stopifnot(has_components(chia.obj))
    
    component.summary = aggregate(chia.obj$Regions$Gene.Representative, 
                                  by=list(Component=chia.obj$Regions$Component.Id),
                                  FUN=sum)
    return(chia.obj$Regions$Component.Id %in% component.summary$Component[component.summary$x >= min.gene])
}

#' Returns a logical vector indicating which nodes are gene representatives.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#'
#' @return A logical vector indicating which nodes are gene representatives.
#' @export
is_gene_representative <- function(chia.obj) {
    return(chia.obj$Regions$Gene.Representative)
}

#' Returns the region annotations for the representative regions for the given genes.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#' @param genenames The genes whose representative annotation must be returned.
#'
#' @return A data-frame containign the annotation for the selected genes.
#' @export
get_gene_annotation <- function(chia.obj, genenames) {
    return(chia.obj$Regions[is_gene_representative(chia.obj) & chia.obj$Regions$SYMBOL %in% genenames,])
}