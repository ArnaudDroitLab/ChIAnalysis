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

#' Obtain the matrix of TF binding.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A matrix indicating which TF binds to which regions.
#' @export
get_tf <- function(chia.obj) {
    stopifnot(has_transcription_factors(chia.obj))
    retval = chia.obj$Regions[,grepl("^TF", colnames(chia.obj$Regions)), drop=FALSE]
    colnames(retval) = get_tf_names(chia.obj)
    return(retval)
}

#' Obtain the name of transcription factors for which annotation is available.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A vector containing the names of the transcription factors for which annotation is available.
#' @export
get_tf_names <- function(chia.obj) {
    stopifnot(has_transcription_factors(chia.obj))
    all.cols = colnames(chia.obj$Regions)
    tf.cols = all.cols[grepl("^TF", all.cols)]
    return(gsub("TF.overlap.", "", tf.cols))
}

#' Obtain the matrix of polymerase binding.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A matrix indicating which polymerases binds to which regions.
#' @export
get_polymerases <- function(chia.obj, overlap.threshold=0) {
    return(ifelse(get_polymerases_percent(chia.obj) > overlap.threshold, 1, 0))
}

#' Obtain the matrix of polymerase binding as a percentage of the ChIA regions.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A matrix indicating which polymerases binds to which regions.
#' @export
get_polymerases_percent <- function(chia.obj) {
    stopifnot(has_polymerases(chia.obj))
    retval = chia.obj$Regions[,grepl("^POL", colnames(chia.obj$Regions)), drop=FALSE]
    colnames(retval) = get_polymerases_names(chia.obj)
    return(retval)
}

#' Obtain the names of polymerases for which annotation is available.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A vector containing the names of the polymerases for which annotation is available.
#' @export
get_polymerases_names <- function(chia.obj) {
    all.cols = colnames(chia.obj$Regions)
    tf.cols = all.cols[grepl("^POL", all.cols)]
    return(gsub("^POL.", "", tf.cols))
}

#' Obtain the matrix of histone marks binding.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A matrix indicating which histone marks binds to which regions.
#' @export
get_histones <- function(chia.obj, overlap.threshold=0) {
    return(ifelse(get_histones_percent(chia.obj) > overlap.threshold, 1, 0))
}

#' Obtain the matrix of histone marks binding as a percentage of the ChIA regions..
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A matrix indicating which histone marks binds to which regions.
#' @export
get_histones_percent <- function(chia.obj) {
    stopifnot(has_histone_marks(chia.obj))
    retval = chia.obj$Regions[,grepl("^HIST", colnames(chia.obj$Regions)), drop=FALSE]
    colnames(retval) = get_histones_names(chia.obj)
    return(retval)
}

#' Obtain the names of histone marks for which annotation is available.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A vector containing the names of the histone marks for which annotation is available.
#' @export
get_histones_names <- function(chia.obj) {
    all.cols = colnames(chia.obj$Regions)
    tf.cols = all.cols[grepl("^HIST", all.cols)]
    return(gsub("^HIST.", "", tf.cols))
}

#' Obtain a matrix of all possible ChIP data sets: transcription factors, polymerase-binding
#' and histone marks.
#'
#' Combines the results of get_tf, get_polymerases and get_histones.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A matrix describing which epigenetic marks appears on which regions.
get.chips <- function(chia.obj, overlap.threshold=0) {
    return(cbind(get_tf(chia.obj), get_polymerases(chia.obj, overlap.threshold), get_histones(chia.obj, overlap.threshold)))
}

#' Obtain the names of ChIP data sets for which annotation is available.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A vector containing the names of the ChIP data sets for which annotation is available.
#' @export
get_chips_names <- function(chia.obj) {
    return(c(get_tf_names(chia.obj), get_polymerases_names(chia.obj), get_histones_names(chia.obj)))
}


#' Obtain the indices of nodes bearing any or all of a set of transcription factors.
#'
#' @param chia.obj A ChIA-PET object.
#' @param tf.names The names of the transcription factor to look for.
#' @param how.many How many of the selected TFs should be present to return a hit? By default,
#'    all TFs must be present.
#'
#' @return A logical vector indicating which regions bind the given transcription factor.
#' @export
nodes_with_tf <- function(chia.obj, tf.names, how.many=length(tf.names)) {
    nodes.with.generic(chia.obj, tf.names, number_of_tfs, how.many)
}

#' Obtain the indices of nodes bearing any or all of a set of transcription factors.
#'
#' @param chia.obj A ChIA-PET object.
#' @param tf.names The names of the transcription factor to look for.
#' @param how.many How many of the selected TFs should be present to return a hit? By default,
#'    all TFs must be present.
#'
#' @return A logical vector indicating which regions bind the given transcription factor.
#' @export
nodes_with_polymerase <- function(chia.obj, polymerase.names, how.many=length(polymerase.names)) {
    nodes.with.generic(chia.obj, polymerase.names, number_of_polymerases, how.many)
}


#' Obtain the indices of nodes bearing any or all of a set of transcription factors.
#'
#' @param chia.obj A ChIA-PET object.
#' @param tf.names The names of the transcription factor to look for.
#' @param how.many How many of the selected TFs should be present to return a hit? By default,
#'    all TFs must be present.
#'
#' @return A logical vector indicating which regions bind the given transcription factor.
#' @export
nodes_with_histone <- function(chia.obj, histone.names, how.many=length(histone.names)) {
    nodes.with.generic(chia.obj, histone.names, number_of_histones, how.many)
}

#' Obtain the indices of nodes bearing any or all of a set of transcription factors.
#'
#' @param chia.obj A ChIA-PET object.
#' @param tf.names The names of the transcription factor to look for.
#' @param how.many How many of the selected TFs should be present to return a hit? By default,
#'    all TFs must be present.
#'
#' @return A logical vector indicating which regions bind the given transcription factor.
#' @export
nodes_with_chip <- function(chia.obj, chip.names, how.many=length(chip.names)) {
    nodes.with.generic(chia.obj, chip.names, number_of_chips, how.many)
}

nodes.with.generic <- function(chia.obj, object.names, number.function, how.many=length(object.names)) {
    numbers = number.function(chia.obj, object.names)
    return(numbers >= how.many)
}

#' Obtain the number of TFs from a given list that bind all regions within a ChIA-PET object.
#'
#' @param chia.obj A ChIA-PET object.
#' @param tf.names The names of the transcription factor to look for.
#'
#' @return An integer vector containing the number of TF in the list binding each region.
#' @export
number_of_tfs <- function(chia.obj, tf.names) {
    number.of.generic(chia.obj, tf.names, get_tf)
}

#' Obtain the number_of_polymerases from a given list that bind all regions within a ChIA-PET object.
#'
#' @param chia.obj A ChIA-PET object.
#' @param polymerase.names The names of the transcription factor to look for.
#'
#' @return An integer vector containing the number_of_polymerases in the list binding each region.
#' @export
number_of_polymerases <- function(chia.obj, polymerase.names) {
    number.of.generic(chia.obj, polymerase.names, get_polymerases)
}

#' Obtain the number of histone marks from a given list that bind all regions within a ChIA-PET object.
#'
#' @param chia.obj A ChIA-PET object.
#' @param histone.names The names of the histone marks to look for.
#'
#' @return An integer vector containing the number of histone marks in the list binding each region.
#' @export
number_of_histones <- function(chia.obj, histone.names) {
    number.of.generic(chia.obj, histone.names, get_histones)
}

#' Obtain the number of histone marks from a given list that bind all regions within a ChIA-PET object.
#'
#' @param chia.obj A ChIA-PET object.
#' @param histone.names The names of the histone marks to look for.
#'
#' @return An integer vector containing the number of histone marks in the list binding each region.
#' @export
number_of_chips <- function(chia.obj, chip.names) {
    number.of.generic(chia.obj, chip.names, get.chips)
}

number.of.generic <- function(chia.obj, object.names, get.function) {
    tf.data = get.function(chia.obj)
    
    # At the beginning, no TFs have been detected.
    results = rep(0, nrow(tf.data))
    
    for(tf in object.names) {
        if(!any(colnames(tf.data)==tf)) {
            warning(paste0("No information for TF ", tf, " found.\n"))
        } else {
            results = results + ifelse(tf.data[,tf] > 0, 1, 0)
        }
    }
    
    return(results)
}

#' Obtain the proportion of regions intersecting a ChIP dataset, for all ChIP datasets.
#'
#' @param chia.obj A ChIA-PET object.
#'
#' @return A named numeric vector containing the proportion of regions intersecting the ChIP datasets
#' @export
get_all_chip_proportion <- function(chia.obj) {
    ret = c()
    for(chip in get_chips_names(chia.obj)) {
        ret[chip] = sum(number_of_chips(chia.obj, chip))/number_of_nodes(chia.obj)
    }
    
    return(ret)
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