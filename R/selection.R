#' Select all nodes which are gene representatives.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A subset of chia.obj containing only nodes which are gene representatives.
#' @export
select_gene_representative <- function(chia.obj) {
    return(chia_vertex_subset(chia.obj, is_gene_representative(chia.obj)))
}                   

#' Select all nodes which are central.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A subset of chia.obj containing only nodes which are central.
#' @export                   
select_central_nodes <- function(chia.obj) {
    return(chia_vertex_subset(chia.obj, chia.obj$Regions$Is.central))
}

#' Select all nodes which are in a factory.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A subset of chia.obj containing only nodes which are in a factory.
#' @export                   
select_factories <- function(chia.obj) {
    return(chia_vertex_subset(chia.obj, chia.obj$Regions$Is.In.Factory))
}

#' Generate a function which selects all nodes whose ID is present in the passed-in chia.obj.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A function which accepts a ChIA object and returns the subset of nodes whose
#'   ID are in the ChIA object passed to this function.
#' @export   
select_from_chia_functor <- function(chia.obj) {
    force(chia.obj)
    return(function(x) {
        return(chia_vertex_subset(x, x$Regions$ID %in% chia.obj$Regions$ID))
    })
}

#' Select all nodes which are in the given components.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A subset of chia.obj containing only nodes which are in the given components.
#' @export  
select_by_components <- function(chia.obj, component.ids) {
    return(chia_vertex_subset(chia.obj, chia.obj$Regions$Component.Id %in% component.ids))
}

#' Generate a function similar to select_by_components where the components.ids parameter is bound.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#'
#' @return A function which accepts a ChIA object and returns the subset of nodes whose
#'   component is in component.ids.
#' @export   
select_by_component_functor <- function(component.ids) {
    force(component.ids)
    function(chia.obj) {
        return(select_by_components(chia.obj, component.ids))
    }
}

#' Select all nodes in the components where the given genes are found.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load_chia}}.
#' @param genenames A vector fo gene names whose components should be selected.
#'
#' @return A subset fo chia.obj containing only the components of the given genes.
#' @export   
select_gene_component <- function(chia.obj, genenames) {
    component.ids = get_gene_annotation(chia.obj, genenames)$Component.Id
    return(select_by_components(chia.obj, component.ids))
}

#' Returns a subset of a ChIA object containing only components classified as gene networks.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#' @param min.gene The minimum number of gene in a component for it to be considered
#'   a gene network.
#'
#' @return A subsetted chia object containing only the gene networks.
#' @export
select_gene_networks <- function(chia.obj, min.gene=2) {
    return(chia_vertex_subset(chia.obj, regions_in_gene_network(chia.obj, min.gene)))
}

#' Returns a subset of a ChIA object containing only nodes with the specified 
#' transcription factors.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#' @param tf.names The transcription factors whose regions must be returned.
#'
#' @return A subset of the ChIA object containing only nodes with the specified 
#'   transcription factors.
#' @export
select_by_tf <- function(chia.obj, tf.names) {
    select_by_chip.generic(chia.obj, tf.names, nodes_with_tf)
}

#' Returns a subset of a ChIA object containing only nodes with the specified 
#' transcription factors.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#' @param polymerase.names The polymerases whose regions must be returned.
#'
#' @return A subset of the ChIA object containing only nodes with the specified 
#'   polymerases.
#' @export
select_by_polymerase <- function(chia.obj, polymerase.names) {
    select_by_chip.generic(chia.obj, polymerase.names, nodes_with_polymerases)
}

#' Returns a subset of a ChIA object containing only nodes with the specified 
#' transcription factors.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#' @param histone.names The histone marks whose regions must be returned.
#'
#' @return A subset of the ChIA object containing only nodes with the specified 
#'   histone marks.
#' @export
select_by_histone <- function(chia.obj, histone.names) {
    select_by_chip.generic(chia.obj, histone.names, nodes_with_histones)
}

#' Returns a subset of a ChIA object containing only nodes with the specified 
#' transcription factors.
#'
#' @param chia.obj The ChIA object whose regions must be assessed.
#' @param chip.names The ChIP whose regions must be returned.
#'
#' @return A subset of the ChIA object containing only nodes with the specified 
#'   ChIPs.
#' @export
select_by_chip <- function(chia.obj, chip.names) {
    select_by_chip.generic(chia.obj, chip.names, nodes_with_chip)
}

select_by_chip.generic <- function(chia.obj, chip.names, nodes.with.function) {
    chia_vertex_subset(chia.obj, nodes.with.function(chia.obj, chip.names))
}
