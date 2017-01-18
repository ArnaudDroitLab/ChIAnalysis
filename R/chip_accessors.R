
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
get_chips <- function(chia.obj, overlap.threshold=0) {
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
    nodes_with_generic(chia.obj, tf.names, number_of_tfs, how.many)
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
    nodes_with_generic(chia.obj, polymerase.names, number_of_polymerases, how.many)
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
    nodes_with_generic(chia.obj, histone.names, number_of_histones, how.many)
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
    nodes_with_generic(chia.obj, chip.names, number_of_chips, how.many)
}

nodes_with_generic <- function(chia.obj, object.names, number.function, how.many=length(object.names)) {
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
    number_of_generic(chia.obj, tf.names, get_tf)
}

#' Obtain the number_of_polymerases from a given list that bind all regions within a ChIA-PET object.
#'
#' @param chia.obj A ChIA-PET object.
#' @param polymerase.names The names of the transcription factor to look for.
#'
#' @return An integer vector containing the number_of_polymerases in the list binding each region.
#' @export
number_of_polymerases <- function(chia.obj, polymerase.names) {
    number_of_generic(chia.obj, polymerase.names, get_polymerases)
}

#' Obtain the number of histone marks from a given list that bind all regions within a ChIA-PET object.
#'
#' @param chia.obj A ChIA-PET object.
#' @param histone.names The names of the histone marks to look for.
#'
#' @return An integer vector containing the number of histone marks in the list binding each region.
#' @export
number_of_histones <- function(chia.obj, histone.names) {
    number_of_generic(chia.obj, histone.names, get_histones)
}

#' Obtain the number of histone marks from a given list that bind all regions within a ChIA-PET object.
#'
#' @param chia.obj A ChIA-PET object.
#' @param histone.names The names of the histone marks to look for.
#'
#' @return An integer vector containing the number of histone marks in the list binding each region.
#' @export
number_of_chips <- function(chia.obj, chip.names) {
    number_of_generic(chia.obj, chip.names, get_chips)
}

number_of_generic <- function(chia.obj, object.names, get.function) {
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
