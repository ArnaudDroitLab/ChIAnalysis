#' Associates centrality scores and boolean centrality markers to regions.
#'
#' @param chia.obj ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param which.measures A vector containing the names of the measures to be used
#'   to assess centrality. Those can be "Degree", "Betweenness" and "Eigenvector".
#' @param weight.attr The anme of the edge attribute to be sued as a weight in 
#'   centrality calculations.
#'
#' @export
#' @return The annotated chia.obj.
associate.centralities <- function(chia.obj, which.measures=c("Degree", "Betweenness", "Eigenvector", "Closeness"), weight.attr=NULL) {
  centralities = calculate.centralities(chia.obj, which.measures, weight.attr)
  chia.obj$Regions <- cbind(chia.obj$Regions, centralities)

  return(chia.obj)
}

#' Calculates centralities for all vertices in a graph.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param which.measures A vector containing the names of the measures to be used
#'   to assess centrality. Those can be "Degree", "Betweenness" and "Eigenvector".
#' @param weight.attr The anme of the edge attribute to be sued as a weight in 
#'   centrality calculations.
#'
#' @return A data-frame containing the centrality scores and markers for all vertices in the graph.
#' @export
calculate.centralities <- function(chia.obj, which.measures=c("Degree", "Betweenness", "Eigenvector", "Closeness"), weight.attr=NULL) {
  # If weights should be used, set them as the weight edge attribute.
  if(!is.null(weight.attr)) {
    stopifnot(weight.attr %in% names(edge_attr(chia.obj$Graph)))
    set_edge_attr(chia.obj$Graph, "weight", value = edge_attr(chia.obj$Graph)[[weight.attr]])
  }

  # Define a matrix and a vector to contain network-wide scores and centrality markers.
  results.matrix = matrix(0, nrow=number.of.nodes(chia.obj), ncol=length(which.measures)+1, dimnames=list(NULL, c(which.measures, "Centrality.score")))
  centrality.marker = rep(FALSE, number.of.nodes(chia.obj))

  # Loop over components to measure centrality. Most of these methods can be applied to
  # discontinuous components, but will give different results which might fail
  # to identify a component specific  central node.
  components.out = components(chia.obj$Graph)
  for(i in 1:components.out$no) {
    cat("Processing", i, "out of", components.out$no, "\n")
  
    # Generate a subgraph for the component.
    which.nodes = components.out$membership==i
    component.subgraph = induced_subgraph(chia.obj$Graph, which.nodes)

    # Get the requested centrality measures.
    measures = list()
    if("Degree" %in% which.measures) {
        measures[["Degree"]] = degree(component.subgraph)
    }

    if("Betweenness" %in% which.measures) {
      measures[["Betweenness"]] = betweenness(component.subgraph, directed = FALSE)
    }

    if("Eigenvector" %in% which.measures) {
      measures[["Eigenvector"]] = eigen_centrality(component.subgraph, directed = FALSE)$vector
    }

    if("Closeness" %in% which.measures) {
      measures[["Closeness"]] = estimate_closeness(component.subgraph, cutoff = 3)
    }    
    
    # Make sure we had at least one valid measure.
    stopifnot(length(measures) > 0)

    # Get combined score and determine if nodes should be marked as central.
    measures[["Centrality.score"]] =  apply(as.data.frame(lapply(measures, scale)), 1, mean)

    # Some networks will have constant centrality on all nodes (Two node networks, ring networks, etc.)
    # This will cause standard deviation to be 0 and centrality to be NAN.
    # Deal with these edge cases by assigning a 0 centrality to them.
    measures[["Centrality.score"]][is.nan(measures[["Centrality.score"]])] <- 0
    if(sd(measures[["Centrality.score"]]) != 0) {
        is.central = measures[["Centrality.score"]] > quantile(measures[["Centrality.score"]], probs = 0.95)
    } else {
        is.central = FALSE
    }

    # Report components' results to the combined matrix/vector of all nodes.
    for(measure in names(measures)) {
        results.matrix[which.nodes, measure] = measures[[measure]]
    }
    centrality.marker[which.nodes] = is.central
  }

  # Return a data frame combining the scores and the centrality marker.
  return(data.frame(results.matrix, Is.central=centrality.marker))
}

#' Associates boolean to regions in fonction of their presence in factories
#'
#' Is.In.Factory is \code{TRUE} if the region is in a network with 3 genes or more.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @return The annotated regions.
#' @export
associate.is.in.factory <- function(regions){
  factories <- as.data.frame(regions)
  factories <- aggregate(Is.Gene.Active~Component.Id, data = factories[factories$Gene.Representative,], FUN = sum)
  factories <- factories$Component.Id[factories$Is.Gene.Active > 1]
  regions$Is.In.Factory <- (regions$Component.Id %in% factories)
  return(regions)
}

#' Associates components ids and sizes to chia data, as returned by \code{\link{load.chia}}.
#'
#' @param chia.obj ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param split Should the data be divided into communities?
#' @param oneByOne Sould the netwoks be treated one by one or as a whole?
#' @param method What method sould be used to split data (ignored if split = \code{FALSE}).
#' @return The annotated chia.obj.
#' @export
associate.components <- function(chia.obj){
  # Add data about components
  components.out <- components(chia.obj$Graph)
  chia.obj$Regions$Component.Id <- components.out$membership
  chia.obj$Regions$Component.size <- components.out$csize[components.out$membership]

  # Add a column with the number of edges in each community
  left.df = as.data.frame(chia.left(chia.obj))
  left.df$Edges <- 1
  number.edges.df <- aggregate(Edges~Component.Id, data = left.df, FUN = sum)
  chia.obj$Regions$Component.edges <- number.edges.df$Edges[match(chia.obj$Regions$Component.Id, number.edges.df$Component.Id)]

  return (chia.obj)
}

#' Associate genes to a \linkS4class{GRanges} object from ChIA-PET data.
#'
#' @param regions A \linkS4class{GRanges} object to annotate.
#' @param expression.data A data frame containing the levels of expression of genes, according to their EMSEMBL id.
#'
#' @return The \linkS4class{GRanges} object with associated genes.
#' @importFrom plyr ddply mutate
#' @export
associate.gene <- function(regions, expression.data=NULL) {
  # Associate a gene to a contact only if it's in the promoter.
  promoter.subset = regions$Simple.annotation == "Promoter"

  # Subset the promoter contacts to only keep the highest degrees
  degree.info = ddply(as.data.frame(regions[promoter.subset]), "ENSEMBL", plyr::mutate, max.degree=max(Degree))
  degree.subset = subset(degree.info, Degree == max.degree)

  # Further subset to keep the one closest to the TSS
  distance.info = ddply(degree.subset, "ENSEMBL", plyr::mutate, min.distance=min(abs(distanceToTSS)))
  distance.subset = subset(distance.info, distanceToTSS == min.distance)

  # If there are still more than one row, just pick the first one.
  gene.subset = ddply(distance.subset, "ENSEMBL", function(x) { return(x[1,]) })

  # Add the gene marker to the annotations.
  regions$Gene.Representative = regions$ID %in% gene.subset$ID

  if(!is.null(expression.data)) {
    index.match = match(regions$ENSEMBL, expression.data$ENSEMBL)
    regions$Expr.mean = expression.data$Mean.FPKM[index.match]
  }

  return(regions)
}

#' Associate named regions to a ChIA-PET object.
#'
#' @param The ChIA-PET object to annotate.
#' @param annotation.regions A \linkS4class{GRanges} object with the named 
#'   regions to be used for annotation.
#' @param annotation.name The name to be given to the new column within the 
#'   ChIA-PET object's annotations.
#'
#' @return The ChIA-PET object, with added annotations.
#' @importFrom GenomicRanges findOverlaps
#' @export
associate.regions.by.name <- function(chia.obj, annotation.regions, annotation.name) {
  # If no name are provided, generate some.
  if(is.null(annotation.regions$name)) {
    annotation.regions$name = paste0(seqnames(annotation.regions), ":", start(annotation.regions), "-", end(annotation.regions))
  }
  
  # Find the first overlapping region.
  overlap.indices = findOverlaps(get.granges(chia.obj), annotation.regions, select="first")
  
  # Associate the region names to the chia object.
  chia.obj$Regions[[annotation.name]] = annotation.regions$name[overlap.indices]
  
  return(chia.obj)
}

#' Adds gene-specific annotations to chia.obj.
#'
#' Given a data frame where the first column are gene symbols, and the second column
#' are values associated with the given gene, this function finds teh gene representatives
#' nodes for all symbols and attaches the given annotation to it.
#'
#' @param chia.obj A chia object to annotate.
#' @param annotation.df A data frame where the first column is a gene symbol, and the
#'   second one is a value to be associated to that gene.
#' @param label The label the annotation should be given in the chia object.
#'
#' @return The ChIA object with the added annotation.
#' @export
add.gene.annotation <- function(chia.obj, annotation.df, label) {
    # Associate the values with their gene.
    values = annotation.df[match(chia.obj$Regions$SYMBOL, annotation.df[,1]), 2]
    
    # Only keep the values for gene representatives.
    values[!chia.obj$Regions$Gene.Representative] <- NA
    
    # Add the values to the chia.obj's regions.
    chia.obj$Regions[,label] = values
    
    return(chia.obj)
}

#' Annotate "chia.obj", given as parameter.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{load.chia}}.
#' @param chia.param An object returned by build.chia.params indicating which data should be used for annotating the ChIA object.
#' @param output.dir The name of the directory where to write the selected annotations.
#' @param verbose If TRUE, console output is suppressed.
#' @param skip.centrality If TRUE, centrality scores are not computed.
#'
#' @return The annotated "\code{chia.obj}".
#'
#' @export
annotate.chia <- function(chia.obj, chia.param, output.dir=".", verbose=TRUE, skip.centrality=FALSE) {
  # Create the output directory.
  dir.create(output.dir, recursive = TRUE, showWarnings=FALSE)

  # If verbose output is turned off, redirect output to a NULL stream.
  cat.sink = ifelse(verbose, "", textConnection(NULL, w))

  # Add an ID to every region.
  chia.obj$Regions$ID = 1:nrow(chia.obj$Regions)

  # Add degree count directly to chia.obj$Regions
  chia.obj$Regions$Degree = degree(chia.obj$Graph)

  # Associate genomic regions
  cat(date(), " : Associating genomic regions...\n",cat.sink)
  tmp = associate.genomic.region(get.granges(chia.obj),
                                 chia.param$genome.build,
                                 output.dir = output.dir,
                                 tssRegion = chia.param$tssRegion)
  chia.obj$Regions <- as.data.frame(tmp)
  chia.obj$Regions$Is.TSS <- chia.obj$Regions$distanceToTSS == 0

  # Associate chromatin states
  if(!is.null(chia.param$input.chrom.state)) {
    cat(date(), " : Associating chromatin states...\n",cat.sink)
    chia.obj$Regions$Chrom.State = chrom.state.match(get.granges(chia.obj), chia.param$input.chrom.state)
    chia.obj$Simple.Chrom.State = simplify.chrom.states(chia.obj$Regions$Chrom.State)
  }

  # Associate transcription factors
  if(!is.null(chia.param$tf.regions)) {
    cat(date(), " : Associating transcription factors...\n",cat.sink)
    tmp = associate.tf(get.granges(chia.obj), chia.param$tf.regions)
    chia.obj$Regions = as.data.frame(tmp)
  }

  # Associate histone marks
  if(!is.null(chia.param$histone.regions)) {
    cat(date(), " : Associating histone marks...\n",cat.sink)
    hist.matrix = region.overlap(get.granges(chia.obj), chia.param$histone.regions)
    colnames(hist.matrix) = paste0("HIST.", colnames(hist.matrix))
    chia.obj$Regions = cbind(chia.obj$Regions, hist.matrix)
  }

  # Associate known polymerase binding
  if(!is.null(chia.param$pol.regions)) {
    cat(date(), " : Associating polymerase II regions...\n",cat.sink)
    pol.matrix = region.overlap(get.granges(chia.obj), chia.param$pol.regions)
    colnames(pol.matrix) = paste0("POL.", colnames(pol.matrix))
    chia.obj$Regions = cbind(chia.obj$Regions, pol.matrix)
  }

  # Associate genes to regions
  cat(date(), " : Associating genes...\n",cat.sink)
  tmp = associate.gene(get.granges(chia.obj), chia.param$expression.data)
  chia.obj$Regions = as.data.frame(tmp)

  if(!is.null(chia.params$tad.regions)) {
    chia.obj = associate.regions.by.name(chia.obj, chia.params$tad.regions, "TAD")
  }
  
  if(!is.null(chia.params$compartments.regions)) {
    chia.obj = associate.regions.by.name(chia.obj, chia.params$compartments.regions, "Compartment")
  }
  
  # Associate tissue specificity and fitness scores if the organism allows it.
  if(chia.param$genome.build=="hg19" || chia.param$genome.build=="hg38") {
    cat(date(), " : Associating tissue specificity...\n",cat.sink)
    tmp = associate.tissue.specificity.human(get.granges(chia.obj))
    chia.obj$Regions = as.data.frame(tmp)
    cat(date(), " : Associating fitness score...\n",cat.sink)
    tmp = associate.fitness.genes(get.granges(chia.obj))
    chia.obj$Regions = as.data.frame(tmp)
  }

  # Associate components ids and sizes
  cat(date(), " : Associating components...\n",cat.sink)
  chia.obj = associate.components(chia.obj)

  # Associate centrality scores
  if(!skip.centrality) {
    cat(date(), " : Associating centrality scores...\n",cat.sink)
    chia.obj = associate.centralities(chia.obj,
                                      which.measures=chia.param$centrality.measures,
                                      weight.attr=chia.param$weight.attr)
  }

  # Associate miscellaneous useful defnitions
  chia.obj$Regions = associate.is.gene.active(chia.obj$Regions)
  chia.obj$Regions = associate.is.in.factory(chia.obj$Regions)

  # If verbose output was turned off, close the dummy stream.
  if(!verbose) {
    close(cat.sink)
  }
  return(chia.obj)
}
