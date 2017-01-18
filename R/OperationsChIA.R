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

regions.to.vertex.attr <- function(chia.obj) {
  vertex_attr(chia.obj$Graph) <- as.data.frame(chia.obj$Regions, stringsAsFactors = FALSE)
  return(chia.obj)
}

vertex.attr.to.regions <- function(graph.obj) {
  region.df = as.data.frame(vertex_attr(graph.obj), stringsAsFactors = FALSE)
  region.df$strand = '*'
  return(GRanges(region.df))
}

#' Read and load a ChIA-PET output file.
#'
#' @param input.chia The path of the file containing the ChIA-PET data.
#' @param excluded.chr Chromosomes to be excluded from the analysis, such as
#'    mitochondrial chromosomes.
#'
#' @return A list of 4 elements: \describe{
#' \item{$Left}{A \linkS4class{GRanges} object containing the information about the "left side" of the ChIA-PET data.}
#' \item{$Right}{A \linkS4class{GRanges} object containing the information about the "right side" of the ChIA-PET data.}
#' \item{$Regions}{A \linkS4class{GRanges} object containing the reduced information of both sides.}
#' \item{$Graph}{A directed \linkS4class{igraph} object picturing every interaction between the left side and the right
#' side of the ChIA-PET data.}}
#'
#' @importFrom plyr ddply
#' @importFrom plyr summarize
#'
#' @export
load_chia <- function(input.chia, excluded.chr=c()) {
    chia.raw = read.table(input.chia, sep="\t")
    chia.raw = chia.raw[,1:7]
    colnames(chia.raw) = c("L.chr", "L.start", "L.end", "R.chr", "R.start", "R.end", "Reads")
    
    chia.raw = chia.raw[!(chia.raw$L.chr %in% excluded.chr) & !(chia.raw$R.chr %in% excluded.chr),]
    
    # Separate and extend on both sides
    split.raw.chia <- function(chia.raw, columns, flank.size=0) {
       result = chia.raw[,columns]
       colnames(result) = c("chr", "start", "end")
       result$start = pmax(0, result$start - flank.size)
       result$end = result$end + flank.size

       return(GRanges(result))
    }

    chia_left.ranges = split.raw.chia(chia.raw, 1:3)
    chia_right.ranges = split.raw.chia(chia.raw, 4:6)

    # Reduce to a single set of coordinates
    single.set = reduce(unlist(GRangesList(chia_left.ranges, chia_right.ranges)))

    # Build a graph.
    # Map back to original contact points
    chia_left.indices = findOverlaps(chia_left.ranges, single.set, select="first")
    chia_right.indices = findOverlaps(chia_right.ranges, single.set, select="first")

    # Find and remove self loops.
    mapped.df = cbind(chia.raw, Left=chia_left.indices, Right=chia_right.indices)
    mapped.df = mapped.df[mapped.df$Left != mapped.df$Right,]

    # Summarize multiple edges.
    max.df = ddply(mapped.df, ~Left*Right, summarize, L.chr=head(L.chr, n=1), L.start=min(L.start), L.end=max(L.end),
                                                      R.chr=head(R.chr, n=1), R.start=min(R.start), R.end=max(R.end), Reads=sum(Reads))

    # Remap IDs: igraph will create as many nodes as max(ids), which will cause problems for us.
    only.unique = sort(unique(c(max.df$Left, max.df$Right)))
    remapped.ids = data.frame(Original=only.unique, Remapped=1:length(only.unique))
    new.left = remapped.ids$Remapped[match(max.df$Left, remapped.ids$Original)]
    new.right = remapped.ids$Remapped[match(max.df$Right, remapped.ids$Original)]
                                                      
    # Create iGraph object and set the original coordinates and the number of supporting reads as edge attributes.
    chia.graph = make_graph(c(rbind(new.left, new.right)), directed=FALSE)
    edge_attr(chia.graph) <- max.df

    # Regions which were only part of a self-loop will have been filtered above
    # and will not be part of the graph object. This will cause a difference between
    # the lengths of single.set and chia.graph.
    # To fix this, we remove all regions which are not in max.df.
    single.set = single.set[only.unique]
    
    chia.obj = list(Regions=as.data.frame(single.set), Graph=chia.graph)
    class(chia.obj) = "ChIA"
    
    return(chia.obj)
}

#' Save annotated ChIA-PET data.
#'
#' Writes files: \itemize{
#'   \item The first one separates the "left" and "right" sides of the ChIA-PET data with annotations.
#'   \item The second set of files is made of components files. They group all interactions in a single component, with
#'        an extra column containing the number of reads supporting the data. Their format is supported by Cytoscape.}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The name of the directory where to save the files.
#'
#' @export
output_annotated_chia <- function(chia.obj, output.dir="output") {
  # Create output directory if it does not exist.
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

  # Export region annotation, putting the ID in the first column.
  chia.data <- as.data.frame(chia.obj$Regions)
  chia.data <- cbind(chia.data$ID, chia.data[,-which(colnames(chia.data) == "ID")])
  write.table(chia.data, file = file.path(output.dir, "Annotated CHIA-PET regions.txt"), row.names = FALSE, sep = "\t", quote=FALSE)

  # Write out annotated interactions by concatening left and right annotations.
  left.df = as.data.frame(chia_left(chia.obj))
  colnames(left.df) <- paste("Left", colnames(left.df), sep=".")
  right.df = as.data.frame(chia_right(chia.obj))
  colnames(right.df) <- paste("Right", colnames(right.df), sep=".")

  write.table(data.frame(left.df, right.df), file=file.path(output.dir, "Interactions.txt"), sep="\t", col.names=TRUE, quote=FALSE)

  # Output Cytoscape components  
  # Create interactions table
  if(is.null(edge_attr(chia.obj$Graph)$Reads)) {
    read.counts = NA
  } else {
    read.counts = edge_attr(chia.obj$Graph)$Reads
  }
  ids <- data.frame(left.df$Left.ID, right.df$Right.ID, read.counts, left.df$Left.Component.Id)
  colnames(ids) <- c("Source", "Target", "Reads", "Component")

  # Export networks in csv files
  cyto.dir = file.path(output.dir, "Cytoscape networks")
  dir.create(file.path(cyto.dir, "Size between 3 and 5 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(cyto.dir, "Size between 6 and 20 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(cyto.dir, "Size between 21 and 50 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(cyto.dir, "Size between 51 and 100 nodes (incl)"), recursive = TRUE, showWarnings=FALSE)
  dir.create(file.path(cyto.dir,"Size over 100 nodes"), recursive = TRUE, showWarnings=FALSE)

  for (i in unique(ids$Component)){
    network <- ids[ids$Component == i,]
    size <- length(unique(c(network$Source, network$Target)))
    if (size > 2){
      if (size < 6){
        dir <- "Size between 3 and 5 nodes (incl)"
      } else if (size < 21){
        dir <- "Size between 6 and 20 nodes (incl)"
      } else if (size < 51){
        dir <- "Size between 21 and 50 nodes (incl)"
      } else if (size < 101){
        dir <- "Size between 51 and 100 nodes (incl)"
      } else {
        dir <- "Size over 100 nodes"
      }
      write.table(network[1:3],
                  file = file.path(cyto.dir, dir, paste0("Create network for components ", i, "(", size, " nodes)", ".csv")),
                  sep = ",", row.names = FALSE, quote=FALSE)
    }
  }
  chia.data <- as.data.frame(chia.obj$Regions)
  chia.data <- cbind(chia.data$ID, chia.data[,-which(colnames(chia.data) == "ID")])
  colnames(chia.data)[1] <- "ID"
  write.table(chia.data, file = file.path(output.dir, "Annotated CHIA-PET regions.txt"), row.names = FALSE, sep = "\t", quote=FALSE)
}


#' Identifies edges crossing community borders.
#'
#' @param input.graph The igraph whose community-crossing edges must be identified,
#' @param method A method returning a \code{community} object, such as igraph::cluster_fast_greedy.
#' @param weight.attr The name fo the edge attribute to be used as edge weight.
#'
#' @return The ids of the edges to be removed.
#' @export
get_crossing_edges <- function(input.graph, method = igraph::cluster_fast_greedy, weight.attr=NULL){
  weights = NULL
  if(!is.null(weight.attr)) {
    if(weight.attr %in% names(edge_attr(input.graph))) {
      weights = edge_attr(input.graph)[[weight.attr]]
    } else {
      warning("The provided weight attribute does not exist.")
    }
  }
  
  communities <- method(as.undirected(input.graph), weights = weights)
  to.delete = crossing(communities, input.graph)
  return(edge_attr(input.graph)$original.id[to.delete])
}

#' Associates components ids and sizes to chia data, as returned by \code{\link{load_chia}}.
#'
#' @param chia.obj ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param split Should the data be divided into communities?
#' @param oneByOne Sould the netwoks be treated one by one or as a whole?
#' @param method What method sould be used to split data (ignored if split = \code{FALSE}).
#' @return The annotated chia.obj.
#' @export
community_split <- function(chia.obj, oneByOne = FALSE, method = igraph::cluster_fast_greedy, weight.attr=NULL) {
  edge_attr(chia.obj$Graph)$original.id = 1:ecount(chia.obj$Graph)
  if (oneByOne){
    # Keep a record of edges marked for deletion, so that we can delete them all
    # at once. This prevents issues with edges being relabeled.
    marked.for.deletion = c()

    # Loop over components one by one.
    components.out = components(chia.obj$Graph)
    for (i in 1:components.out$no){
      # Get the component subgraph.
      component.subgraph = induced_subgraph(chia.obj$Graph, components.out$membership==i)

      # Split it into communities and record teh deleted edges.
      crossing.edges = get_crossing_edges(component.subgraph, method = method, weight.attr=weight.attr)
      marked.for.deletion = c(marked.for.deletion, crossing.edges)
    }

    # Delete removed edges in the original chia object.
    chia.obj$Graph = delete_edges(chia.obj$Graph, marked.for.deletion)
  } else {
    # Split it into communities and delete the necessary edges.
    crossing.edges = get_crossing_edges(chia.obj$Graph, method = method, weight.attr=weight.attr)
    chia.obj$Graph = delete_edges(chia.obj$Graph, crossing.edges)
  }

  # Remove the original.id edge attribute, since it won't be needed anymore.
  edge_attr(chia.obj$Graph)$original.id = NULL

  # Update the degree attribute of regions if it is present.
  chia.obj$Regions$Degree = degree(chia.obj$Graph)
  
  # Also update components, if present.
  if(has_components(chia.obj)) {
    chia.obj = associate_components(chia.obj)
  }
  return(chia.obj)
}

#' Subset a CHIA object.
#'
#' @param chia.obj The chia object to be subset.
#' @param indices The indices of the vertices to be kept.
#'
#' @return A chia object containing only the selected vertices.
#' @export
chia_vertex_subset <- function(chia.obj, indices) {
    if(mode(indices) == "logical") {
      indices = which(indices)
    }
    
    return(list(Regions = chia.obj$Regions[indices,],
                Graph = induced_subgraph(chia.obj$Graph, indices)))
}

#' Subset a chia object by keeping all components where at least one node is selected.
#'
#' @param chia.obj The chia object to be subset.
#' @param indices The indices of the vertices to be kept.
#'
#' @return A chia object containing only the selected vertices.
#' @export
chia_component_subset <- function(chia.obj, indices, min.selection=1) {
    stopifnot(has_components(chia.obj))
    
    selected.nodes.per.component = table(chia.obj$Regions$Component.Id[indices])
    selected.components = names(selected.nodes.per.component)[selected.nodes.per.component >= min.selection]
    return(chia_vertex_subset(chia.obj, chia.obj$Regions$Component.Id %in% selected.components))
}

#' Applies a function to a set chia object subsets.
#'
#' For each column of node.categories, generates a subset of chia.obj and applies chia.function.
#' The results of all calls are returned in a named list.
#'
#' @param chia.obj The chia object on which the functions must be applied.
#' @param chia.function The function to be applied to all given subsets of the chia object.
#' @param node.categories A data-frame with a number of rows equal to the number of regions
#'   in chia.obj, where each column represents a set of indices indicating which nodes belong
#'   in a givenc ategory.
#'
#' @return A list containing th results of the calls to chia.function.
#' @export
category_apply <- function(chia.obj, chia.function, node.categories, ...) {
  # Calculate metrics for all categories
  result.list = list()
  for(node.category in names(node.categories)) {
    graph.subset = chia_vertex_subset(chia.obj, node.categories[[node.category]])
    result.list[[node.category]] = chia.function(graph.subset, ...)
  }
  
  return(result.list)
}

#' Apply a metric function to each component in the ChIA object.
#'
#' @param chia.obj The chia object whose compoennts should have their metrics evaluated.
#' @param metric.function The metric function to apply to all components.
#' @return A list containing the measures metrics for all components.
#' @export
apply_by_component <- function(chia.obj, metric.function, ...) {
    return(category_apply(chia.obj, metric.function, categorize_by_component(chia.obj), ...))
}

#' Apply a metric function returning a single metric to all components of a ChIA object,
#' and return teh resulting values as a single multi-value metric.
#'
#' @param metric.function The metric function to apply to all components.
#' @param metric.name The name for the resulting combined metric.
#' @return A list containing the measures metrics for all components.
#' @export
apply_single_metric_by_component <- function(metric.function, metric.name) {
    function(chia.obj) {
        results = list(unlist(apply_by_component(chia.obj, function(x) { metric.function(x) })))
        names(results) = metric.name
        return(results)
    }
}


#' Analyze ChIA-PET data and produce graphs.
#'
#' @param input.chia The file containing processed ChIA-PET data.
#' @param output.dir The name of the directory where output should be saved.
#' @return The annotated chia.obj.
#' @export
process_chia_pet <- function(input.chia, chia.param, output.dir="output", verbose=TRUE) {
    # Create output directory.
    dir.create(file.path(output.dir), recursive=TRUE, showWarnings=FALSE)

    # Load interaction data.
    chia.obj = load_chia(input.chia)

    # Annotate the ChIA object.
    chia.obj <- annotate_chia(chia.obj, chia.param, output.dir=output.dir, verbose=verbose)

    # Analyze the ChIA network.
    analyze_chia_pet(chia.obj, output.dir)

    # Return teh created object.
	return(chia.obj)
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


#' Contracts an interaction network, keeping only the vertices in indices.
#'
#' Suppose we have a ring network, 1->2->3->4->5->1.
#' Now suppose the passed indices are 2, 4 and 5.
#' Thus, nodes 1 and 3 would be removed, and new edges from 2 to 4 and 5 to
#' 2 would be added. The new network would be:
#' 2->4->5->2.
#'
#' @param chia.obj The ChIA object to be contracted.
#' @param chip.names The indices of the regions to keep in the output chia object.
#'
#' @return A contracted ChIA object.
#' @import igraph
#' @export
chia_contract <- function(chia.obj, indices) {
    # Get the vertex IDs we're going to keep.
    kept.ids = paste0("ChIA-", as.character(chia.obj$Regions$ID[indices]))
    graph.obj = chia.obj$Graph
   
    # Attribute names to nodes to help working with igraph methods.
    vertex_attr(graph.obj) <- list(name=paste0("ChIA-", as.character(chia.obj$Regions$ID)))
   
    one.substitution = TRUE
    while(one.substitution) {
      one.substitution = FALSE
      # Shortcut for vertez names
      vertex.names = vertex_attr(graph.obj)$name
      
      # Initially, all nodes are part of their own set.
      vertex.mapping = vertex.names
      names(vertex.mapping) = vertex.names
      
      # Loop over IDs, marking all neighbors which are not in the input set.
      for(i in kept.ids) {
          # Identify neighbors
          neighborhood = neighbors(graph.obj, i)   
         
          # Remove any identifiers which is in the input set.
          neighbor.ids  = setdiff(names(neighborhood), kept.ids)
         
          if(length(neighbor.ids)>0) {
              vertex.mapping[neighbor.ids] = i
              one.substitution = TRUE
          }
      }
      
      # Convert the vertex mapping to new numeric IDs.
      integer.mapping = 1:length(unique(vertex.mapping))
      names(integer.mapping) = unique(vertex.mapping)
      
      temp.graph = contract(graph.obj, integer.mapping[vertex.mapping])
      new.names = unlist(lapply(vertex_attr(temp.graph)$name, function(x) { ifelse(length(x)>1, intersect(x, kept.ids), x) }))
      vertex_attr(temp.graph) = list(name=new.names)
      
      graph.obj = temp.graph
    }
    
    # Match vertices back to their original regions.
    new.regions.indices = match(gsub("ChIA-", "", vertex_attr(graph.obj)$name), chia.obj$Regions$ID)
    vertex_attr(graph.obj) <- list()
    
    new.chia = list(Graph=simplify(graph.obj),
                    Regions=chia.obj$Regions[new.regions.indices,])

    return(new.chia)
}