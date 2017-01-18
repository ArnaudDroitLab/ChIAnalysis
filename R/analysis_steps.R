#' Analyze the chromatin state of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to its chromatin states.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned
#'    by \code{\link{annotate_chia}}
#' @param chia.params Original chia.params used to build chia.obj, used to 
#'    determine genomic backgrounds.
#' @param output.dir The name of the directory where to save the graphs.
analyze.chromatin.states <- function(chia.obj, chia.params=NULL, output.dir="output") {
    if(!has_chrom_state(chia.obj)) {
        warning("No chromatin states to analyze!")
    } else {
        facet.col = c("Chrom.State"=6, "Simple.Chrom.State"=2)
        param.name = c("Chrom.State"="input.chrom.state", "Simple.Chrom.State"="simple.chrom.state")
        for(cs in c("Chrom.State", "Simple.Chrom.State")) {
            # Write table of proportions of chromatin states/annotation types
            state.proportions = table(chia.obj$Regions[[cs]])/length(chia.obj$Regions)
            write.table(data.frame(State=names(state.proportions), Proportion=as.vector(state.proportions)),
                        file.path(output.dir, paste0(cs, " summary.txt")), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

            # Plot proportions of chromatin states as a function of connectivity.
            connectivity.df <- categorize_by_connectivity(chia.obj)
            chia_plot_metrics(chia.obj, level_counts, connectivity.df, "Connectivity", "Proportion of nodes in category",
                         graph.type = "line", facet.cols = facet.col[cs],
                         file.out = file.path(output.dir, paste0("Proportion of ", cs, " as a function of connectivity category.pdf")),
                         variable.name = cs)
                         
            # Generate contact_heatmaps.
            contact_heatmap(chia.obj, cs, cs, output.dir=output.dir)
            
            if(has_components(chia.obj)) {
                size.categories <- categorize_by_components_size(chia.obj)
                chia_plot_metrics(chia.obj, level_counts, size.categories, 
                  x.lab = "Size category", y.lab = "Proportion",
                  graph.type = "line", facet.cols = facet.col[cs],
                  file.out = file.path(output.dir, paste0("Proportion of ", cs, " as a function of size category.pdf")),
                  variable.name = cs, proportion = TRUE)
            }
            
            if(!is.null(chia.params)) {
                # Chromatin state enrichment of the whole network.
                region.enrichment(get_granges(chia.obj), 
                                  chia.params[[param.name[cs] ]], 
                                  file.out=file.path(output.dir, paste0(cs, " enrichment of regions in network.pdf")))
            }
        }
    }
}

#' Analyze the annotation of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the available annotation.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned 
#'    by \code{\link{annotate_chia}}
#' @param chia.params Original chia.params used to build chia.obj, used to 
#'    determine genomic backgrounds.
#' @param output.dir The name of the directory where to save the graphs.
analyze.annotation <- function(chia.obj, chia.params=NULL, output.dir="output") {
    if(!has_gene_annotation(chia.obj)) {
        warning("No gene annotation to analyze!")
    } else {
        # Plot genomic region s connectivity.
        connectivity.df <- categorize_by_connectivity(chia.obj)
        chia_plot_metrics(chia.obj, level_counts, connectivity.df, "Connectivity", "Proportion of nodes in category",
                   graph.type = "line", facet.rows = 3, facet.cols = 3,
                   file.out = file.path(output.dir, "Proportion of genomic location as a function of connectivity category.pdf"),
                   variable.name = "Simple.annotation")
        
        # Plot genomic region contact map.
        contact_heatmap(chia.obj, "Simple.annotation", "genomic location", output.dir)
        
        # Plot genomic region vs component size.
        if(has_components(chia.obj)) {
            size.categories <- categorize_by_components_size(chia.obj)
            chia_plot_metrics(chia.obj, level_counts, size.categories,
                 x.lab = "Size category", y.lab = "Proportion", graph.type = "line", facet.rows = 3,
                 file.out = file.path(output.dir, "Proportion of genomic location as a function of size category.pdf"),
                 variable.name = "Simple.annotation", proportion = TRUE)
        }
        
        # Perform network regions enrichment within genomic regions.
        if(!is.null(chia.params)) {
            region.enrichment(get_granges(chia.obj), 
                              chia.params$genomic.regions, 
                              file.out=file.path(output.dir, "Genomic region enrichment of regions in networks.pdf"))         
        }
    }
}

#' Analyze the expression of ChIA-PET data
#'
#' Produce a plot to the expression of a node according to its degree.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The name of the directory where to save the graphs.
analyze.expression <- function(chia.obj, chia.params=NULL, output.dir="output") {
    if(!(has_degree(chia.obj) && has_expression_levels(chia.obj))) {
        warning("No expression levels to analyze!")
    } else {
        # Plot expression as a function of connectivity.
        log.expression = function(x) {
            log2(x$Regions$Expr.mean[x$Regions$Gene.Representative]+1)
        }
        chia_plot_metrics(chia.obj, log.expression, categorize_by_connectivity(chia.obj), graph.type="boxplot", 
                          file.out=file.path(output.dir, "Boxplot of expression as a function of connectivity.pdf"))

        # Plot expression inside network vs outside network.
        if(!is.null(chia.params)) {
            genomewide.expression.vs.network(chia.obj, chia.params, output.dir)
        }
    }
}

#' Analyze the gene specificity of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to gene specificity: \describe{
#' \item{Tau vs degree.pdf}{A plot of the Tau of the nodes according to their degree.}
#' \item{Boxplot of degrees by expression category.pdf}{A boxplot the degrees of the nodes in each expression category.}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The name of the directory where to save the graphs.
analyze.gene.specificity <- function(chia.obj, output.dir="output") {
    if(!has_gene_specificity(chia.obj)) {
        warning("No gene specificity to analyze!")
    } else {
        # Plot Tau and category vs degree
        tissue.specificity.df = with(chia.obj$Regions[chia.obj$Regions$Gene.Representative,],
                                    data.frame(Degree=degree(chia.obj$Graph),
                                               Tau=chia.obj$Regions$Expression.Tau,
                                               Category=chia.obj$Regions$Expression.Category))

        # Plot category vs degree
        ggplot(tissue.specificity.df, aes(x=Category, y=log2(Degree))) +
            geom_boxplot() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
        ggsave(file.path(output.dir, "Boxplot of connectivity by expression category.pdf"))
    }
}

chip.enrichment.vs.background <- function(all.chip, output.dir, suffix) {
    if(length(all.chip) > 0) {
        plot.width = (0.35 * length(all.chip)) + 2.45
        if(!is.null(chia.params$input.chrom.state)) {
            multiple.region.enrichment(all.chip, chia.params$input.chrom.state, plot.width=plot.width, 
                                       file.prefix=file.path(output.dir, paste0("Chromatin states enrichments ", suffix)))
            
            multiple.region.enrichment(all.chip, chia.params$simple.chrom.state, plot.width=plot.width, 
                                       file.prefix=file.path(output.dir, paste0("Simple chromatin states enrichments", suffix)))
        }
        
        multiple.region.enrichment(all.chip, chia.params$genomic.regions, plot.width=plot.width, 
                                   file.prefix=file.path(output.dir, paste0("Genomic regions enrichments", suffix)))        
    }
}

#' Analyze the TF of ChIA-PET data
#' Produce a plot of the presence of TF at the connection points of ChIA-PET data, according to the
#' connectivity of the nodes.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @importFrom reshape2 melt
analyze.tf <- function(chia.obj, chia.params=NULL, output.dir="output") {
    if(has_transcription_factors(chia.obj)) {
        # Plot TF/chip presence against node connectivity.
        connectivity.df <- categorize_by_connectivity(chia.obj)
        metrics = chia_plot_metrics(chia.obj, metric.function=get_all_chip_proportion, 
                                    node.categories=connectivity.df, graph.type="line", 
                                    file.out=file.path(output.dir, "TF presence on contact point by connectivity.pdf"))
        write.table(metrics$Metrics, file=file.path(output.dir, "TF presence on contact point by connectivity.txt"), 
                    sep="\t", col.names=TRUE, row.names=FALSE)
                    
        # Perform TF/chip enrichment against genomic background, for both all TF/ChIP regions
        # and TF/ChIP regions within network.
        if(!is.null(chia.params)) {
            all.chip = c(chia.params$tf.regions, chia.params$pol.regions, chia.params$histone.regions)
            chip.enrichment.vs.background(all.chip, output.dir, "(All)")
            
            in.network.chip = GRangesList()
            for(chip in get_chips_names(chia.obj)) {
                in.network.chip[[chip]] = get_granges(select_by_chip(chia.obj, chip))
            }
            chip.enrichment.vs.background(in.network.chip, output.dir, "(In network)")
        }
    } else {
        warning("No transcription factor to analyze!")
    }
}

#' Analyze the components of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the topology of the nodes in the whole set
#' of regions and in the separated components: \describe{
#' \item{Log2(size of component) vs Log2(Number of components).pdf}{A plot of the number_of_components according to the size of the components.}
#' \item{Component table.txt}{A file with the information about TSS for each component.}
#' \item{Proportion of TSS in low connectivity nodes.pdf}{A plot of the proportion of TSS in nodes with low connectivity (component with 5 nodes or less).}
#' \item{Proportion of TSS in high connectivity nodes.pdf}{A plot of the proportion of TSS in nodes with high connectivity (component with more than 5 nodes).}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
#' @export
analyze_components <- function(chia.obj, output.dir="output") {
  # Analyze components
  if(has_components(chia.obj)) {
    component.table <- chia.obj$Regions[,c("Component.Id", "Component.size")]
    component.df <- data.frame(Size = unique(component.table)$Component.size, Number = 1)
    component.df <- aggregate(Number~Size, data = component.df, FUN = sum)

    ggplot(component.df) + geom_point(mapping=aes(x=log2(Size), y=log2(Number)))
    ggsave(file.path(output.dir, "Log2(size of component) vs Log2(Number of components).pdf"))

    annotate.component <- function(x) {
      result.df = data.frame(NumberOfTSS=sum(x$distanceToTSS==0),
                             Size=nrow(x),
                             data.frame(as.list((table(x$Simple.annotation)/nrow(x))), check.names=FALSE))

      if(!is.null(x$Chrom.State)) {
        result.df = cbind(result.df, data.frame(as.list((table(x$Chrom.State)/nrow(x))), check.names=FALSE))
      }

      return(result.df)
    }

    component.table = ddply(as.data.frame(chia.obj$Regions), "Component.Id", annotate.component)
    write.table(component.table, file=file.path(output.dir, "Component table.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

    component.table = ddply(component.table, "Size", function(x) { return(cbind(x, NumberOfComponentsOfThisSize=nrow(x)))})
    proportion.per.size = ddply(component.table, ~NumberOfTSS * Size, function(x) { return(nrow(x)/x$NumberOfComponentsOfThisSize[1])})
    proportion.per.size[order(proportion.per.size$Size, proportion.per.size$NumberOfTSS),]

    ggplot(subset(proportion.per.size, Size <= 5), aes(x=NumberOfTSS, y=V1)) + geom_bar(stat="identity") + facet_wrap(~Size)
    ggsave(file.path(output.dir, "Proportion of TSS in low connectivity nodes.pdf"))

    ggplot(subset(component.table, Size >5), aes(x=NumberOfTSS/Size)) + geom_histogram()
    ggsave(file.path(output.dir, "Proportion of TSS in high connectivity nodes.pdf"))
    
    if(has_gene_annotation(chia.obj)) {
        tss.metric.function <- functor_constructor(boolean_count, "Is.TSS", proportion=TRUE)
        chia_plot_metrics(chia.obj, apply_single_metric_by_component(tss.metric.function, "TSS proportion"),
                          categorize_by_components_size(chia.obj), graph.type="boxplot",
                          file.out=file.path(output.dir, "Boxplot of the proportion of TSS in fonction of the size of the networks.pdf"))
    }
  }
}

#' Analyze the topology of ChIA-PET data
#'
#' Analyzes ChIA-PET data and produces graphs according to the topology of the nodes in the whole set
#' of regions and in the separated components: \describe{
#' \item{Histogram of number of edges.pdf}{Histogram of the number of edges in the whole set of data.}
#' \item{Scatter plot of degrees of the left node versus the right node.pdf}{Scatter plot showing the relation between the left and right vertices' degree.}
#' \item{Degree vs cluster width.pdf}{A plot of the degree of the nodes in one cluster according to the size of its cluster.}
#' \item{Log2(size of component) vs Log2(Number of components).pdf}{A plot of the number_of_components according to the size of the components.}
#' \item{Component table.txt}{A file with the information about TSS for each component.}
#' \item{Proportion of TSS in low connectivity nodes.pdf}{A plot of the proportion of TSS in nodes with low connectivity (component with 5 nodes or less).}
#' \item{Proportion of TSS in high connectivity nodes.pdf}{A plot of the proportion of TSS in nodes with high connectivity (component with more than 5 nodes).}}
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The name of the directory where to save the graphs.
#'
analyze.generic.topology <- function(chia.obj, output.dir="output") {
  # Plot an histogram of the number of edges.
  #hist(log2(degree(chia.obj$Graph)), breaks=seq(0, 300, by=5))
  ggplot(chia.obj$Regions) + geom_histogram(aes(Degree)) + scale_y_log10() + scale_x_log10()
  ggsave(file.path(output.dir, "Histogram of number of edges.pdf"))

  # Plot a scatter plot showing the relation between the left and right vertices' degree.
  chia.df = get.data.frame(chia.obj$Graph)
  count.per.index = table(c(chia.df$from, chia.df$to))
  left.count = count.per.index[as.character(chia.df$from)]
  right.count = count.per.index[as.character(chia.df$to)]

  ggplot(data.frame(Left=as.vector(left.count), Right=as.vector(right.count)), aes(x=Left, y=Right)) +
    geom_point(alpha=0.01)
  ggsave(file.path(output.dir, "Scatter plot of degrees of the left node versus the right node.pdf"))

  # Degree vs cluster size
  cluster.size.df = data.frame(Degree=degree(chia.obj$Graph), Width=chia.obj$Regions$width)
  ggplot(cluster.size.df, aes(x=Width, y=Degree)) + geom_point() + scale_x_log10() + scale_y_log10()
  ggsave(file.path(output.dir, "Degree vs cluster width.pdf"))
}

#' Performs all possible analysis regarding the central nodes.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The directory where output should be saved.
#'
#' @export
analyze_central_nodes <- function(chia.obj, output.dir=".") {
  if(has_centrality(chia.obj)) {
    centrality.categories <- categorize_by_centrality(chia.obj)
  
    if(has_chrom_state(chia.obj)) {
      chia_plot_metrics(chia.obj, level_counts, centrality.categories, graph.type = "heatmap",
                   file.out = file.path(output.dir, "Heatmap of centrality vs chromatin states.pdf"), variable.name = "Chrom.State")
    }
    
    if(has_gene_annotation(chia.obj)) {
      chia_plot_metrics(chia.obj, level_counts, centrality.categories, graph.type = "heatmap",
                   file.out = file.path(output.dir, "Heatmap of centrality vs genomic locations.pdf"), variable.name = "Simple.annotation")
    }
    
    if(has_transcription_factors(chia.obj)) {
      results = chia_plot_metrics(chia.obj, calculate_tf_presence, centrality.categories, graph.type = "heatmap",
        x.lab = "Centrality category", y.lab = "Proportion", 
        file.out = file.path(output.dir, "Proportion of TF as a function of centrality.pdf"),
        proportion = TRUE)
      
        write.table(results$Metrics, file=file.path(output.dir, "Proportion of TF as a function of centrality.txt"), 
                    sep="\t", row.names=FALSE, col.names=TRUE)
    }
  }
}

#' Performs analyses which do not fit in any of the other sub-analyses.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param output.dir The directory where output should be saved.
#'
#' @export
analyze_misc <- function(chia.obj, output.dir=".") {
    if(has_fitness(chia.obj)) {
        chia_plot_metrics(chia.obj, subset_counts, categorize_by_connectivity(chia.obj), x.lab = "Connectivity",
                          graph.type = "histogram", file.out = file.path(output.dir, "Histogram of the number of essential genes by connectivity category.pdf"),
                          all.conditions = "(Gene.Representative & (Fitness > 0.5))", proportion = FALSE)
    }
}

#' Performs all possible analysis steps on a ChIA-PET object.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate_chia}}.
#' @param chia.params The chia.params object which was used to annotate this one. Used
#'   for calculating in-network TF/chromatin state enrichments vs the entire genome.
#' @param output.dir The directory where output should be saved.
#' @param verbose Set to TRUE to get progress output to the console.
#' @param label The label to add to the heatmap title (name of the variable).
#'
#' @export
analyze_chia_pet <- function(chia.obj, chia.params=NULL, output.dir=".", verbose=TRUE) {
    # If verbose output is turned off, redirect output to a NULL stream.
    cat.sink = ifelse(verbose, "", textConnection(NULL, w))

    # Output the results of the annotation.
    output_annotated_chia(chia.obj, output.dir)

    # Perform further in-depth analysis of the networks.
    cat(date(), " : Analyzing network topologies...\n",cat.sink)
    analyze.generic.topology(chia.obj,  file.path.create(output.dir, "Topology"))

    cat(date(), " : Analyzing network components...\n",cat.sink)
    analyze_components(chia.obj, file.path.create(output.dir, "Components"))
    
    cat(date(), " : Analyzing genomic annotations...\n",cat.sink)
    analyze.annotation(chia.obj, chia.params, file.path.create(output.dir, "Genomic regions"))

    cat(date(), " : Analyzing chromatin states...\n",cat.sink)
    analyze.chromatin.states(chia.obj, chia.params, file.path.create(output.dir, "Chromatin states"))

    cat(date(), " : Analyzing gene expression...\n",cat.sink)
    analyze.expression(chia.obj, chia.params, file.path.create(output.dir, "Expression"))

    cat(date(), " : Analyzing transcription factor overlaps...\n",cat.sink)
    analyze.tf(chia.obj, chia.params, file.path.create(output.dir, "Transcription factors"))

    cat(date(), " : Analyzing gene specificity...\n",cat.sink)
    analyze.gene.specificity(chia.obj, file.path.create(output.dir, "Expression"))

    cat(date(), " : Analyzing central nodes...\n",cat.sink)
    analyze_central_nodes(chia.obj, output.dir)
    
    cat(date(), " : Analyzing miscellaneous...\n",cat.sink)
    analyze_misc(chia.obj, output.dir)
    
    # Close dummy verbose stream.
    if(!verbose) {
        close(cat.sink)
    }
}

file.path.create <- function(...) {
    retval = file.path(...)
    dir.create(retval, recursive=TRUE, showWarnings=FALSE)
    return(retval)
}