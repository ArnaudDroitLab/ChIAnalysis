#' Create a heatmap by connectivity group.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param variable.name The name of the variable according to which the heatmap should be computed.
#' @param label The name to give to the variable name in the resulting heatmap.
#' @param output.dir The name of the directory where to save the heatmaps.
#' @param log.scale Should the color scale be logged?
#' @param proportions Should proportions of contacts be represented rather than the absolute numbers?
#' @param stretch.scale Should the scale extend from 0 to 100?
#' @param show.labels Should each tile contain a label with its number/proportion?
#' @param file.name The name of the file where the heatmap should be saved.
#'
#' @return The matrix used to create the heatmap.
#'
#' @importFrom reshape2 melt
#' @export
contact.heatmap <- function(chia.obj, variable.name, label=NULL, output.dir=NULL, 
                            log.scale=FALSE, proportions=FALSE, stretch.scale=proportions, show.labels=FALSE, 
                            file.name=file.path(output.dir, paste0("Contact heatmap for ", label, ".pdf"))) {
  type.df = data.frame(Left=chia.left(chia.obj)[,variable.name],
                       Right=chia.right(chia.obj)[,variable.name],
                       stringsAsFactors=FALSE)

  var.levels = levels(chia.left(chia.obj)[,variable.name])
  results.matrix = matrix(NA, nrow=length(var.levels), ncol=length(var.levels), dimnames=list(var.levels, var.levels))

  for(i in 1:(length(var.levels))) {
    for(j in i:length(var.levels)) {
      count = sum(type.df$Left==var.levels[i] & type.df$Right==var.levels[j] |
                    type.df$Left==var.levels[j] & type.df$Right==var.levels[i],
                  na.rm=TRUE)
      results.matrix[i, j] = count
    }
  }
  
  if(proportions) {
    results.matrix = results.matrix / nrow(type.df)
  }
  
  scale.name = ifelse(proportions, "Percentage of contacts", "No of contacts")
  if(log.scale) {
    results.matrix = log10(results.matrix+1)
    scale.name = paste0("log10(", scale.name, ")")
  }
  
  results.df = melt(results.matrix, varnames=c("Var1", "Var2"))
  results.df = results.df[!is.na(results.df$value),]

  results.df$Var1 = factor(results.df$Var1, levels = rev(var.levels))
  results.df$Var2 = factor(results.df$Var2, levels = var.levels)
  if(proportions) {
    plot.obj = ggplot(results.df, aes(y=Var1, x=Var2, fill=value*100))
  } else {
    plot.obj = ggplot(results.df, aes(y=Var1, x=Var2, fill=value))
  }
  
  if(stretch.scale) {
    limits=c(0,100)
  } else {
    limits=c(min(results.df$Value), max(results.df$Value))
  }      
  
  plot.obj = plot.obj +  
               geom_tile(color="black") +
               labs(x=NULL, y=NULL) +
               scale_fill_gradient(low="white", high="red", name=scale.name, limits=limits) +
               theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.key.size = unit(1, "cm"))

  if(show.labels) {
    # Can't sue the same albel size for 4 rows as we do for 18.
    # 2 is a good size for 18 rows, 8 is a good size for 4. Extrapolate
    # a linear relation from that, and you get:
    label.size=(-3/7)*length(unique(results.df$Var1))+(68/7)
    if(proportions) {
      plot.obj = plot.obj + geom_text(mapping=aes(y=Var1, x=Var2, label=sprintf("%.0f%%", value*100)), size=label.size)
    } else if(log.scale) {
      plot.obj = plot.obj + geom_text(mapping=aes(y=Var1, x=Var2, label=round(10^value)), size=label.size)
    } else {
      plot.obj = plot.obj + geom_text(mapping=aes(y=Var1, x=Var2, label=value), size=label.size)
    }
  }
  ggsave(file.name)
  write.table(results.matrix, file=paste0(file.name, ".txt"), sep="\t", row.names=TRUE, col.names=TRUE)

  return(results.matrix)
}

#' Produces boxplots for every transcription factor, in relation to its 3D connectivity
#'
#' For every transcription factor in the ChIP data, creates a boxplot of the force of the ChIP-seq signal
#' in function of the contact frequence of the region.
#'
#' @param chip.data A \linkS4class{GRangesList} containing the regions of the ChIp-seq data, with signal values.
#' @param biosample The biosample identifier from ENCODE. Valid examples are GM12878, K562 or MCF-7.
#' @param genome.build The name of the chosen annotation ("hg38", "hg19").
#' @param chia.obj Annotated ChIA-PET data, as returned by \link{analyze.chia.pet} or \link{annotate.chia}.
#' @param output.dir The directory where to write the boxplots.
#' @param TSS Should only the TSS regions be kept?
#' @param tssRegion tssRegion A vector with the region range to TSS.
#'
#' @importFrom cowplot plot_grid
#' @export
boxplot.per.tf <- function(chip.data, biosample, genome.build, chia.obj, output.dir, TSS = TRUE, tssRegion = c(-3000, 3000)) {

  # Extract ChIA-PET regions
  chia.data <- get.granges(chia.obj)

  # Exctract all TF
  if (genome.build == "hg19"){
    TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if (genome.build == "hg38"){
    TXDB <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if (genome.build == "mm9"){
    TxDb <- TxDxDb.Mmusculus.UCSC.mm9.knownGene::TxDxDb.Mmusculus.UCSC.mm9.knownGene
  } else if (genome.build == "mm10"){
    TxDb <- TxDxDb.Mmusculus.UCSC.mm10.knownGene::TxDxDb.Mmusculus.UCSC.mm10.knownGene
  } else {
    stop("Error : genome.build is invalid.")
  }

  tss.regions <- genes(TxDb)
  tss.regions <- promoters(tss.regions, abs(tssRegion[1]), tssRegion[2])


  # Function to create boxplot with histogram
  create.boxplot <- function(chip.data, chia.data, label, label.x, label.y, output.dir, tss.regions, TSS=TRUE){
    if (TSS){
      chip.data <- chip.data[chip.data$distanceToTSS == 0]
      chia.data <- chia.data[chia.data$distanceToTSS == 0]
    }
    indices <- GenomicRanges::findOverlaps(chip.data, chia.data)
    if (length(indices) != 0) {

      signal.degree.df <- data.frame(Signal = log2(chip.data$signalValue[indices@from]),
                                     Degree = chia.data$Degree[indices@to])
      signal.degree.df$CutDegree <- cut(signal.degree.df$Degree, breaks = c(1, 5, 10, 20, 40, Inf), right = FALSE)
      if (length(chip.data$signalValue[-indices@from]) != 0){
        signal.degree.df <- rbind(data.frame(Signal = log2(chip.data$signalValue[-indices@from]),
                                             Degree = 0, CutDegree = "0"), signal.degree.df)
      }
      box <- ggplot(signal.degree.df) + geom_boxplot(aes(CutDegree, Signal)) + ylab(label.y) + xlab(label.x) + ggtitle(label)


      if (TSS) {
        tss.indices <- findOverlaps(tss.regions, chia.data)
        tss.degree.df <- data.frame(Region = tss.indices@from, Degree = chia.data$Degree[tss.indices@to])
        tss.degree.df$CutDegree <- cut(tss.degree.df$Degree, breaks = c(1, 5, 10, 20, 40, Inf), right = FALSE)
        if (length(chip.data$signalValue[-indices@from]) != 0){
          tss.degree.df <- rbind(data.frame(Region = c(1:length(tss.regions))[-tss.indices@from],
                                            Degree = 0, CutDegree = "0"), tss.degree.df)
        }

        mapping.df <- data.frame(x.pos = levels(tss.degree.df$CutDegree), y.pos = as.vector(table(signal.degree.df$CutDegree)),
                                 label = paste0(round(as.vector(table(signal.degree.df$CutDegree) / table(tss.degree.df$CutDegree))*100, digits = 1), "%"))
        hist <- ggplot() +
          geom_bar(aes(CutDegree), data = tss.degree.df, fill = "light blue") +
          geom_bar(aes(CutDegree), data = signal.degree.df) +
          geom_text(data = mapping.df, aes(x = x.pos, y = y.pos, label = label), vjust = -1) +
          xlab(label.x) +
          ggtitle(paste("Histogram of ", label.x))


      } else {
        hist <- ggplot(signal.degree.df) + geom_bar(aes(CutDegree)) + xlab(label.x) + ggtitle(paste("Histogram of ", label.x))
      }

      plot_grid(box, hist, nrow = 2, align = "v")
      ggsave(file.path(output.dir, paste0(label, ".pdf")), height = 14, width = 7)
    }
  }

  # Create plot for every TF
  dir.create(file.path(output.dir, biosample), recursive = TRUE, showWarnings=FALSE)
  for (tf in names(chip.data)){
    cat("Factor : ", tf, "\n")
    chip.subset <- chip.data[tf][[1]]
    chip.subset <- annotate.chip(chip.subset, input.chrom.state = NULL, tf.regions = NULL, histone.regions = NULL,
                                 pol.regions = NULL, expression.levels = NULL, genome.build = genome.build,
                                 biosample = biosample, tssRegion = tssRegion, output.dir = "output/annotations",
                                 label = paste0(tf, " annotated"))

    if (TSS){
      create.boxplot(chip.subset, chia.data,
                     paste0("Boxplot of log2(Signal) in fct of Degree of ", tf ," at TSS"),
                     "Contact frequency at TSS", "log2(Signal)", file.path(output.dir, biosample), tss.regions, TSS = TSS)
    } else {
      create.boxplot(chip.subset, chia.data,
                     paste0("Boxplot of log2(Signal) in fct of Degree of ", tf),
                     "Contact frequency", "log2(Signal)", file.path(output.dir, biosample), tss.regions, TSS = TSS)
    }

  }
}

#' Histogram of the proportion of "essential genes", by connectivity categories.
#'
#' @param chia.obj A list containing the ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param essential.genes The ID's of the genes to anlayze. The ID's refer the the ChIA-PET ID's of the nodes.
#' @param label.x The title of the histogram.
#' @param output.dir The name of the directory where to save the plot.
#' @export
histogram.essential.genes <- function(chia.obj, essential.genes, label.x, output.dir){
  # Convert to data.frame
  chia.annotated.df <- chia.obj$Regions
  # Cut the size into categories
  chia.annotated.df$CutDegree <- cut(chia.annotated.df$Degree, breaks = c(1, 2, 6, 21, Inf), right = FALSE,
                                     labels = c("Singles", "Low", "Intermediate", "High"))
  # Keep only the ids that represent "essential genes"
  chia.essential <- chia.annotated.df[chia.annotated.df$ID %in% essential.genes,]
  # Prepare the mapping to add the percentage label over the bars of the histogram
  mapping.df <- data.frame(x.pos = levels(chia.essential$CutDegree), y.pos = as.vector(table(chia.essential$CutDegree)),
                           label = paste0(round(as.vector(table(chia.essential$CutDegree) / table(chia.annotated.df$CutDegree[chia.annotated.df$Gene.Representative]))*100, digits = 1), "%"))
  # Create the histogram
  ggplot(chia.essential) +
    geom_bar(aes(CutDegree)) +
    geom_text(data = mapping.df, aes(x = x.pos, y = y.pos, label = label), vjust = -1) +
    xlab("connectivity") +
    ggtitle(paste("Histogram of ", label.x))
  ggsave(file.path(output.dir, paste0("Histogram of ", label.x, ".pdf")))
}

#' Plot heatmaps in fonction of a ChIA-PET annotation
#'
#' Creates heatmaps of the proportion of a variable from the ChIA-PET annotation, given as parameter, in fonction of the networks.
#' To use the transcription factors as variable, write "TF". This will create the same heatmaps, but with presence or absence of each factor instead of proportions.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param size.limit The networks must have a size over size.limit.
#' @param variable.name Name of the column to use, or "TF" for transcription factors.
#' @param label The label to add to the heatmap title (name of the variable).
#' @param output.dir The name of the directory where output should be saved.
#'
#' @importFrom plyr ddply
#' @importFrom reshape2 dcast
#' @importFrom NMF aheatmap
#'
#' @export
chia.plot.network.heatmap <- function(chia.obj, size.limit, variable.name, label, output.dir) {
  dir.create(output.dir, recursive = TRUE, showWarnings=FALSE)

  presence.by.tf <- function(chia.df){
    tf.col <- grep("TF.overlap", colnames(chia.df))
    tf.presence <- apply(chia.df[,tf.col], 2, sum)
    tf.presence <- ifelse(tf.presence > 0, 1, 0)
    names(tf.presence) <- sub("TF.overlap.", "", colnames(chia.df)[tf.col])
    return(tf.presence)
  }

  if (variable.name == "TF") {
    per.component.cs = ddply(as.data.frame(chia.obj$Regions[chia.obj$Regions$Component.size>size.limit]), ~Component.Id,
                             presence.by.tf)

  } else {
    per.component.cs = ddply(as.data.frame(chia.obj$Regions[chia.obj$Regions$Component.size>size.limit]), ~Component.Id,
                             function(x) { as.data.frame(table(x[, variable.name])/nrow(x)) })
    # Un-melt it.
    per.component.cs = dcast(per.component.cs, Component.Id~Var1)
  }

  component.ids = per.component.cs$Component.Id

  per.component.cs.matrix = per.component.cs[,-1]
  component.sizes = chia.obj$Regions$Component.size[match(component.ids, chia.obj$Regions$Component.Id)]
  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix), annCol=list("log2(Network size)"=log2(component.sizes)), color="-Blues")
  dev.off()


  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes no clustering.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix[order(component.sizes),]), annCol=list("log2(Network size)"=log2(component.sizes)[order(component.sizes)]), Colv=NA, color="-Blues")
  dev.off()

  per.component.cs.matrix[per.component.cs.matrix > 0.4] = 0.3
  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes limit 0.3.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix), annCol=list("log2(Network size)"=log2(component.sizes)), color="-Blues")
  dev.off()


  pdf(file.path(output.dir, paste0("Heatmap of the ", label, " - networks with more than ", size.limit, " nodes no clustering limit 0.3.pdf")), width=7, height=7)
  aheatmap(t(per.component.cs.matrix[order(component.sizes),]), annCol=list("log2(Network size)"=log2(component.sizes)[order(component.sizes)]), Colv=NA, color="-Blues")
  dev.off()
}

# Compare expression of genes inside and outside of networks.
genomewide.expression.vs.network <- function(chia.obj, chia.params, output.dir) {
    annot = select.annotations(chia.params$genome.build)
    
    # Get all genes which could have an expression value.
    all.possible.genes = AnnotationDbi::keys(annot$TxDb, "GENEID")
    
    # Get the ENSEMBL ids associated with the Entrez ID
    ensembl.ids = mapIds(annot$OrgDb, as.character(all.possible.genes), c("ENSEMBL"), "ENTREZID")
    symbols = mapIds(annot$OrgDb, as.character(all.possible.genes), c("SYMBOL"), "ENTREZID")
    expression.data = data.frame(Entrez=names(ensembl.ids), ENSEMBL=ensembl.ids, SYMBOL=symbols)
    
    # Associate expression data
    expression.data$FPKM = chia.params$expression.data$FPKM[match(expression.data$ENSEMBL, chia.params$expression.data$ENSEMBL)]
    
    # Keep only those with EnsemblIDs and at least SOME expression
    expression.data = expression.data[!is.na(expression.data$ENSEMBL) & !is.na(expression.data$FPKM) & expression.data$FPKM > 0,]
    
    in.chia = expression.data$ENSEMBL %in% chia.obj$Regions$ENSEMBL[chia.obj$Regions$Gene.Representative]
    in.chia.low.connect = expression.data$ENSEMBL %in% chia.obj$Regions$ENSEMBL[chia.obj$Regions$Gene.Representative & chia.obj$Regions$Degree < 20]
    in.chia.high.connect = expression.data$ENSEMBL %in% chia.obj$Regions$ENSEMBL[chia.obj$Regions$Gene.Representative & chia.obj$Regions$Degree >= 20]
    
    expression.data$In.ChIA = in.chia
    expression.data$In.ChIA = in.chia.low.connect
    expression.data$In.ChIA = in.chia.high.connect
    
    plot.df = rbind(data.frame(Category="All", FPKM=expression.data$FPKM),
                    data.frame(Category="Outside of networks", FPKM=expression.data$FPKM[!in.chia]),
                    data.frame(Category="In networks (all)", FPKM=expression.data$FPKM[in.chia]),
                    data.frame(Category="In networks (low connectivity)", FPKM=expression.data$FPKM[in.chia.low.connect]),
                    data.frame(Category="In networks (high connectivity)", FPKM=expression.data$FPKM[in.chia.high.connect]))
    plot.df$Category = factor(plot.df$Category, levels=c("All", "Outside of networks", "In networks (all)", "In networks (low connectivity)", "In networks (high connectivity)"))
                    
    plot.subset = plot.df[plot.df$Category %in% c("Outside of networks", "In networks (all)"),]
    ggplot(plot.subset, aes(x=Category, y=log2(FPKM))) + geom_boxplot()
    ggsave(file.path(output.dir, "Expression inside networks vs outside networks.pdf"), width=7, height=7)
    
    plot.subset = plot.subset[plot.subset$FPKM >= 1,]
    ggplot(plot.subset, aes(x=Category, y=log2(FPKM))) + geom_boxplot()
    ggsave(file.path(output.dir, "Expression of active genes inside networks vs outside networks.pdf"), width=7, height=7)
    
    return(plot.df)
}

#' Represents contact-frequencies between different categories of elements
#' as a stylized graph.
#'
#' Creates a graph where each node represents a node category from the original chia.obj,
#' and edge width represents the frequency of those contacts. Node size is
#' also proportional to the number of nodes in the given category.
#'
#' @param chia.obj A list containing the annotated ChIA-PET data, as returned by \code{\link{annotate.chia}}.
#' @param col.name The name of the column on which the contact graph should be based.
#' @param file.name The name of the output file.
#'
#' @import igraph
#'
#' @export
fancy.topology <- function(chia.obj, col.name, file.name=paste0("Abstract inter-category contact graph for ", col.name, ".pdf")) {
    # Remove any node with unknown chromatin state.
    chia.subset = chia.vertex.subset(chia.obj, !is.na(chia.obj$Regions[[col.name]]))
    
    # Get list of "left" and "right" items.
    left = chia.left(chia.subset)[[col.name]]
    right = chia.right(chia.subset)[[col.name]]

    # Cross-tabulate, and group elements where only order differ.
    cross.table = table(data.frame(Left=as.integer(left), Right=as.integer(right)))
    for(i in 1:(ncol(cross.table) - 1)) {
        for(j in (i+1):nrow(cross.table)) {
            cross.table[j, i] = cross.table[j, i] + cross.table[i, j]
            cross.table[i, j] = 0
        }
    }
    
    # Remove zero elements (So only X -> Y remain, and all Y -> X are lost)
    cross.df = melt(cross.table)
    cross.df = cross.df[cross.df$value != 0,]

    # Make a graph
    graph.obj = make_graph(c(rbind(cross.df$Left, cross.df$Right)), directed = FALSE)    
    edge_attr(graph.obj) <- data.frame(Weight=cross.df$value)
    vertex_names = levels(chia.obj$Regions[[col.name]])
    vertex_weights = as.integer(table(chia.obj$Regions[[col.name]])[vertex_names])
    vertex_attr(graph.obj) <- data.frame(Name=vertex_names, Weight=vertex_weights)

    pdf(file.name, width=21, height=14)
    plot(graph.obj, 
         vertex.label = as.character(vertex_attr(graph.obj)$Name),
         vertex.size = vertex_attr(graph.obj)$Weight / 1000,
         edge.width = edge_attr(graph.obj)$Weight / 1000,
         layout=layout_in_circle)
    dev.off()
    
    # Generate tables to accompany graph.
    colnames(cross.table) = levels(left)
    rownames(cross.table) = levels(left)
    
    sink(paste0(file.name, ".txt"))
    print(list(EdgeWeights=cross.table,
               VertexWeights=table(chia.obj$Regions[[col.name]])))
    sink(NULL)
}
