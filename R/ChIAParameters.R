#' Create an object storing parameters for ChIA processing.
#'
#' @param input.chrom.state A GRanges object containing chromatin states.
#' @param biosample The biosample identifier from ENCODE. Valid examples are
#'   GM12878, K562 or MCF-7.
#' @param genome.build The name of genomic build all genomic coordinates are 
#'   provided in ("hg38", "hg19", etc).
#' @param tf.regions A \linkS4class{GRangesList} object containing regions 
#'   where transcription factors bind.
#' @param histone.regions A \linkS4class{GRangesList} object containing regions
#'   where histone marks are found.
#' @param pol.regions A \linkS4class{GRangesList} object containing polymerase
#'   binding regions.
#' @param expression.data A data frame containing the levels of expression of 
#'   genes, according to their EMSEMBL id.
#' @param tad.regions A GRanges object containign TAD boundaries.
#' @param compartments.regions A GRanges object containing compartment
#'   boundaries.
#' @param tssRegion A 2-element integer vector indicating how far from the TSS
#'   the "Promoter" region should extend.
#' @param centrality.measures The centrality measures to use. A vector 
#'   containing at least one of "Degree", "Betweenness", "Eigenvector",
#'   "Closeness"
#' @param weight.attr The name of the edge attribute to be used as weight when
#'   calculating centralities or splitting the network into its component 
#'   communities.
#' @param simple.chrom.state A simplified version of the chromatin states. If
#'   NULL, those are inferred from input.chrom.state when possible.
#' @param genomic.regions A partition of the genome into different kind of 
#'   genomic regions, such as "Intron", "Exon", "Promoter", etc.
#' @param gene.annotations A named list providing per-gene annotations. Each 
#'   element of the list should be a named vector, where each element is named
#'   after the gene it annotates.
#'
#' @return An environment that can be passed to \code{\link{process_chia_pet}}, 
#'   \code{\link{annotate_chia.pet}} or \code{\link{analyze_chia_pet}}.
#' @export
build_chia_params <- function(input.chrom.state = NULL, biosample = NULL, genome.build = NULL, tf.regions = NULL,
                             histone.regions = NULL, pol.regions = NULL, expression.data = NULL, tad.regions = NULL,
                             compartments.regions = NULL, tssRegion = c(-3000, 3000), centrality.measures=c("Degree"),
                             weight.attr=NULL, simple.chrom.state=NULL, genomic.regions=NULL, gene.annotations=NULL) {
    chia.params = new.env()
    
    chia.params$biosample = biosample
    chia.params$genome.build = genome.build

    chia.params$input.chrom.state = input.chrom.state
    chia.params$tf.regions = tf.regions
    chia.params$histone.regions = histone.regions
    chia.params$pol.regions = pol.regions
    chia.params$expression.data = expression.data
    chia.params$tad.regions = tad.regions
    chia.params$compartments.regions = compartments.regions
    
    chia.params$tssRegion = tssRegion
    chia.params$centrality.measures = centrality.measures
    chia.params$weight.attr = weight.attr
    chia.params$gene.annotations = gene.annotations
    
               
    # Calculate simplified chromatin states.
    chia.params = simplify.param.chrom.states(chia.params)
    
    # Calculate genomic region partition
    if(is.null(genomic.regions) && !is.null(genome.build)) {
        annot = select.annotations(genome.build)
        genomic.regions = partition.genomic.regions(annot$TxDb, input.chrom.state, annot$BSGenome)
    }
    chia.params$genomic.regions = genomic.regions

    return(chia.params)
}

#' Retrieve whatever data it can from ENCODE and add it to the ChIA parameter object.
#'
#' @param chia.params The ChIA parameters object for which additional annotations must be fetched from ENCODE.
#'
#' @return An environment that can be passed to process_chia_pet, annotate_chia.pet or analyze_chia_pet.
#' @export
add_encode_data <- function(chia.params) {
    biosample = chia.params$biosample
    genome.build = chia.params$genome.build
    
    # If biosample is provided, download missing annotations from ENCODE.
    if(!is.null(biosample) && !is.null(genome.build)) {
        # Download transcription factors
        if (is.null(chia.params$tf.regions)) {
            chia.params$tf.regions <- download.encode.chip(biosample, genome.build)$narrow$Regions
        }

        # Download histone marks
        if (is.null(chia.params$histone.regions)) {
            chia.params$histone.regions <- download.encode.histones(biosample, genome.build)$broad$Regions
        }

        # Download PolII regions.
        if (is.null(chia.params$pol.regions)) {
            chia.params$pol.regions <- download.encode.polymerases(biosample, genome.build)$narrow$Regions
        }

        # Download expression data
        if (is.null(chia.params$expression.data)) {
            chia.params$expression.data <- download.encode.rna(biosample, genome.build)$Expression
            chia.params$expression.data$ENSEMBL = gsub("\\.\\d+$", "", chia.params$expression.data$gene_id)
            chia.params$expression.data$FPKM = log2(chia.params$expression.data$Mean.FPKM + 1)
        }

        # Download chromatin states
        if (is.null(chia.params$input.chrom.state) && genome.build=="hg19") {
            chia.params$input.chrom.state <- download.chromatin.states(biosample, genome.build)
            chia.params = simplify.param.chrom.states(chia.params)
        }
    }
    
    return(chia.params)
}

simplify.param.chrom.states <- function(chia.params) {
    # Make sure we have full chromatin states, and no simplified chromatin states already assigned.
    if(is.null(chia.params$simple.chrom.state) && !is.null(chia.params$input.chrom.state)) {
        # Make sure that chromatin states map the types for our mapping.
        if(all(chia.params$input.chrom.state$name %in% names(simple.chrom.state.map()))) {
            # Calculate simplified states and assign them.
            simple.chrom.state = chia.params$input.chrom.state
            simple.chrom.state$name = simplify.chrom.state(simple.chrom.state$name)
            chia.params$simple.chrom.state = simple.chrom.state
        }
    }
    
    # Return the modified chia.params (Useless since it's an environment, 
    # but useful if we change it back to a list.)
    chia.params
}