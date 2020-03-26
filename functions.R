# Better rounding function than R's base round
myround <- function(x) { trunc(x + 0.5) }


# Scales reads by 
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding 
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "floor") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- myround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

# Modified igraph functions

plot_network_custom <- function (g, physeq = NULL, type = "samples", color = NULL, shape = NULL, 
                                 point_size = 4, alpha = 1, label = "value", hjust = 1.35, 
                                 line_weight = 0.5, line_color = color, line_alpha = 0.4, 
                                 layout.method = layout.fruchterman.reingold, title = NULL, label_size = 3) 
{
  if (vcount(g) < 2) {
    stop("The graph you provided, `g`, has too few vertices. \\n         Check your graph, or the output of `make_network` and try again.")
  }
  if (type %in% c("taxa", "species", "OTUs", "otus", "otu")) {
    type <- "taxa"
  }
  edgeDF <- data.frame(get.edgelist(g))
  edgeDF$id <- 1:length(edgeDF[, 1])
  vertDF <- layout.method(g)
  colnames(vertDF) <- c("x", "y")
  vertDF <- data.frame(value = get.vertex.attribute(g, "name"), 
                       vertDF)
  if (!is.null(physeq)) {
    extraData <- NULL
    if (type == "samples" & !is.null(sample_data(physeq, 
                                                 FALSE))) {
      extraData = data.frame(sample_data(physeq))[as.character(vertDF$value), 
                                                  , drop = FALSE]
    }
    else if (type == "taxa" & !is.null(tax_table(physeq, 
                                                 FALSE))) {
      extraData = data.frame(tax_table(physeq))[as.character(vertDF$value), 
                                                , drop = FALSE]
    }
    if (!is.null(extraData)) {
      vertDF <- data.frame(vertDF, extraData)
    }
  }
  graphDF <- merge(reshape2::melt(edgeDF, id = "id"), vertDF, 
                   by = "value")
  p <- ggplot(vertDF, aes(x, y))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), axis.text.x = element_blank(), 
                              axis.text.y = element_blank(), axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), axis.ticks = element_blank(), 
                              panel.border = element_blank())
  p <- p + geom_point(aes_string(color = color, shape = shape), 
                      size = point_size, na.rm = TRUE)
  if (!is.null(label)) {
    p <- p + geom_text(aes_string(label = label), size = label_size, 
                       hjust = hjust, na.rm = TRUE)
  }
  p <- p + geom_line(aes_string(group = "id", color = line_color), 
                     graphDF, size = line_weight, alpha = line_alpha, na.rm = TRUE)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}
