#' classHeatmap
#'
#' @name classHeatmap
#' @rdname classHeatmap
#' @author Jason T. Serviss
#' @keywords classHeatmap
#' @param data A data frame with columns \emph{gene} and \emph{class} indicating
#'    the class specific gene expression.
#' @param counts.log A matrix of log normalized counts.
#' @param classes A tibble with columns \emph{sample} and \emph{class}
#'    indicating the sample ID and class (cell type) of the sample.
#' @export
#' @import ggplot2
#' @import tidyverse
#' @importFrom ggthemes theme_few
#' @importFrom viridis scale_fill_viridis

NULL

classHeatmap <- function(data, counts.log, classes) {
  
  #cluster to get gene order
  uClass <- unique(data$class)
  geneOrd <- function(c) {
    IDs <- filter(data, class == c)$gene
    if(length(IDs) == 1) {
      IDs
    } else {
      mat <- counts.log[rownames(counts.log) %in% IDs , ]
      cr <- cor(t(mat), method = "pearson")
      my.dist <- as.dist(1 - cr)
      c <- hclust(my.dist)
      c[[4]][c[[3]]]
    }
  }
  
  order <- map(uClass, geneOrd) %>%
    setNames(uClass) %>%
    namedListToTibble(.) %>%
    rename(class = names, gene = variables) %>%
    mutate(order = 1:n())

  #get counts and normalize
  cx <- counts.log[rownames(counts.log) %in% data$gene, colnames(counts.log) %in% classes$sample] %>%
    apply(., 1, function(x) (x - min(x)) / (max(x) - min(x))) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene") %>%
    gather(sample, value, -gene) %>%
    as_tibble()

  #order samples
  n <- cx %>%
    left_join(classes, by = "sample") %>%
    rename(sampleClass = class) %>%
    arrange(sampleClass) %>%
    mutate(plotSample = factor(sample, levels=unique(sample)))

  #order genes
  m <- n %>%
    left_join(order, by = "gene") %>%
    rename(geneClass = class) %>%
    arrange(order) %>%
    mutate(plotGene = factor(gene, levels=unique(gene)))
  
  #plot
  p <- ggplot(m, aes(plotSample, plotGene)) +
  geom_tile(aes(fill = value)) +
  facet_grid(geneClass ~ sampleClass, scales = "free") +
  scale_fill_viridis() +
  theme_few() +
  theme(
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    strip.text.x = element_text(angle = 90),
    strip.text.y = element_text(angle = 0)
  ) +
  labs(x = "Samples", y = "Genes") +
  guides(
    fill = guide_colourbar(
      title = "z-score",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 10
    )
  )
  
  p
  return(p)
}
