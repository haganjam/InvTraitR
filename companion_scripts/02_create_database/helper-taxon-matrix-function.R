#' @title Taxonomic hierarchy framework
#' @description The different taxonomic levels and the distances (weights)
#'  between them for the GBIF taxonomic backbone
#' @details This method relies on having a taxonomic framework which describes
#'  how different species are taxonomically related to each other. The GBIF
#'  taxonomic backbone uses the seven major Linnaen ranks (i.e. kingdom,
#'  phylum, class, order, family, genus and species). It also includes
#'  sub-species and varieties
#'  (https://data-blog.gbif.org/post/gbif-backbone-taxonomy/).
#'  Details of the backbone can be found at:
#'  https://data-blog.gbif.org/post/gbif-backbone-taxonomy/.
#'  Then, we assign numerical distances between taxonomic levels where one unit
#'  reflects on whole major Linnaen rank. Since we only consider ranks from
#'  order to species, the distances were chosen as: 1, 1 and 1/4
#'  (for genus to species). We create the distance matrix using igraph. This
#'  link was useful for setting this up:
#'  https://kevintshoemaker.github.io/NRES-746/graph.theory.html
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @return igraph object representing taxonomic distances among levels

# load the igraph library
library(igraph)

# build the taxon matrix
tax_mat <- data.frame(
  from = c("order", "family"),
  to = c("family", "genus"),
  weights = c(1, 1)
)

# build igraph from the input matrix
tax_g <- graph_from_data_frame(d = tax_mat, directed = FALSE)

# calculate taxononomic distance between each taxonomic level
tax_d <- distances(
  tax_g,
  v = V(tax_g),
  to = V(tax_g),
  weights = tax_mat$weights,
  mode = c("all"),
  algorithm = c("bellman-ford")
)
