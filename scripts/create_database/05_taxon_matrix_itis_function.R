#'
#' @title ITIS taxonomic hierarchy framework
#'
#' @description The different taxonomic levels and the distances (weights) between them
#'
#' @details This method relies on having a taxonomic framework which describes how different
#' species are taxonomically related to each other. For the ITIS database, this consists of 
#' 11 different levels: "order", "suborder", "infraorder", "section", "subsection", "superfamily", 
#' "family", "subfamily", "tribe", "subtribe". For more information on the ITIS taxonomic framework
#' see: source: (pg. 12): https://www.itis.gov/pdf/ITIS_ConceptualModelEntityDefinition.pdf. In addition,
#' we need to decide what the distances are between these different taxonomic levels. These can
#' be set to anything but we decided on: 1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/4, 1/4, 1/4, 1/4. This
#' reflects that one taxonomic level i.e. family to genus represents a distance of one. All steps
#' in between are then simply proportional. We create the distance matrix using igraph. This link was
#' useful for setting this up: https://kevintshoemaker.github.io/NRES-746/graph.theory.html
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @return igraph object representing taxonomic distances among levels 
#'

# load the igraph library
library(igraph)

# build the taxon matrix
tax.mat <- data.frame(from = c("order", "suborder", "infraorder", "section", "subsection", "superfamily", "family", "subfamily", "tribe", "subtribe"),
                      to = c("suborder", "infraorder", "section", "subsection", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus"),
                      weights = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/4, 1/4, 1/4, 1/4))

# build igraph from the input matrix
tax.g <- graph_from_data_frame(d = tax.mat, directed=FALSE)

# calculate taxononomic distance between each taxonomic level
tax.itis <- distances(
  tax.g,
  v = V(tax.g),
  to = V(tax.g),
  weights = tax.mat$weights,
  mode = c("all"),
  algorithm = c("bellman-ford")
)

### END
