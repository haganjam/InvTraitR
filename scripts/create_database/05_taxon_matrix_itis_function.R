
# taxon matric from ITIS

# source: (pg. 12): https://www.itis.gov/pdf/ITIS_ConceptualModelEntityDefinition.pdf

# useful igraph link: https://kevintshoemaker.github.io/NRES-746/graph.theory.html

# load the igraph library
library(igraph)

# build the taxon matrix
tax.mat <- data.frame(from = c("order", "suborder", "infraorder", "section", "subsection", "superfamily", "family", "subfamily", "tribe", "subtribe"),
                      to = c("suborder", "infraorder", "section", "subsection", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus"),
                      weights = c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6, 1/4, 1/4, 1/4, 1/4))

# build igraph from the input matrix
tax.g <- graph_from_data_frame(d = tax.mat, directed=FALSE)

# calculate taxononomic distance between each taxonomic level
tax.d <- distances(
  tax.g,
  v = V(tax.g),
  to = V(tax.g),
  weights = tax.mat$weights,
  mode = c("all"),
  algorithm = c("bellman-ford")
)

### END
