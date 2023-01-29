#' @title Special taxon names
#' @description A vector of special taxon names that are not searched using
#'  the regular search framework
#' @details Some taxa are often not properly identified in freshwater field
#'  research. For example, different worm taxa, tardigrads, rotifers etc. are
#'  often specified as simply: Oligochaeta, Rotifera, Nematoda etc. Thus,
#'  these special names can not be used in the more general taxonomic
#'  similarity framework and are rather matched to a specific set of equations
#'  developed for those taxa. We provide several of these for different
#'  geographical locations in the database so that they can be matched as
#'  accurately as possible.
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' @return vector of special taxon names
special_taxon_names <- function() {
  c(
    "Rotifera",
    "Tardigrada",
    "Nematoda",
    "Platyhelminthes",
    "Turbellaria",
    "Annelida",
    "Oligochaeta"
  )
}
