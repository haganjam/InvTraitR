
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# load relevant libraries
library(taxize)
library(readr)
library(here)

# next to do

# does it only make sense to go down rather than up??


# load the equation and data input databases
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]
in.dat <- readxl::read_xlsx(here("raw_data/variable_input_data.xlsx"))
in.dat <- in.dat[!is.na(in.dat$equation_id),]

# taxon list
tax.list <- equ.dat[, c("equation_id", "equation_target_taxon") ]
tax.list <- split(tax.list, tax.list$equation_id)

x <- tax.list[1:5]

tax.list <- 
  
  lapply(x, function(y) {
  
  x.name <- y$equation_target_taxon
  
  x.test <- get_tsn(x.name, accepted = TRUE)
  
  if (attr(x.test, "match") == "not found") {
    
    x.df <- NULL
    
  } else {
    
    x.tsn <- get_tsn(x.name, accepted = TRUE)
    
    x.df <- taxize::itis_hierarchy(tsn = x.tsn[[1]], "full")
    
  }
  
  if (nrow(x.df) == 1) {
    
    x.df <- NULL
    
  } else {
    
    x.df$equation_id <- y$equation_id
    
    x.df <- x.df[, c("equation_id", "rankname", "taxonname", "tsn")]
    
    # add a column to identify the focal species and the taxonomic rank
    x.df$focal_species = if_else(x.df$taxonname == x.name, 1, 0)
    x.df$taxon_rank <- rep(1:length(unique(x.df$rankname)), rle(x.df$rankname)$lengths)
    
  }
  
  return(x.df)
  
})

tax.list[!sapply(tax.list, is.null)]

x 

y <- get_uid(sci_com = "Mesostigmata", ask = FALSE)
y[[1]]

classification(sci_id = y[[1]], db = "ncbi")
downstream(y[[1]], db = 'ncbi', downto = 'infraorder', intermediate = TRUE)
tax_rank(sci_id = y[[1]], db = "ncbi")

rank_ref

downstream(sci_id = y[[1]], db = "ncbi")

x.tsn <- get_tsn("Mesostigmata", accepted = TRUE, ask = FALSE)
x.df <- taxize::itis_hierarchy(tsn = x.tsn[[1]], "full")




get_n

gnr_resolve(sci = "Acari")

x.tsn <- get_tsn("Acari")
x.df <- taxize::itis_hierarchy(tsn = x.tsn[[1]], "full")
x.df %>% View()



