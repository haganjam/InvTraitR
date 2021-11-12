
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# load relevant libraries
library(taxize)
library(readr)
library(here)

# decide on a set of taxonomic steps i.e. phylum to sub-phylum, that can be 0.33 or something
# cut everything off at phylum
# which (phylum). Phlym to nrow

# use classification to get upstream taxonomy
# we use downstream to get children

# load the equation and data input databases
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]
in.dat <- readxl::read_xlsx(here("raw_data/variable_input_data.xlsx"))
in.dat <- in.dat[!is.na(in.dat$equation_id),]

# taxon list
tax.list <- equ.dat[, c("equation_id", "equation_target_taxon") ]
tax.list <- split(tax.list, tax.list$equation_id)

x <- tax.list[1:5]
x
x 

# figure out how to use these taxonomic backbones

# taxonomic ranks:
tax.ranks <- c('phylum', 'subphylum','class','subclass', 'superorder','order','suborder','infraorder','superfamily','family', 'subfamily','tribe','subtribe','genus','subgenus','section','subsection', 'species') 
tax.ranks

y <- get_boldid(sci = "Daphnia magna", ask = FALSE)
y[[1]]

x <- classification(sci_id = y[[1]], db = "bold")
downstream(y[[1]], db = 'bold', downto = 'family', intermediate = TRUE)
tax_rank(sci_id = y[[1]], db = "ncbi")

y <- get_tsn(sci = "Daphnia magna", ask = FALSE)
y[[1]]

z <- classification(sci_id = y[[1]], db = "itis")
downstream(y[[1]], db = 'bold', downto = 'family', intermediate = TRUE)

y <- get_gbifid(sci = "Daphnia magna", ask = FALSE)
y[[1]]

u <- classification(sci_id = y[[1]], db = "gbif")

y <- get_eolid(sci_com = "Daphnia magna", ask = TRUE)
y[[1]]

v <- classification(sci_id = y[[1]], db = "eol")


library(dplyr)
full_join(z$`83884`[, -3], x$`26983`[, -3], by = c("rank"))

rank_ref

downstream(sci_id = y[[1]], db = "ncbi")

x.tsn <- get_tsn("Mesostigmata", accepted = TRUE, ask = FALSE)
x.df <- taxize::itis_hierarchy(tsn = x.tsn[[1]], "full")




get_n

gnr_resolve(sci = "Acari")

x.tsn <- get_tsn("Acari")
x.df <- taxize::itis_hierarchy(tsn = x.tsn[[1]], "full")
x.df %>% View()




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





