
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# load relevant libraries
library(here)

rm(list = ls() )

# load important functions
source(here("scripts/functions/taxon_query_functions.R"))

# load the equation and data input databases
equ.dat <- readxl::read_xlsx(here("raw_data/equation_data.xlsx"))
equ.dat <- equ.dat[!is.na(equ.dat$equation_id),]

# taxon list
tax.list <- equ.dat[, c("equation_id", "equation_target_taxon", "life_stage")]
tax.list <- split(tax.list, tax.list$equation_id)

# run the function for a few taxa on the list
x.out <- 
  lapply(tax.list[1:10], function(x) {
    
    Get_taxonomic_info(x.name = x[["equation_target_taxon"]], 
                       equ.id = x[["equation_id"]],
                       life_stage = x[["life_stage"]],
                       data.base = "gbif", 
                       rank.diff = 1, tries = 2, ask_or_not = FALSE
                       )
    
  })

### END
