
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
tax.list <- equ.dat[, c("equation_id", "equation_target_taxon") ]
tax.list <- split(tax.list, tax.list$equation_id)

# try the taxon query function
Get_taxonomic_info(x.name = "Lumbriculidae", 
                   data.base = "itis", 
                   rank.diff = 1, tries = 5, ask_or_not = FALSE,
                   equ.dat = equ.dat,
                   taxon_var = "equation_target_taxon",
                   equation_id = "equation_id"
                   ) 

### END
