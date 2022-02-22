
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# load relevant libraries
library(here)
library(dplyr)

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
x.samp <- sample(1:length(tax.list), 1)
x.samp <- 34
print(tax.list[x.samp] )

x.out <- 
  lapply(tax.list[x.samp], function(x) {
    
    Get_taxonomic_info(x.name = x[["equation_target_taxon"]], 
                       equ.id = x[["equation_id"]],
                       life_stage = x[["life_stage"]],
                       data.base = "itis", 
                       rank.diff = 2, tries = 5, ask_or_not = FALSE
                       )
    
  })

# check this output
y <- get_taxon_id(database_function = "itis", taxon_name = "Tipulidae", ask_or_not = FALSE, tries = 5)
z <- downstream(sci_id = y[[1]], downto = "species", db = "itis", intermediate = TRUE)
head(z$`118840`)
View(z[[1]])

u <- 
  z$`118840`$intermediate %>%
  bind_rows() %>%
  filter(rankname != "tribe")
View(u)


# try a different approach
u1 <- downstream(sci_id = y[[1]], downto = "subfamily", db = "itis", intermediate = FALSE)
u1$`118840`$tsn[1]



taxize::classification(sci_id = y[[1]], db = "itis")
x.out[[1]]

### END
