
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
x.samp <- sample(1:length(tax.list), 5)
x.out <- 
  lapply(tax.list[x.samp], function(x) {
    
    Get_taxonomic_info(x.name = x[["equation_target_taxon"]], 
                       equ.id = x[["equation_id"]],
                       life_stage = x[["life_stage"]],
                       data.base = "itis", 
                       rank.diff = 2, tries = 5, ask_or_not = FALSE
                       )
    
  })

x.out[[1]]

# process the output into a queriable list

rank.difference <- 1
foc_tax <- "Saltiseiidae"
life_stage <- NA

equ_id <- NA
z <- NA
i <- 1

n <- sample(x = 1:length(x.out), length(x.out), replace = FALSE)
x <- x.out[n]

while ( is.na(equ_id) & ( class(z) != "numeric" ) ) {
  
  if (i == length(x.out)){
    break
  }
  
  # extract information from list
  eq <- x[[i]][["equation_id"]]
  sy <- x[[i]][["synonymns"]]
  ls <- x[[i]][["life_stage"]]
  
  ti <- x[[i]][["taxonomic_information"]]
  
  foc_tax_search <- foc_tax
  
  # if it is a synonymn then replace it with the actual focal name
  if ( foc_tax_search %in% sy ) {
    
    foc_tax_search <- ti[ti$focal_taxa == 1,]$name
    
  }
  
  # search the taxonomic database
  if (foc_tax_search %in% ti[["name"]]) {
    
    z <- abs( ti[ti[["focal_taxa"]] == 1, ]$rank_number - ti[ti[["name"]] == foc_tax_search, ]$rank_number )
    pm <- ifelse(ti[ti[["focal_taxa"]] == 1, ]$name == foc_tax_search, TRUE, FALSE)
    equ_id <- ifelse(z <= rank.difference, eq, NA)
    
  } else {
    
    z <- NA
    pm <- NA
    equ_id <- NA
    
  }
  
  # check if it is the correct life-stage
  if ( !is.na(life_stage) & (life_stage == ls)  ) {
    
    equ_id <- equ_id
    
  } else if ( is.na(life_stage) & is.na(ls) ) {
    
    equ_id <- equ_id
    
  } else {
    
    equ_id <- NA
    
  }
  
  i <- i + 1
  
}

### END
