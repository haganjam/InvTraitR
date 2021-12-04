
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# Function to search the taxonomic database and pull out the correct equation ID and taxonomic distance

# args

# focal_taxa_name - name of the taxon that the database will search for
# life_stage - name of the life stage for the taxon name of interest
# rank.difference - acceptable diffference in rank between focal taxon and equation in the database
# length_only - whether the function should only search for equations based on length
# taxon_database - the taxon database that should be searched
# var_input_data - database associated with the taxon_database that has details about the equations


focal_taxa_name = "Sinantherina"
life_stage = NA
rank.difference = 1
length_only = FALSE
taxon_database <- x.out
var_input_data <- equ.dat

get_mass_from_length( focal_taxa_name = "Hexarthra libica",
                      life_stage = NA,
                      rank.difference = 1,
                      length_in = 2.3 )

length_in = 2.3

taxon_db <- taxon_database

  # if only length equations are requested, then we subset them out first
  if (length_only == TRUE) {
    
    v <- lapply(split(var_input_data, var_input_data$equation_id), function(y) {
      
      ifelse(!any(y$size_measurement != "body_length"), unique(y$equation_id), NA)
      
    })
    v <- unlist(v, use.names = FALSE)
    
    tx.data.ids <- unlist(lapply(taxon_db, function(x) { x[["equation_id"]] } ), use.names = FALSE)
    tx.data.ids[tx.data.ids %in% v[!is.na(v)]]
    
    taxon_db <- taxon_db[tx.data.ids %in% v[!is.na(v)]]
    
  }

u <- 
  lapply(taxon_db, function(x) { 
  
  z <- x[["taxonomic_information"]]$name
  w <- x[["synonymns"]]
  
  (focal_taxa_name %in% z) | (focal_taxa_name %in% w)
  
  })

taxon_db <- taxon_db[unlist(u, use.names = FALSE)]

suitable_equations <- 
  
  lapply(taxon_db, function(x) {

  # extract information from list
  eq <- x[["equation_id"]]
  sy <- x[["synonymns"]]
  lis <- x[["life_stage"]]
  
  ti <- x[["taxonomic_information"]]
  
  foc_tax_search <- focal_taxa_name
  
  # if it is a synonymn then replace it with the actual focal name
  if ( foc_tax_search %in% sy ) {
    
    foc_tax_search <- ti[ti$focal_taxa == 1,]$name
    
  }
  
  # search the taxonomic database
  if (foc_tax_search %in% ti[["name"]]) {
    
    rank_distance <- abs( ti[ti[["focal_taxa"]] == 1, ]$rank_number - ti[ti[["name"]] == foc_tax_search, ]$rank_number )
    pm <- ifelse(ti[ti[["focal_taxa"]] == 1, ]$name == foc_tax_search, TRUE, FALSE)
    equ_id <- ifelse(near_equal( x = rank_distance , y = rank.difference , tol = 1.5e-8 , mode = "ne.lt"), 
                     eq, 
                     NA)
    
  } else {
    
    rank_distance <- NA
    pm <- NA
    equ_id <- NA
    
  }
  
  # check if it is the correct life-stage
  if ( !is.na(life_stage) & (life_stage == lis)  ) {
    
    equ_id <- equ_id
    
  } else if ( is.na(life_stage) & is.na(lis) ) {
    
    equ_id <- equ_id
    
  } else {
    
    equ_id <- NA
    
  }
  
  # package the output
  output <- 
    list("equation_id" = equ_id,
         "perfect_match" = pm,
         "taxonomic_difference" = rank_distance)
  
  return(output)
  
} )

suitable_equations

tax_difference <- 
  
  lapply(suitable_equations, function(z) {
  
  z[["taxonomic_difference"]]
  
}) %>%
  unlist(use.names = FALSE)

if ( tax_difference == tax_difference ) {
  
  x <- sample(1:length(tax_difference), 1, replace = FALSE)
  y <- 1:length(tax_difference) == x
  
} else {
  
  y <- tax_difference == min(tax_difference)
  
}

id.out <- suitable_equations[[y]][["equation_id"]]

# calculate the mass from the length

# load the equation and data input databases
var.dat <- readxl::read_xlsx(here("raw_data/variable_input_data.xlsx"))
var.dat <- var.dat[!is.na(var.dat$equation_id),]

var.in <- var.dat[var.dat$equation_id == id.out,]

un_var <- unique(var.in[["variable"]])
for(i in 1:length(un_var)) {
  
  x <- var.in[var.in[["variable"]] == un_var[i], ]
  
  assign(x = un_var[i], value = length_in)
  
}

# implement the equation
parsed_eq <- parse(text = equ.dat[equ.dat[["equation_id"]] == id.out,][["equation"]] )
eval(parsed_eq)

# return this as a variable
c(suitable_equations, mass_out_g = eval(parsed_eq))

### END

