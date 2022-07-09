
# Search the database

# load relevant spatial libraries
library(sp)
library(sf)
library(raster)
library(dplyr)
library(here)
library(igraph)
library(bdc)
library(taxadb)

# write functions for:

# 1. cleaning the input data target taxa robustly using the bdc package

# - data.frame scale

# 2. get habitat data for the target taxa

# - data.frame scale

# 3. calculate taxonomic distance habitat for each target taxa

# - taxa by taxa

# - think about where the taxa come from... what if we implement different traits?

# 4. merge all this information to make a decision

#'
#' @title extract_genus()
#' 
#' @description If a binomial taxa is supplied, extract the genus name only
#' 
#' @details This method works almost exclusively with genera and not species. This is because
#' once the genus is known, the relationship among the species within that genus is known. It
#' makes the method much more computationally efficient. This function is used to then extract
#' the genus from a species name. The genus name is used to place the species in the correct
#' taxonomic framework 
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param binomial - binomial character string separated by a space (e.g. "Loxodonta loxodonta")
#' 
#' @return string with the genus name
#' 

Extract_Genus <- function(binomial) {
  
  # test if the input is a string
  test_1 <- function(x) {
    
    assertthat::is.string(x)
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(deparse(call$x), " is not a character string")
    
  }
  
  assertthat::assert_that(test_1(x = binomial))
  
  # test if the input is of length 1
  test_2 <- function(x) {
    
    assertthat::are_equal(length(x), 1)
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is not of length = 1")
    
  }
  
  assertthat::assert_that(test_2(x = binomial))
  
  # check if input has any special characters
  test_3 <- function(x) {
    
    !grepl(pattern = '[[:punct:]]', x = x)
    
  }
  
  assertthat::on_failure(test_3) <- function(call, env){
    
    paste0(deparse(call$x), " contains special characters")
    
  }
  
  assertthat::assert_that(test_3(x = binomial))
  
  # split the binomial into separate parts
  binomial.1st <- unlist( strsplit(x = binomial, split = " ", fixed = TRUE) )
  
  # calculate the length of the split object
  binomial.l <- length(binomial.1st)
  
  # if the resulting object has a length of greater than 1, then extract first element
  if (binomial.l > 1) {
    
    binomial <- binomial.1st[1]
    
  } 
  
  # add a word count attribute
  attr( binomial, "n") <- binomial.l
  
  # return the modified name
  return(binomial)
  
}

Extract_Genus(binomial = "Loxodonta africanus")


#'
#' @title Get_tax_dist()
#' 
#' @description Get taxonomic distances of target names relative to the taxa databases
#' 
#' @details This function takes a data.frame with target_taxa and searches a set of taxonomic
#' graphs from gbif, itis and col to find taxa that are closely related. The names are first
#' cleaned and harmonised using the bdc package. Then, prep-prepared graphs that
#' contain taxonomic information are searched and distances between taxa in our database
#' and the target taxon are identified.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param data - data.frame with a column containing target taxon names
#' @param database - taxonomic database to use: gbif, itis, col
#' @param target_taxon - character string with the column name containing the taxon names
#' 
#' @return data.frame with target taxon names, database names and relevant taxonomic distances
#' 

# write a function to get the taxonomic distances for database entries
Get_tax_dist <- function(data, database = "gbif", target_taxon) {
  
  t.dat <- 
    data %>%
    rename(target_taxon = target_taxon)
  
  # clean the names for typos etc.
  cleaned_names <- bdc_clean_names(sci_names = t.dat[["target_taxon"]], save_outputs = FALSE)
  
  # write some code to remove the output file
  unlink("Output", recursive=TRUE)
  
  # check if any names were changed
  if ( !any(cleaned_names$scientificName != cleaned_names$names_clean) ) {
    message("No names were changed")
  }
  
  # replace the names in tax.dat with these cleaned names
  t.dat$target_taxon <- cleaned_names$names_clean
  
  t.dat <- 
    t.dat %>%
    mutate(row_id = 1:n()) %>%
    dplyr::select(row_id, target_taxon)
  
  # create the local database
  td_create(
    provider = database,
    overwrite = FALSE)
  
  # harmonise the names to the gbif database
  harm.tax <- 
    bdc_query_names_taxadb(sci_name = t.dat$target_taxon,
                           db = database,
                           rank_name = "Animalia",
                           rank = "kingdom"
    )
  
  # write some code to remove the output file
  unlink("Output", recursive=TRUE)
  
  harm.tax <- 
    harm.tax %>%
    mutate(db_taxon_higher_rank = ifelse(is.na(order) & is.na(family), NA, 
                                         ifelse(is.na(order) & !is.na(family), "family", "order") ) ) %>%
    mutate(db_taxon_higher = ifelse(is.na(order) & is.na(family), NA, 
                                    ifelse(is.na(order), family, order) ) ) %>%
    mutate(db_higher_rank_source = database) %>%
    mutate(row_id = 1:n()) %>%
    dplyr::select(row_id, original_search, scientificName, acceptedNameUsageID, db_higher_rank_source, db_taxon_higher_rank, db_taxon_higher)
  
  # remove the names that we were not able to resolve
  harm.tax <- 
    harm.tax %>%
    filter(!(is.na(scientificName) |is.na(db_taxon_higher_rank) | is.na(db_taxon_higher) ) ) %>%
    rename(target_taxon = original_search)
  
  # join these data to the tax.dat data
  t.clean <- right_join(t.dat, harm.tax, by = c("row_id", "target_taxon") )
  
  # check that the join worked correctly
  if ( nrow(harm.tax) == nrow(t.clean) ) {
    message("Join worked correctly")
  }
  
  # remove the row_id column
  t.clean <- 
    t.clean %>%
    dplyr::select(-row_id)
  
  # check if any of these names could be found and if not, return a blank data.frame
  if(nrow(t.clean) == 0) {
    
    t.dat$db_source <- database
    return( t.dat[, c("target_taxon", "db_source") ] )
    
  }
  
  dist.list <- vector("list", length = nrow(t.clean))
  for(i in 1:nrow(t.clean)) {
    
    # extract the target name
    targ.higher <- t.clean[i,][["db_taxon_higher"]]
    clean.targ <- t.clean[i,][["scientificName"]]
    targ <- t.clean[i,][["target_taxon"]]
    
    # check if there are any relevant names to check
    if (is.na(targ.higher) | is.na(clean.targ) ) {
      
      dist.df <- 
        tibble(target_taxon = targ, 
               clean_target_name = clean.targ,
               database_name = NA,
               tax_distance = NA)
      
      dist.list[[i]] <- dist.df
      
    } else {
      
      # get the correct database
      if(database == "gbif") {
        
        tax.graph <- gbif_htm[which(names(gbif_htm) == targ.higher)][[1]]
        v.x <- V(tax.graph)
        
        tax.td <- 
          gbif_td %>%
          filter(db_taxon_higher == targ.higher)
        
      } else if (database == "itis") {
        
        tax.graph <- itis_htm[which(names(itis_htm) == targ.higher)][[1]]
        v.x <- V(tax.graph)
        
        tax.td <- 
          itis_td %>%
          filter(db_taxon_higher == targ.higher)
        
      } else if (database == "col") {
        
        tax.graph <- col_htm[which(names(col_htm) == targ.higher)][[1]]
        v.x <- V(tax.graph)
        
        tax.td <- 
          col_td %>%
          filter(db_taxon_higher == targ.higher)
        
      } else {
        
        stop("Choose appropriate database: gbif, itis, col")
        
      }
      
      # get taxonomic distance between target name and database names
      tax.dist <- 
        lapply(tax.td$scientificName, function(y) {
          
          # extract genus for species-level names
          tn <- extract_genus(clean.targ)
          sn <- extract_genus(y)
          
          x <- 
            distances(tax.graph, 
                      v.x[which(attr(v.x, "names") == tn)],
                      v.x[which(attr(v.x, "names") == sn)],
                      mode = c("all"),
                      algorithm = c("bellman-ford"))
          
          if(length(x) == 0) {
            x <- 0
          } else {
            x <- x[[1]]
          }
          
          # extra distance for species level: 0.5
          sp.l <- sum(ifelse(c(attr(tn, "n"), attr(sn, "n")) > 1, 0.5, 0))
          
          dist.df <- 
            tibble(target_taxon = targ, 
                   clean_target_name = clean.targ,
                   database_name = y,
                   tax_distance = x + sp.l
            )
          
          return(dist.df)
          
        } )
      
      # bind into a data.frame
      tax.dist <- bind_rows(tax.dist, .id = "id_rep")
      
      # add a data.base identifier
      tax.dist$db_source <- database
      
      tax.dist <- 
        tax.dist %>%
        dplyr::select(target_taxon, db_source, clean_target_name, id_rep, database_name, tax_distance)
      
      dist.list[[i]] <- tax.dist
      
    }
    
  }
  
  dist.list <- bind_rows(dist.list)
  
  return(dist.list)
  
}


# load the taxon data
if (!exists("itis_td")) {
  itis_td <- readRDS(file = here("database/itis_taxon_database.rds"))
}

if (!exists("itis_htm")) {
  itis_htm <- readRDS(file = here("database/itis_higher_taxon_matrices.rds"))
}

if (!exists("gbif_td")) {
  gbif_td <- readRDS(file = here("database/gbif_taxon_database.rds"))
}

if (!exists("gbif_htm")) {
  gbif_htm <- readRDS(file = here("database/gbif_higher_taxon_matrices.rds"))
}

if (!exists("col_td")) {
  col_td <- readRDS(file = here("database/col_taxon_database.rds"))
}

if (!exists("col_htm")) {
  col_htm <- readRDS(file = here("database/col_higher_taxon_matrices.rds"))
}

# load the freshwater habitat map data
if (!exists("fw_er")) {
  fw_er <- readRDS(file = here("database/freshwater_ecoregion_data.rds"))
}

# load the freshwater habitat map
if (!exists("fw_map")) {
  fw_map <- readRDS(file = here("database/freshwater_ecoregion_map.rds"))
}

# load the equation data
if (!exists("equ.dat")) {
  equ.dat <- readRDS(file = here("database/equation_database.rds"))
}
head(equ.dat)


# test data.frame
t.dat <- tibble(name = c("Coloburiscus", "Zelandoperla"),
                life_stage = c("nymph", "adult"),
                length_mm = c(10, 20),
                lat = c(-25.18, -25),
                lon = c(147.69, 145.5)
)


# loop over all the relevant data.bases then extract the possible names
Tax_dist_selector <- function(data, target_taxon, max_taxdist = 2.5) {
  
  tax.df <- 
    lapply(c("gbif", "itis", "col"), function(x) {
    
    tax <- Get_tax_dist(data = data, 
                        database = x, 
                        target_taxon = target_taxon)
    
    return(tax)
    
  } )
  
  # bind the rows together
  tax.df <- bind_rows(tax.df)
  
  # keep the rows where taxonomic distance is low enough
  tax.df <- 
    tax.df %>%
    filter(tax_distance < max_taxdist)
  
  return(tax.df)
  
}

x <- Tax_dist_selector(data = t.dat, target_taxon = "name", max_taxdist = 2.5)


# 


