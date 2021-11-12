
# Pipeline to draw allometric equations from taxonomic (and geographic) information

# load relevant libraries
library(taxize)
library(readr)
library(here)
library(dplyr)

rm(list = ls() )

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

# we will use three databases
data.base <- "bold"

# set-up the taxonomic ranks that will be considered
tax.ranks <- c('class','subclass', 
               'superorder','order','suborder',
               'superfamily','family', 'subfamily',
               'genus', 'species') 

# write this into a data.frame
df.tax <- data.frame(rank = tax.ranks)

# add a quantitative hierarchical taxonomic information
rank_number <- c(1, (1+(1/3)), (1+(2/3)), 2, (2+(1/3)), (2+(2/3)), 3, (3+(1/3)), (3+(2/3)), (4+(2/3)))
df.tax$rank_number <- rank_number

# arguments for the function
rank.diff = 2

# set a name
x.name <- "Dytiscidae"

# figure out what to do when there are no children

# write the x.taxa id depending on which database is selected

if (data.base == "bold") {
  
  x.taxa <- get_boldid(sci = x.name, ask = FALSE)
  
} else if (data.base == "itis") {
  
  x.taxa <- get_tsn(sci_com = x.name, ask = FALSE)
  
} else if (data.base == "gbif") {
  
  x.taxa <- get_gbifid(sci = x.name, ask = FALSE)
  
}

if (attr(x.taxa, "match") == "not found") {
  
  x.df <- NULL
  
} else {
  
  x.taxa <- x.taxa[[1]]
  
  # get the taxonomic rank
  x.rank <- tax_rank(sci_id = x.taxa, db = data.base)[[1]]
  
  # for this rank, we exact the numeric rank
  x.rank.num <- df.tax[df.tax$rank == x.rank, "rank_number"]
  
  # get the rank to get children
  down.to.rank <- df.tax[ (rank_number > x.rank.num) & (rank_number < foc.rank + x.rank.num), ]
  
  # get all species down 
  x.down <- downstream(sci_id = x.taxa, db = data.base, 
                       downto = down.to.rank$rank[nrow(down.to.rank)], 
                       intermediate = TRUE)
  x.down <- x.down[[1]][[2]]
  
  # get the upwards classification
  x.class <- classification(sci_id = x.taxa, db = data.base)
  x.class <- x.class[[1]]
  
  # get all species up using the classification
  x.up <- df.tax[(rank_number > (x.rank.num - rank.diff) ) & (rank_number <= x.rank.num), "rank"]
  
  # merge the upstream and downstream species
  x.merge <- rbind(x.class[x.class$rank %in% x.up, ], bind_rows(x.down))
  
  # merge with the df.tax
  x.df <- right_join(df.tax, df.merge, by = "rank")
  
  # add a focal taxa column
  x.df$focal_taxa <- ifelse(x.df$name == x.name, 1, 0)
  
}

# add the equation label
x.df










# get database id
y <- get_boldid(sci = x.name, ask = FALSE)
y[[1]]

# get the classification
x <- classification(sci_id = y[[1]], db = "bold")
x

# get the taxonomic rank
x.rank <- tax_rank(sci_id = y[[1]], db = "bold")
x.rank[[1]]

foc.rank <- df.tax[df.tax$rank == x.rank[[1]], "rank_number"]
down.to.rank <- df.tax[ (df.tax$rank_number > foc.rank) & (df.tax$rank_number < foc.rank + rank.diff), ]
down.to.rank

df.down <- downstream(sci_id = y[[1]], db = "bold", downto = down.to.rank$rank[nrow(down.to.rank)], intermediate = TRUE)

x[1][[1]]

rank.up <- df.tax[(df.tax$rank_number > (foc.rank - 2) ) & (df.tax$rank_number <= foc.rank), "rank"]

df.merge <- rbind(filter( x[1][[1]], rank %in% rank.up), bind_rows(df.down[[1]][[2]]))
df.merge

df.x <- right_join(df.tax, df.merge, by = "rank")

# add a focal taxa column
df.x$focal_taxa <- ifelse(df.x$name == x.name, 1, 0)

df.x






y <- get_tsn(sci = "Daphnia magna", ask = FALSE)

z <- classification(sci_id = y[[1]], db = "itis")
z

left_join(df.tax, z$`85214`, by = "rank")

y <- get_gbifid(sci = x.name, ask = FALSE)
y[[1]]

u <- classification(sci_id = y[[1]], db = "gbif")
u

left_join(df.tax, u$`5741657`, by = "rank")


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





