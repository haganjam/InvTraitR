
# Graph-based approach to similarity

# load relevant libraries
library(igraph)

# good link on using igraph
# https://kevintshoemaker.github.io/NRES-746/graph.theory.html

# graph-based approach to similarity
df_ver <- data.frame(name = c(1, 2, 3, 4, 5))
df_edge <- data.frame(from = c(1, 2, 2, 2), to = c(2, 3, 4, 5),
                      weights = c(1, 0.5, 0.5, 0.5))

g <- graph_from_data_frame(df_edge, directed=FALSE, vertices=df_ver)
print(g, e=TRUE, v=TRUE)
d <- distances(
  g,
  v = V(g),
  to = V(g),
  mode = c("all"),
  weights = df_edge$weights,
  algorithm = c("bellman-ford")
)

# example code

# set the input name
equ.name <- "Tanytarsus"

# get taxon id: database specific function
#
y <- get_taxon_id(database_function = "gbif", taxon_name = equ.name, ask_or_not = FALSE, tries = 5)
#

# get upwards classification: general function
#
y.c <- classification(y)[[1]]
y.c$ranknum <- 1:nrow(y.c)

# get target rank
equ.rank <- y.c[y.c$name == in.name,][["ranknum"]]

# get order rank
ord.rank <- y.c[y.c$rank == "order", ][["ranknum"]]

if (equ.rank < ord.rank) {
  
  message("equation above the rank of order")
  equ.suitable <- FALSE
  
} else { equ.suitable <- TRUE }
#

# stop process if the equation is not suitable
#

# get taxon id: database specific function
high_id <- get_taxon_id(database_function = "gbif", 
                        taxon_name =  y.c[y.c$rank == "order", ][["name"]], 
                        ask_or_not = TRUE, tries = 5)

x <- get_gbifid2(sci = y.c[y.c$rank == "order", ][["name"]], ask = TRUE)


# deal with the subscript out of bounds error
# this seems to happen when there are too many taxa
# will need to go up to the next taxnonomic level as an error fix

z <- downstream(sci_id = high_id[[1]], downto = "genus", db = "gbif", intermediate = TRUE)
z$`3343`$intermediate

View(z[[1]])

# set appropriate taxonomic weights
tax.ranks <- c('class','subclass', 
               'superorder','order','suborder',
               'superfamily','family', 'subfamily',
               'genus', 'species') 

u <- 
  z$`127917`$intermediate %>%
  bind_rows()
View(u)

# sequential sorting process
x <- 
  u %>%
  select(from = parentname, to = taxonname)
head(x)

g <- graph_from_data_frame(d = x, directed=FALSE)
print(g, e=TRUE, v=TRUE)
plot(g)
d <- distances(
  g,
  v = V(g),
  to = V(g),
  mode = c("all"),
  algorithm = c("bellman-ford")
)
d[1:50, 1:50]

# taxon name database
db.name <- "Pediciini"
db.name

# target name
tar.name <- "Limnophila costata"
tar.name

# this calculation across all the different distance matrices
d[which(row.names(d) == db.name), which(colnames(d) == tar.name) ]

# rank of all the different names



 ### END

