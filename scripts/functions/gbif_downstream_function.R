
# get downstream taxa gbif
ord.id <- "1005"
ord.name <- "Flosculariaceae"
ord.rank <- 4

rank.vec <- y.c[y.c$ranknum > ord.rank, ][["rank"]]

library(taxize)
p <- downstream(sci_id = ord.id, downto = rank.vec[1], db = "gbif", intermediate = TRUE)
p1 <- p[[1]][[1]]
p1$parentname <- ord.name
p1

o <- downstream(sci_id = p1$key[1], downto = rank.vec[2], db = "gbif", intermediate = TRUE)
o1 <- o[[1]][[1]]
o1

p1 <- p[[1]]$intermediate[[1]]
p1$parentname <- ord.name
p1