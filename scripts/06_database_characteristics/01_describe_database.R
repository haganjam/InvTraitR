
# Describe the characteristics of the database

# load the relevant libraries
library(here)
library(dplyr)
library(readr)
library(igraph)
library(ggplot2)
library(sp)
library(sf)
library(raster)

# load the function plotting theme
source(here("scripts/02_function_plotting_theme.R"))

# describe the taxonomic coverage

# need to decide on some higher classification scheme so that we can quantify the
# coverage... How detailed should it be?

# load the taxon databases from the different taxonomic backbones
tax.dat <- 
  
  lapply(c("col", "gbif", "itis"), function(x) {
  
  y <- readRDS(file = paste0(here::here("database"), "/", x, "_taxon", "_database.rds"))
  y[y$group1 == "Acari", ][["group1"]] <- "Collembola"
  y[y$group2 == "Cladocera", ][["group2"]] <- "Diplostraca"
  
  return(y)
  
} )
tax.dat <- bind_rows(tax.dat)
View(tax.dat)

# select the relevant columns for describing the taxonomic distribution in the database
names(tax.dat)
tax.dat <- 
  tax.dat %>%
  dplyr::select(group1, group2, db_taxon, db_higher_rank_source, id, scientificName, family, db_taxon_higher_rank, db_taxon_higher)
View(tax.dat)

# summarise the taxon data
tax.dat %>%
  group_by(group1) %>%
  summarise(n_family = length(unique(family)),
            n_taxa = length(unique(scientificName)))

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  group_by(group1) %>%
  summarise(n_order = length(unique(db_taxon_higher)))

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  View()

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  pull(db_taxon_higher) %>%
  unique()

tax.dat %>%
  filter(db_taxon_higher_rank == "order") %>%
  group_by(db_higher_rank_source, group1, group2, db_taxon_higher) %>%
  summarise(n_family = length(unique(family))) %>%
  View()

tax.dat %>%
  filter(group1 == "Acari") %>%
  View()


# describe the geographical coverage

# load the habitat data
hab.dat <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver2/habitat_database.xlsx")
head(hab.dat)

# get the unique lat-lon coordinates
hab.dat <- 
  hab.dat %>%
  filter( lat_dd != "NA" ) %>%
  mutate(longitude = as.numeric(lon_dd),
         latitude = as.numeric(lat_dd))

# plot the points over the natural earth map
hab.sum <- 
  hab.dat %>%
  mutate(location = as.character(paste0(lat_dd, lon_dd )) ) %>%
  group_by(database, location) %>%
  summarise(N = length(unique(id)),
            lat = mean(latitude),
            lon = mean(longitude))
p1 <- 
  ggplot() + 
  geom_sf(data = spData::world, fill = "white", col = "black", size = 0.35) +
  coord_sf(ylim = c(-70,90), xlim = c(-170, 180)) +
  geom_jitter(data = hab.sum,
              mapping = aes(x = lon, y = lat, size = N),
              colour = "black",
              fill = "cornflowerblue",
              width = 0.05,
              shape = 21) +
  ylab(NULL) +
  xlab(NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom")
p1

ggsave(filename = "figures/map.png", dpi = 400,
       units = "cm", width = 20, height = 12)

### END
