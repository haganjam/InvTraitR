
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
source("companion_scripts/helper-plot-theme.R")

# load the taxon databases from the different taxonomic backbones
tax_dat <-
  lapply(c("col", "gbif", "itis"), function(x) {
    y <- readRDS(file = paste0(
      here::here("database"),
      "/",
      x,
      "_taxon",
      "_database.rds"
    ))

    return(y)
  })
tax_dat <- bind_rows(tax_dat)

# fix a duplicate name with slightly different spelling from the different
# taxonomic backbones
tax_dat$order <- ifelse(tax_dat$order == "Unionoida", "Unionida", tax_dat$order)

# summarise the taxon data
tax_sum <- 
  tax_dat %>%
  group_by(db_source, group1, group2, order) %>%
  summarise(
    n_family = length(unique(family)),
    n_taxa = length(unique(scientificName))
  ) %>%
  ungroup() %>%
  group_by(group1, group2, order) %>%
  summarise(n_family_m = mean(n_family),
            n_family_sd = sd(n_family),
            n_taxa_m = mean(n_taxa),
            n_taxa_sd = sd(n_taxa)) %>%
  mutate(group2 = ifelse(is.na(group2) | group2 == "NA", "Other", group2),
         order = ifelse(is.na(order)  | order == "NA", "Other", order))

ggplot(data = tax_sum, 
         aes(x = order, y = n_family_m, fill = n_taxa_m)) +
  geom_col() +
  facet_grid(~group1, 
             scales = "free_x",
             space = "free_x", 
             switch = "x") +
  theme_meta() +
  xlab(NULL) +
  ylab("Number of families") +
  scale_fill_viridis_c(option = "C") +
  theme(strip.placement = "outside",                      
        strip.background = element_rect(fill = "white"),  
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.3, size = 7),
        legend.position = "top")



# describe the geographical coverage

# load the habitat data
hab.dat <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver2/habitat_database.xlsx")
head(hab.dat)

# get the unique lat-lon coordinates
hab.dat <-
  hab.dat %>%
  filter(lat_dd != "NA") %>%
  mutate(
    longitude = as.numeric(lon_dd),
    latitude = as.numeric(lat_dd)
  )

# plot the points over the natural earth map
hab.sum <-
  hab.dat %>%
  mutate(location = as.character(paste0(lat_dd, lon_dd))) %>%
  group_by(database, location) %>%
  summarise(
    N = length(unique(id)),
    lat = mean(latitude),
    lon = mean(longitude)
  )
p1 <-
  ggplot() +
  geom_sf(data = spData::world, fill = "white", col = "black", size = 0.35) +
  coord_sf(ylim = c(-70, 90), xlim = c(-170, 180)) +
  geom_jitter(
    data = hab.sum,
    mapping = aes(x = lon, y = lat, size = N),
    colour = "black",
    fill = "cornflowerblue",
    width = 0.05,
    shape = 21
  ) +
  ylab(NULL) +
  xlab(NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom"
  )
p1

ggsave(
  filename = "figures/map.png", dpi = 400,
  units = "cm", width = 20, height = 12
)
