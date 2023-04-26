
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
library(cowplot)

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

# change the Platyhelminthes and Annelida label to Platy for plotting
tax_sum$group1 <- ifelse(tax_sum$group1 == "Platyhelminthes", "Platy.", tax_sum$group1)
tax_sum$group1 <- ifelse(tax_sum$group1 == "Annelida", "Annel.", tax_sum$group1)

p1 <- 
  ggplot(data = tax_sum, 
         aes(x = order, y = n_family_m, fill = n_taxa_m)) +
  geom_col(width = 0.6) +
  facet_grid(~group1, 
             scales = "free_x",
             space = "free_x", 
             switch = "x") +
  theme_meta() +
  guides(fill = guide_colourbar(title.position = "left", 
                                title.vjust = 1.2,
                                frame.colour = "black", 
                                ticks.colour = NA,
                                barwidth = 5,
                                barheight = 0.3)) +
  xlab(NULL) +
  ylab("Number of families") +
  labs(fill = "Number of taxa") +
  scale_fill_viridis_c(option = "C") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12)) +
  theme(strip.placement = "outside",                      
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 10),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.3, size = 9)) +
  theme(legend.position = c(0.6, 1),
        plot.margin = unit(c(1.8, 1, 0.5, 0.5), "lines"),
        legend.direction="horizontal",
        legend.justification=c(1, 0), 
        legend.key.width=unit(1, "lines"), 
        legend.key.height=unit(1, "lines"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10))
plot(p1)
  
ggsave(filename = "figures/fig_A.png", p1, dpi = 400,
       units = "cm", width = 20, height = 10)


# describe the geographical coverage

# load the habitat data
hab_dat <- readRDS("database/freshwater_ecoregion_data.rds")
head(hab_dat)

# get the unique lat-lon coordinates
hab_dat <-
  hab_dat %>%
  filter(lat_dd != "NA") %>%
  mutate(
    longitude = as.numeric(lon_dd),
    latitude = as.numeric(lat_dd)
  )

# plot the points over the natural earth map
hab_sum <-
  hab_dat %>%
  mutate(location = as.character(paste0(lat_dd, lon_dd))) %>%
  group_by(database, location) %>%
  summarise(
    N = length(unique(id)),
    lat = mean(latitude),
    lon = mean(longitude)
  )

p2 <-
  ggplot() +
  geom_sf(data = spData::world, fill = "white", col = "black", size = 0.35) +
  coord_sf(ylim = c(-55, 80), xlim = c(-150, 170)) +
  geom_jitter(
    data = hab_sum,
    mapping = aes(x = lon, y = lat, size = N),
    colour = "#ec7853",
    alpha = 0.75,
    width = 0.05,
    shape = 16
  ) +
  ylab(NULL) +
  xlab(NULL) +
  scale_size_continuous(range = c(1, 3.5)) +
  guides(colour = guide_legend(override.aes = list(size = 2.5,
                                                   alpha = 1,
                                                   shape = 16))) +
  theme_meta() +
  theme(
    axis.text = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key = element_rect(fill = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
  )
plot(p2)

ggsave(filename = "figures/fig_B.png", p2, dpi = 400,
       units = "cm", width = 15, height = 10)

# load the habitat metadata
hab_meta <- readRDS("database/freshwater_ecoregion_metadata.rds")
head(hab_meta)

# what are the unique realms?
unique(hab_meta$realm)

# bind the tax_dat to the hab_dat
tax_hab <- 
  full_join(tax_dat, 
            dplyr::select(hab_dat, id, realm), by = "id")

# count the unique orders per realm
N <- length(unique(paste0(tax_dat$group1, tax_dat$order)))

tax_hab <- 
  tax_hab %>%
  mutate(group1_order = paste0(group1, order)) %>%
  group_by(realm) %>%
  summarise(n = length(unique(group1_order))/N) %>%
  ungroup() %>%
  filter(!is.na(realm))

# add the missing 
tax_hab <- 
  bind_rows(tax_hab,
            hab_meta %>%
              dplyr::select(realm) %>%
              distinct() %>%
              filter( !(realm %in% tax_hab$realm) ) %>%
              mutate(n = 0)
            )

# plot these data  
p3 <- 
  ggplot(data = tax_hab,
         mapping = aes(x = realm, y = n)) +
  geom_col(width = 0.5, fill = "#ec7853") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab(NULL) +
  ylab("Proportion orders represented") +
  theme_meta()
plot(p3)

ggsave(filename = "figures/fig_B1.png", p3, dpi = 400,
       units = "cm", width = 6, height = 10)

### END
