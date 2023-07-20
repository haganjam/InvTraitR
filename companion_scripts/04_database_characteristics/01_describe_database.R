
# Describe the characteristics of the database

# load the relevant libraries
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
tax_dat <- dplyr::bind_rows(tax_dat)

# fix a duplicate name with slightly different spelling from the different
# taxonomic backbones
tax_dat$order <- ifelse(tax_dat$order == "Unionoida", "Unionida", tax_dat$order)

# how many unique taxon names are represented in the equation database
length(unique(tax_dat$db_taxon))
length(unique(tax_dat$family))
length(unique(tax_dat$order))

# dplyr::summarise the taxon data
tax_sum <- 
  tax_dat |>
  dplyr::group_by(db_source, group1, group2, order) |>
  dplyr::summarise(
    n_family = length(unique(family)),
    n_taxa = length(unique(scientificName))) |>
  dplyr::ungroup() |>
  dplyr::group_by(group1, group2, order) |>
  dplyr::summarise(n_family_m = mean(n_family),
                   n_family_sd = sd(n_family),
                   n_taxa_m = mean(n_taxa),
                   n_taxa_sd = sd(n_taxa)) |>
  dplyr::mutate(group2 = ifelse(is.na(group2) | group2 == "NA", "Other", group2),
                order = ifelse(is.na(order)  | order == "NA", "Other", order))

# sum up by group1
tax_sum |>
  dplyr::group_by(group1) |>
  dplyr::summarise(n_taxa = sum(n_taxa_m))

# change the Platyhelminthes and Annelida label to Platy for plotting
tax_sum$group1 <- ifelse(tax_sum$group1 == "Platyhelminthes", "Platy.", tax_sum$group1)
tax_sum$group1 <- ifelse(tax_sum$group1 == "Annelida", "Annel.", tax_sum$group1)

p1 <- 
  ggplot(data = tax_sum, 
         aes(x = order, y = n_family_m, fill = group1)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.1)) +
  geom_text(aes(x = order, y = n_family_m + 0.5, label = round(n_taxa_m, 0) ),
            size = 3) +
  facet_grid(~group1, 
             scales = "free_x",
             space = "free_x", 
             switch = "x") +
  theme_meta() +
  xlab(NULL) +
  ylab("Number of families") +
  scale_fill_manual(values = wesanderson::wes_palette("Darjeeling1", n = 6, "continuous")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 13), breaks = seq(0, 12, 2)) +
  theme(strip.placement = "outside",                      
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size = 11),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.3, size = 10),
        axis.ticks.x = element_blank(),
        legend.position = "none")
plot(p1)
  
ggsave(filename = "figures/fig_2.png", p1, dpi = 400,
       units = "cm", width = 20, height = 10)


# describe the geographical coverage

# load the habitat data
hab_dat <- readRDS("database/freshwater_ecoregion_data.rds")
head(hab_dat)

# get the unique lat-lon coordinates
hab_dat <-
  hab_dat |>
  dplyr::filter(lat_dd != "NA") |>
  dplyr::mutate(
    longitude = as.numeric(lon_dd),
    latitude = as.numeric(lat_dd)
  )

# plot the points over the natural earth map
hab_sum <-
  hab_dat |>
  dplyr::mutate(location = as.character(paste0(lat_dd, lon_dd))) |>
  dplyr::group_by(database,realm, location) |>
  dplyr::summarise(
    N = length(unique(id)),
    lat = mean(latitude),
    lon = mean(longitude)
  )

# how many unique geographical locations?
length(unique(hab_sum$location))

# set-up the colour palette
pal1 <- c("#556A5B", "#50A45C", "#F2AD00", "#F69100", "#5BBCD6", "#C49647", "#FF0000")

p2 <-
  ggplot() +
  geom_sf(data = spData::world, fill = "white", col = "black", size = 0.35) +
  coord_sf(ylim = c(-55, 80), xlim = c(-150, 170)) +
  geom_jitter(
    data = hab_sum |> dplyr::mutate(realm = factor(realm)),
    mapping = aes(x = lon, y = lat, size = N, fill = realm ),
    colour = "black",
    alpha = 0.6,
    width = 0.05,
    shape = 21
  ) +
  ylab(NULL) +
  xlab(NULL) +
  scale_fill_manual(values = pal1[c(2, 4, 5, 7)]) +
  scale_size_continuous(range = c(2, 5.5)) +
  guides(size = guide_legend(override.aes = list(alpha = 0.5,
                                                 shape = 21,
                                                 colour = "black",
                                                 fill = "white",
                                                 linewidth = 2)),
         fill = "none") +
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

ggsave(filename = "figures/fig_3a.png", p2, dpi = 400,
       units = "cm", width = 15, height = 10)

# load the habitat metadata
hab_meta <- readRDS("database/freshwater_ecoregion_metadata.rds")
head(hab_meta)

# what are the unique realms?
unique(hab_meta$realm)

# bind the tax_dat to the hab_dat
tax_hab <- 
  dplyr::full_join(tax_dat, 
                   dplyr::select(hab_dat, id, realm), by = "id")

# count the unique orders per realm
N <- length(unique(paste0(tax_dat$group1, tax_dat$order)))

tax_hab <- 
  tax_hab |>
  dplyr::mutate(group1_order = paste0(group1, order)) |>
  dplyr::group_by(realm) |>
  dplyr::summarise(n_equ = length(unique(id)),
                   n_prop = length(unique(group1_order))/N) |>
  dplyr::ungroup() |>
  dplyr::filter(!is.na(realm))

# add the missing 
tax_hab <- 
  dplyr::bind_rows(tax_hab,
                   hab_meta |>
                     dplyr::select(realm) |>
                     distinct() |>
                     dplyr::filter( !(realm %in% tax_hab$realm) ) |>
                     dplyr::mutate(n_prop = 0)
                   )

# plot these data  
p3 <- 
  ggplot(data = tax_hab,
         mapping = aes(x = realm, y = n_prop, fill = realm)) +
  geom_col(width = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.75)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_fill_manual(values = pal1) +
  xlab(NULL) +
  ylab("Proportion orders represented") +
  theme_meta() +
  theme(legend.position = "none")
plot(p3)

ggsave(filename = "figures/fig_3b.png", p3, dpi = 400,
       units = "cm", width = 6, height = 10)

### END
