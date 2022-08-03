
# Test the method for calculating species biomasses

# load the relevant libraries
library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

# load the use-scripts
source(here("scripts/02_function_plotting_theme.R"))
source(here("scripts/03_use_database/01_search_database_ver2.R"))

# check if a figure folder exists
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}

# path to where raw data is stored
path_rd <- "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/"

# test 1: actual body size data and actual length data from the literature

# dataset1

# load the literature test data
test1.a <- read_csv(paste0(path_rd, "test_data_Eklof_2016_Dumont_1975.csv") )
str(test1.a)

# remove cases where life-stages are NA
test1.a <- dplyr::filter(test1.a, !is.na(Life_stage), Life_stage != "NA")
names(test1.a)

# equalise the names
test1.a <- 
  test1.a %>%
  select(Reference, Taxa, Life_stage, lat, lon, Length_mm, Dry_weight_mg)

# dataset2
test1.b <- read_csv(paste0(path_rd, "test_data_Gonzalez_2008.csv") )
str(test1.b)

# remove cases where life-stages are NA
test1.b <- dplyr::filter(test1.b, !is.na(life_stage), life_stage != "NA")
names(test1.b)

# equalise the names
test1.b <- 
  test1.b %>%
  select(author_year, taxon, life_stage, lat, lon, mean_length_mm, mean_mass_g)

# dataset3
test1.c <- read_csv(paste0(path_rd, "test_data_Rudolph_2014.csv") )
str(test1.c)

# remove cases where life-stages are NA
test1.c <- dplyr::filter(test1.c, !is.na(life_stage), life_stage != "NA")
names(test1.c)

# equalise the names
test1.c <- 
  test1.c %>%
  select(author_year, taxon, life_stage, lat, lon, average_body_length_mm, average_mass_mg)

test1.dat <- 
  
  lapply(list(test1.a, test1.b, test1.c), function(x) {
  
  names(x) <- c("author_year", "taxon", "life_stage", "lat", "lon", "body_length_mm", "biomass_mg")
  
  return(x)
  
} )

# bind the rows
test1.dat <- bind_rows(test1.dat)
View(test1.dat)

# use method to get biomass data
test1.output <- 
  Get_Trait_From_Taxon(data = test1.dat, 
                       target_taxon = "taxon", 
                       life_stage = "life_stage", 
                       body_size = "body_length_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       workflow = "workflow2",
                       trait = "equation", 
                       max_tax_dist = 4,
                       gen_sp_dist = 0.5
  )

# get names that were not found
test1.output %>%
  filter(is.na(id)) %>%
  pull(taxon) %>%
  unique()

# remove rows where the dry-biomass is not there
test1.output <- 
  test1.output %>%
  filter(!is.na(dry_biomass_mg) )

# check how many unique taxa are left
length(unique(test1.output$taxon))

# null model i.e. randomly chosen equations

# load the equation data
equ.dat <- readRDS(file = paste0(here::here("database"), "/", "equation", "_database.rds"))
head(equ.dat)

equ.dat <- 
  equ.dat %>%
  select(equation_id, equation)

equ.null <- 
  lapply(1:1000, function(y) {
  
  sapply(test1.output$body_length_mm, function(x) {
    
    var1 <- x
    equ <- equ.dat[sample(1:nrow(equ.dat), 1),][["equation"]]
    
    eval(parse(text = equ))
    
  })
  
})  

# combine into a matrix
equ.null.m <- do.call("rbind", equ.null)

# calculate percentile intervals
equ.null.m <- apply(equ.null.m, 2, function(x) quantile(x = x, c(0.025, 0.975)) )

# add the percentile intervals to the test data
test1.output$lower_PI <- apply(equ.null.m, 2, function(x) x[1] )
test1.output$upper_PI <- apply(equ.null.m, 2, function(x) x[2] )

# plot dry weight inferred versus actual dry weight
p1 <- 
  ggplot() +
  geom_ribbon(data = test1.output,
              mapping = aes(x = log10(biomass_mg),
                            ymin = log10(lower_PI),
                            ymax = log10(upper_PI)),
              alpha = 0.1) +
  geom_point(data = test1.output,
             mapping = aes(x = log10(biomass_mg), 
                           y = log10(dry_biomass_mg), colour = author_year),
             alpha = 0.5) +
  ylab("Estimated dry biomass (mg, log10)") +
  xlab("Measured dry biomass (mg, log10)") +
  geom_abline(intercept = 0, slope = 1, 
              colour = "#ec7853", linetype = "dashed", size = 1) +
  scale_colour_viridis_d(option = "C", begin = 0, end = 0.9) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)

# observed correlation
cor.test(test1.output$biomass_mg, test1.output$dry_biomass_mg)

# null correlations
null_cor <- sapply(equ.null, function(x) cor(test1.output$biomass_mg, x) )
quantile(null_cor, c(0.025, 0.975))

# calculate percentage
test1.output <- 
  test1.output %>%
  mutate(error_perc = (abs(biomass_mg - dry_biomass_mg)/biomass_mg)*100 )

# calculate the error for each randomisation
equ.null.error <- 
  
  lapply(equ.null, function(x) {
  
  (abs(test1.output$biomass_mg - x)/test1.output$biomass_mg)*100
  
} )

# unlist all randomisations
equ.null.error <- unlist(equ.null.error)

# plot the distributions of the true error and the randomised error
error_hist <- bind_rows(tibble(Method = "FreshInvTraitR",
                               error_perc = test1.output$error_perc),
                        tibble(Method = "Random",
                               error_perc = equ.null.error))

# get the 90 % quantiles of the error for both
p2 <- 
  error_hist %>%
  group_by(Method) %>%
  mutate(low_PI = quantile(error_perc, 0.025),
         upper_PI = quantile(error_perc, 0.975)) %>%
  filter(error_perc > low_PI, error_perc < upper_PI) %>%
  
  ggplot(data = .,
       mapping = aes(x = log( error_perc ), colour = Method, fill = Method)) +
  geom_density(alpha = 0.6) +
  scale_colour_manual(values = c("#0c1787", "#fadb25")) +
  scale_fill_manual(values = c("#0c1787", "#fadb25")) +
  xlab("Absolute deviation (%, loge)") +
  ylab("Density") +
  geom_vline(xintercept = mean( log( test1.output$error_perc) ), 
             linetype = "dashed", size = 0.5,
             colour = "#0c1787") +
  geom_vline(xintercept = mean(log(equ.null.error) ), 
             linetype = "dashed", size = 0.5,
             colour = "#fadb25") +
  theme_meta()
plot(p2)

p12 <- ggarrange(p1, p2, ncol = 2,nrow = 1,
                 labels = c("a", "b"),
                 widths = c(1, 1.3),
                 font.label = list(size = 11, color = "black", face = "plain")
                 )
plot(p12)

# export Fig. 4
ggsave(filename = here("figures/fig_4.pdf"), plot = p12, 
       units = "cm", width = 20, height = 10, dpi = 300)

# what is this like on an absolute scale? Ten-fold difference in average error
mean(test1.output$error_perc)
sd(test1.output$error_perc)
mean(equ.null.error)
sd(equ.null.error)

max(test1.output$error_perc)

# plot how the error changes with taxonomic distance
test1.output.s <- 
  test1.output %>%
  group_by(tax_distance) %>%
  summarise(mean_error = mean(error_perc, na.rm = TRUE),
            sd_error = sd(error_perc, na.rm = TRUE))

test1.output %>%
  filter(tax_distance <= 1) %>%
  summarise(mean_error = mean(error_perc),
            sd_error = sd(error_perc), 
            n = length(unique(taxon)))

p3 <- 
  ggplot() +
  ylab("Absolute deviation (%)") +
  xlab("Taxonomic distance") +
  geom_quasirandom(data = test1.output,
              mapping = aes(x = tax_distance, y = error_perc), 
              width = 0.15, alpha = 0.15) +
  geom_point(data = test1.output.s,
             mapping = aes(x = tax_distance, y = mean_error),
             size = 2.5) +
  geom_errorbar(data = test1.output.s,
                mapping = aes(x = tax_distance, 
                              ymin = mean_error - sd_error,
                              ymax = mean_error + sd_error),
                width = 0.1) +
  theme_meta()
plot(p3)

# export Fig. 5
ggsave(filename = here("figures/fig_5.pdf"), plot = p3, 
       units = "cm", width = 10, height = 10, dpi = 300)


# test 2: equations selected by expert

# load data from Dolmans (2022, unpublished)
test2.a <- read_csv("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data_Dolmans_2022.csv")
test2.a$lat <- NA
test2.a$lon <- NA

# check the summary statistics
summary(test2.a)
str(test2.a)

# add an author_year column
test2.a$author_year <- "Dolmans_2022"

# equalise the columns
names(test2.a)

test2.a <- 
  test2.a %>%
  select(author_year, Focal_taxon, Life_stage, lat, lon, length_mm, Biomass_mg)

# load the data from Gorman et al. (2017, Nature Climate Change)
test2.b <- read_csv("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database_ver2/test_data_Gorman_2017.csv")

# check the summary statistics
summary(test2.b)
str(test2.b)

# add an author_year column
test2.b$author_year <- "Gorman_2017"

# equalise the columns
names(test2.b)

test2.b <- 
  test2.b %>%
  select(author_year, species, life_stage, lat_dd, lon_dd, length, mass)

# make the names the same between the two datasets
test2.dat <- 
  lapply(list(test2.a, test2.b), function(x) {
  
  names(x) <- c("author_year", "taxon", "life_stage", "lat", "lon", "length_mm", "mass_mg")
  
  return(x)
  
} )

# collapse into a single data.frame
test2.dat <- bind_rows(test2.dat)
View(test2.dat)

# test the method
test2.output <- 
  Get_Trait_From_Taxon(data = test2.dat, 
                       target_taxon = "taxon", 
                       life_stage = "life_stage", 
                       body_size = "length_mm",
                       latitude_dd = "lat", 
                       longitude_dd = "lon",
                       trait = "equation", 
                       workflow = "workflow2",
                       max_tax_dist = 3,
                       gen_sp_dist = 0.5
  )

# get names that were not found
test2.output %>%
  filter(is.na(id)) %>%
  pull(taxon) %>%
  unique()

# remove rows where the weight is not there
test2.output <- 
  test2.output %>%
  filter(!is.na(dry_biomass_mg) )

test2.output <- 
  test2.output %>%
  select(author_year, taxon, life_stage, length_mm, scientificName, db.scientificName, 
         tax_distance, id,
         mass_mg, dry_biomass_mg, life_stage_match) %>%
  mutate(dry_biomass_mg = round(dry_biomass_mg, 5),
         error_perc = (abs(mass_mg - dry_biomass_mg)/mass_mg)*100)

# check how many unique taxa there are
length(unique(test2.output$taxon))

p4 <- 
  ggplot(data = test2.output,
       mapping = aes(x = log10(mass_mg), 
                     y = log10(dry_biomass_mg),
                     colour = author_year) ) +
  geom_point() +
  ylab("Estimated dry biomass (mg, log10)") +
  xlab("Expert dry biomass (mg, log10)") +
  geom_abline(intercept = 0, slope = 1, 
              colour = "#ec7853", linetype = "dashed", size = 1) +
  theme_meta()
plot(p4)

# check the large errors
test2.output %>%
  filter(error_perc > quantile(test2.output$error_perc, 0.95)) %>%
  View()

ggsave(filename = here("figures/fig_6.pdf"), plot = p4, 
       units = "cm", width = 10, height = 10, dpi = 300)

cor.test(log10(test2.output$mass_mg), log10(test2.output$dry_biomass_mg) )

# what about the error?
mean(test2.output$error_perc)
sd(test2.output$error_perc)
round(range(test2.output$error_perc), 3)

test2.output %>%
  filter(error_perc < 200) %>%
  summarise(mean = mean(error_perc),
            sd = sd(error_perc))

### END
