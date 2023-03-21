# test the method for calculating species biomasses

# load the relevant libraries
library(dplyr)
library(readr)
library(mobsim)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(igraph)
library(assertthat)

# load the use-scripts
source("companion_scripts/02_function_plotting_theme.R")

# load the required functions
source("R/clean_taxon_names.R")
source("R/get_habitat_data.R")
source("R/select_traits_tax_dist.R")
source("R/special_names.R")
source("R/get_trait_from_taxon.R")
source("R/helpers.R")

# check if a figure folder exists
if (!dir.exists("figures")) {
  dir.create("figures")
}

# load datasets compiled from the literature
test_names <- list.files("database/")
test_names <- test_names[grepl(pattern = "test_a_", x = test_names)]

# load the datasets and name them dat1 ... datN 
for(i in 1:length(test_names)) {
  
  x <- read_csv(paste0("database/", test_names[i]))
  assign(x = paste0("dat", i), value = x)
  
}

# check the datasets
lapply(list(dat1, dat2, dat3, dat4, dat5), head)

# extract the relevant columns from each dataset

# dat1
dat1 <- 
  dat1 %>%
  select(Reference, order, Taxa, lat, lon, Life_stage, Length_mm, Dry_weight_mg)

# dat2
dat2 <- 
  dat2 %>%
  select(Reference, order, Taxa, lat, lon, Life_stage, Length_mm, Dry_weight_mg)

# dat3
dat3 <-
  dat3 %>%
  select(author_year, order, taxon, lat, lon, life_stage, mean_length_mm, mean_mass_g)

# dat4
dat4 <- 
  dat4 %>%
  select(reference, order, taxa, lat, lon, life_stage, length_mm, dry_weight_mg)

# dat5
dat5 <- 
  dat5 %>%
  select(author_year, order, taxon, lat, lon, life_stage, average_body_length_mm, average_mass_mg)

# standardise the column names in these data
dat <- 
  
  lapply(list(dat1, dat2, dat3, dat4, dat5), function(x) {
  
  names(x) <- c("author_year", "order", "taxon", "lat", "lon", "life_stage", "length_mm", "obs_dry_biomass_mg")
  return(x)
  
})

# bind into a single large data.frame
dat <- bind_rows(dat)

# how many species do we have here
length(unique(dat$taxon))

# remove cases where life-stages are NA
dat <- dplyr::filter(dat, !is.na(life_stage), life_stage != "NA")

# check the summary statistics
summary(dat)

# check why there are negative obs dry biomass values
dat %>%
  filter(obs_dry_biomass_mg < 0 | obs_dry_biomass_mg == 0)

# filter these negative values
dat <- 
  dat %>%
  filter( !(obs_dry_biomass_mg < 0) )

# sample from these data to make sure we don't pseudoreplicate too much
head(dat)
dat <- 
  dat %>%
  group_by(author_year, taxon) %>%
  sample_n(size = ifelse(min(n()) < 5, min(n()), 5), replace = FALSE) %>%
  ungroup()

# use method to get biomass data
output <-
  get_trait_from_taxon(
    data = dat,
    target_taxon = "taxon",
    life_stage = "life_stage",
    body_size = "length_mm",
    latitude_dd = "lat",
    longitude_dd = "lon",
    workflow = "workflow2",
    trait = "equation",
    max_tax_dist = 3,
    gen_sp_dist = 0.5
  )

# get names that were not found
output %>%
  filter(is.na(id)) %>%
  pull(taxon) %>%
  unique()

# remove rows where the dry-biomass is not there
output <-
  output %>%
  filter(!is.na(dry_biomass_mg))

# check how many unique taxa are left
length(unique(output$taxon))
nrow(output)

# make a author_year - taxon combination column
output$group <- paste(output$author_year, output$taxon, sep = "_")

# what's the minimum number in the output author_year column
output %>%
  group_by(group) %>%
  summarise(n = n())

# plot dry weight inferred versus actual dry weight
p1 <-
  ggplot()+
  geom_smooth(
    data = output,
    mapping = aes(
      x = log10(obs_dry_biomass_mg),
      y = log10(dry_biomass_mg), group = group, colour = author_year),
    alpha = 1, linewidth = 0.25, method = "lm", se = FALSE) +
  geom_point(
    data = output,
    mapping = aes(
      x = log10(obs_dry_biomass_mg),
      y = log10(dry_biomass_mg), colour = author_year),
    alpha = 1, shape = 1, size = 2) +
  ylab("Estimated dry biomass (mg, log10)") +
  xlab("Measured dry biomass (mg, log10)") +
  geom_abline(
    intercept = 0, slope = 1,
    colour = "#ec7853", linetype = "dashed", linewidth = 1) +
  scale_colour_viridis_d(option = "C", begin = 0, end = 0.9) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)

# observed correlation
cor.test((output$obs_dry_biomass_mg), (output$dry_biomass_mg))

# calculate percentage for actual data
output <-
  output %>%
  mutate(error_perc = (abs(obs_dry_biomass_mg - dry_biomass_mg) / obs_dry_biomass_mg) * 100)

# null model i.e. order-level equations
order_null <- read_csv("database/test_order_level_null_model.csv")

# create a data.frame to collect the null biomass estimations
null_error <- tibble(
  tax_distance = ">4",
  obs_dry_biomass_mg = output$obs_dry_biomass_mg
)

# calculate the order-level biomass from the input body lengths
null_error$order_biomass_mg <-
  mapply(function(x, y) {
    
    var1 <- y
    equ <- parse(text = order_null[order_null[["order"]] == x, ][["equation"]])
    eval(equ)
    
  }, output$order, output$length_mm, USE.NAMES = FALSE)

cor.test(null_error$obs_dry_biomass_mg, null_error$order_biomass_mg)

# calculate percentage for order-level equations
null_error <-
  null_error %>%
  mutate(error_perc = (abs(obs_dry_biomass_mg - order_biomass_mg) / obs_dry_biomass_mg) * 100) %>%
  select(tax_distance, error_perc)

# combine the null error with actual biomass data
error_plot <-
  output %>%
  mutate(tax_distance = as.character(tax_distance)) %>%
  select(tax_distance, error_perc) %>%
  bind_rows(., null_error)

# add a row with taxonomic distance of 2
error_plot <-
  bind_rows(
    tibble(
      tax_distance = c("1.5", "2"),
      error_perc = NA
    ),
    error_plot
  )

error_plot$tax_distance <- factor(error_plot$tax_distance,
  levels = c("0", "0.5", "1", "1.5", "2", "2.5", "3", ">4")
)

# what about the average error?
mean(output$error_perc)
sd(output$error_perc)

# plot how the error changes with taxonomic distance
error_plot_sum <-
  error_plot %>%
  group_by(tax_distance) %>%
  summarise(
    mean_error = mean(error_perc, na.rm = TRUE),
    sd_error = sd(error_perc, na.rm = TRUE),
    n = n()
  )
print(error_plot_sum)

output %>%
  filter(tax_distance != ">4") %>%
  filter(as.numeric(tax_distance) <= 1) %>%
  summarise(
    mean_error = mean(error_perc),
    sd_error = sd(error_perc),
    n = length(unique(taxon))
  )

p2 <-
  ggplot() +
  geom_hline(yintercept = 50, linetype = "dashed") +
  annotate(geom = "text", x = 8.35, y = 1800, label = "c", size = 5) +
  ylab("Deviation (%)") +
  xlab("Taxonomic distance") +
  geom_vline(xintercept = 7.5) +
  geom_quasirandom(
    data = error_plot,
    mapping = aes(x = tax_distance, y = error_perc),
    width = 0.15, alpha = 0.1
  ) +
  geom_point(
    data = error_plot_sum,
    mapping = aes(x = tax_distance, y = mean_error),
    size = 2,
    colour = "red"
  ) +
  geom_errorbar(
    data = error_plot_sum,
    mapping = aes(
      x = tax_distance,
      ymin = mean_error - sd_error,
      ymax = mean_error + sd_error
    ),
    width = 0.1,
    colour = "red"
  ) +
  scale_x_discrete(breaks = c("0", "1", "2", "3", ">4")) +
  theme_meta()
plot(p2)

# combine the figures
p12 <-
  ggarrange(p1, p2,
    ncol = 2, nrow = 1, labels = c("a", "b"),
    label.y = 0.975,
    font.label = list(size = 12, color = "black", face = "plain")
  )

# export Fig. 4
ggsave(
  filename = here("figures/fig_4.pdf"), plot = p12,
  units = "cm", width = 20, height = 10, dpi = 300
)


# check the quality of the equations and the r2 value
plot(output$r2_match, output$error_perc)

# can we use the metadata information to improve estimates
names(output)

mod <- 
  output %>%
  select(order, 
         dry_biomass_mg, tax_distance,
         life_stage_match, r2_match, body_size_range_match, realm_match,
         major_habitat_type_match,
         obs_dry_biomass_mg)
head(mod)

# check the order levels
unique(mod$order)

# get the complete cases
mod <- mod[complete.cases(mod), ]
nrow(mod)

# predict the observed dry biomass using these factors
View(mod)

# remove the Odonates because they seem extreme
mod <- 
  mod %>%
  filter(order != "Odonata")

# fit a linear model

# null model
lm_null <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg), data = mod)
summary(lm_null)
AIC(lm_null)

# fit a series of linear models with different predictors

# the simplest model can be an interaction between order and estimated dry biomass
# this can then give us an error distribution for each outputted value

# but including body-size range match, r2 match, tax_distance and habitat match
# I think this will considerably improve our predictions

# order
lm1 <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg)*order, data = mod)
summary(lm1)
AIC(lm1)

# r2 of equation
lm2 <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg) + r2_match, data = mod)
summary(lm2)
AIC(lm2)

# taxonomic distance
lm3 <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg) + tax_distance, data = mod)
summary(lm3)
AIC(lm3)

# habitat match
lm4 <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg) + realm_match + major_habitat_type_match,
          data = mod)
summary(lm4)
AIC(lm4)

# body size range
lm5 <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg) + body_size_range_match,
          data = mod)
summary(lm5)
AIC(lm5)

# combine the best variables

# tax distance and r2 match
lm5 <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg) + r2_match*tax_distance*order,
          data = mod)
summary(lm5)
AIC(lm5)

# tax distance and r2 match plus body size range match
lm6 <- lm(log10(obs_dry_biomass_mg) ~ log10(dry_biomass_mg) + r2_match*tax_distance*order + body_size_range_match,
          data = mod)
summary(lm6)
AIC(lm6)

# plot these model predictions
lmx <- lm5
pred <- 10^(predict(lmx))

# plot the predictions without model corrections
plot(mod$dry_biomass_mg, mod$obs_dry_biomass_mg)
abline(0, 1, col = "red")

# plot the predictions with model corrections
plot((pred), mod$obs_dry_biomass_mg)
abline(0, 1, col = "red")

plot(pred[pred<0.05], mod$obs_dry_biomass_mg[pred<0.05] )
abline(0, 1, col = "red")

# calculate absolute error without model correction
abs_nomod <- abs(mod$dry_biomass_mg - mod$obs_dry_biomass_mg)
mean(abs_nomod)
hist(abs_nomod)
round(range(abs_nomod), 2)

mean((abs_nomod/mod$obs_dry_biomass_mg)*100)

# calculate absolute error with model correction
abs_mod <- (abs(pred - mod$obs_dry_biomass_mg))
mean(abs_mod)
hist(abs_mod)
round(range(abs_mod), 2)

mean((abs_mod/mod$obs_dry_biomass_mg)*100)




# test 1b: community-level simulations

# function to simulate vectors of species abundances from a
# log-normal distribution
# core function is from the mobsim() package
Sim_regional <- function(nspec = 10, nind = 10000, cv_abun = 1, plot = TRUE) {
  # simulate a species pool
  reg_vector <- sim_sad(
    s_pool = nspec,
    n_sim = nind,
    sad_type = "lnorm",
    sad_coef = list(cv_abund = cv_abun)
  )
  attr(reg_vector, "class") <- NULL
  names(reg_vector) <- NULL

  # add species names to these species
  specnames <- paste("sp", 1:length(reg_vector), sep = "")
  names(reg_vector) <- specnames

  # plot a histogram of the regional species pool abundance distribution
  if (plot) {
    hist(reg_vector)
  }

  # add a species number attribute to the reg_vector
  attr(reg_vector, "specnumber") <- nspec

  return(reg_vector)
}

# get a species list
sp_list <- unique(test1.output$taxon)

# set-up a parameter list
par.space <-
  expand.grid(
    REP = 1:10,
    SR = seq(6, 16, 2),
    IND = seq(50, 500, 50)
  )

nrow(par.space)

output_df <- vector("list", length = nrow(par.space))
for (i in 1:nrow(par.space)) {
  # generate a vector species abundances
  sp_abun <- Sim_regional(
    nspec = par.space[i, ][["SR"]],
    nind = rpois(n = 1, par.space[i, ][["IND"]]),
    cv_abun = 1,
    plot = FALSE
  )

  # remove species names from the vector of species abundances
  names(sp_abun) <- NULL

  # get a set of n species from the species list
  n_sp <- sample(sp_list, length(sp_abun), replace = FALSE)

  df_in <-
    test1.output %>%
    filter(taxon %in% n_sp) %>%
    group_by(taxon) %>%
    sample_n(size = 1) %>%
    ungroup() %>%
    select(taxon, tax_distance, biomass_mg, dry_biomass_mg) %>%
    arrange(biomass_mg)

  # add the abundance data to the df_in vector
  df_in$abun <- (sp_abun)

  # calculate the true and estimated biomasses
  df_in <-
    df_in %>%
    mutate(
      true_com_biomass = biomass_mg * abun,
      est_com_biomass = dry_biomass_mg * abun
    )

  # generate an output table
  output_df[[i]] <- tibble(
    total_abun = sum(sp_abun),
    SR = length(sp_abun),
    tax_distance = mean(df_in$tax_distance),
    max_tax_distance = max(df_in$tax_distance),
    true_com_biomass = sum(df_in$true_com_biomass),
    est_com_biomass = sum(df_in$est_com_biomass)
  )
}

# bind output into a data.frame
output_df <- as_tibble(do.call("rbind", output_df))

# calculate the absolute deviation (%)
output_df <-
  output_df %>%
  mutate(error_perc = (abs(true_com_biomass - est_com_biomass) / true_com_biomass) * 100)

p1 <-
  ggplot(
    data = output_df,
    mapping = aes(
      x = log10(true_com_biomass),
      y = log10(est_com_biomass), colour = SR
    )
  ) +
  geom_point(alpha = 0.5) +
  ylab("Est. community dry biomass (mg, log10)") +
  xlab("True community dry biomass (mg, log10)") +
  geom_abline(
    intercept = 0, slope = 1,
    colour = "#ec7853", linetype = "dashed", size = 1
  ) +
  scale_colour_viridis_c(option = "C", begin = 0, end = 0.9) +
  guides(color = guide_colourbar(
    title.position = "left",
    title.vjust = 1.2,
    frame.colour = "black",
    ticks.colour = NA,
    barwidth = 8,
    barheight = 0.3
  )) +
  theme_meta() +
  theme(legend.position = "top")

# what is the average error?
mean(output_df$error_perc)
sd(output_df$error_perc)

# what's the correlation
cor.test(log10(output_df$true_com_biomass), log10(output_df$est_com_biomass))

# export Fig. 5
ggsave(
  filename = here("figures/fig_5.pdf"), plot = p1,
  units = "cm", width = 10, height = 10, dpi = 300
)

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

# there is a problem with the Scatella data i.e. same mass for different lengths
test2.b <-
  test2.b %>%
  filter(species != "Scatella")

# make the names the same between the two datasets
test2.dat <-
  lapply(list(test2.a, test2.b), function(x) {
    names(x) <- c(
      "author_year", "taxon",
      "life_stage", "lat",
      "lon", "length_mm",
      "mass_mg"
    )

    return(x)
  })

# collapse into a single data.frame
test2.dat <- bind_rows(test2.dat)
View(test2.dat)

# test the method
test2.output <-
  get_trait_from_taxon(
    data = test2.dat,
    target_taxon = "taxon",
    life_stage = "life_stage",
    body_size = "length_mm",
    latitude_dd = "lat",
    longitude_dd = "lon",
    trait = "equation",
    workflow = "workflow2",
    max_tax_dist = 3.5,
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
  filter(!is.na(dry_biomass_mg))

# calculate error percentages of the measurements
test2.output <-
  test2.output %>%
  mutate(
    dry_biomass_mg = round(dry_biomass_mg, 5),
    error_perc = (abs(mass_mg - dry_biomass_mg) / mass_mg) * 100
  )

# check how many unique taxa there are
length(unique(test2.output$taxon))

# check how many data points there are
nrow(test2.output)

p1 <-
  ggplot(
    data = test2.output,
    mapping = aes(
      x = log10(mass_mg),
      y = log10(dry_biomass_mg),
      colour = author_year
    )
  ) +
  geom_point(alpha = 0.75) +
  ylab("Estimated dry biomass (mg, log10)") +
  xlab("Expert dry biomass (mg, log10)") +
  geom_abline(
    intercept = 0, slope = 1,
    colour = "#ec7853", linetype = "dashed", size = 1
  ) +
  scale_colour_manual(values = c("#0c1787", "#fadb25")) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)

# which points are outliers?
test2.output %>%
  filter(error_perc > quantile(test2.output$error_perc, 0.95)) %>%
  mutate(mass_mg = log10(mass_mg)) %>%
  View()

# summarise the error by taxonomic distance
test2.output.s <-
  test2.output %>%
  group_by(tax_distance) %>%
  summarise(
    mean_error = mean(error_perc, na.rm = TRUE),
    sd_error = sd(error_perc, na.rm = TRUE)
  )
print(test2.output.s)

p2 <-
  ggplot() +
  ylab("Deviation (%)") +
  xlab("Taxonomic distance") +
  geom_hline(yintercept = 50, linetype = "dashed") +
  geom_quasirandom(
    data = test2.output,
    mapping = aes(x = tax_distance, y = error_perc),
    width = 0.15, alpha = 0.1
  ) +
  geom_point(
    data = test2.output.s,
    mapping = aes(x = tax_distance, y = mean_error),
    size = 2,
    colour = "red"
  ) +
  geom_errorbar(
    data = test2.output.s,
    mapping = aes(
      x = tax_distance,
      ymin = mean_error - sd_error,
      ymax = mean_error + sd_error
    ),
    width = 0.1,
    colour = "red"
  ) +
  theme_meta()
plot(p2)

p12 <- ggarrange(p1, p2,
  ncol = 2, nrow = 1,
  labels = c("a", "b"),
  widths = c(1, 1),
  font.label = list(size = 11, color = "black", face = "plain")
)
plot(p12)

# check the outliers
ggsave(
  filename = here("figures/fig_6.pdf"), plot = p12,
  units = "cm", width = 20, height = 10, dpi = 300
)

cor.test(log10(test2.output$mass_mg), log10(test2.output$dry_biomass_mg))

# what about the error?
mean(test2.output$error_perc)
sd(test2.output$error_perc)
round(range(test2.output$error_perc), 3)

test2.output %>%
  filter(error_perc < 200) %>%
  summarise(
    m = mean(error_perc),
    sd = sd(error_perc)
  )
