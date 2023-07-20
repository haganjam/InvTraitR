
# generate a test database that uses many different types of equations

# load the relevant libraries
library(dplyr)
library(ggplot2)
library(rstan)

# load the use-scripts
source("companion_scripts/helper-plot-theme.R")
source("companion_scripts/03_accuracy_analysis/helper-miscellaneous.R")

# load the required functions
source("R/clean_taxon_names.R")
source("R/get_habitat_data.R")
source("R/select_traits_tax_dist.R")
source("R/special_names.R")
source("R/helpers.R")
source("R/get_trait_from_taxon.R")

# load libraries required for those function
library(igraph)
library(assertthat)

# load the test data
dat <- readRDS("database/test_a_data_compilation.rds")

output <-
  get_trait_from_taxon(
    data = dat,
    target_taxon = "taxon",
    life_stage = "life_stage",
    body_size = "length_mm",
    latitude_dd = "lat",
    longitude_dd = "lon",
    workflow = "workflow1",
    trait = "equation",
    max_tax_dist = 4,
    gen_sp_dist = 0.5
  )

# remove the entries without an equation id
output <- 
  output |>
  dplyr::filter(!is.na(id))

# only keep entries with matching life-stages
output <- 
  output |>
  dplyr::filter(life_stage_match == TRUE)

# set-up a vector to capture the dry biomass values
dry_biomass_mg <- vector(length = nrow(output))

# loop over all the rows
for(i in 1:nrow(output )) {
  
  # get the ith row of data
  L <- unlist(output[i, "length_mm"], use.names = FALSE)
  model <- output[i,]$equation_form
  log_base <- output[i,]$log_base
  a <- output[i,]$a
  b <- output[i,]$b
  CF <- output[i,]$lm_correction
  scale <- output[i,]$dry_biomass_scale
  
  # evalulate the equation
  if( any(is.na(c(a, b))) )  {
    
    dry_biomass_mg[i] <- NA
    
  } else if (model == "model1") {
    
    # calculate the raw prediction on the log-scale
    x <- a + (b*logb(x = L, base = log_base))
    
    # convert to the natural scale
    x <- (log_base^x)
    
    # apply the correction factor
    dry_biomass_mg[i] <- ifelse(!is.na(CF), x*CF, x)*scale
    
  } else if (model == "model2") {
    
    # calculate the raw prediction
    dry_biomass_mg[i] <- a*(L^b)*scale
    
  }
  
}

# add this dry biomass estimate to the data.frame
output[["dry_biomass_mg"]] <- dry_biomass_mg

# create a data.frame for modelling the error
output_df <- 
  output |>
  dplyr::select(row, order, taxon, db_scientificName,
         tax_distance, body_size_range_match, equation_form, r2_match, n,
         realm_match, major_habitat_type_match, ecoregion_match,
         length_mm, obs_dry_biomass_mg, dry_biomass_mg
         )

# get the percentage prediction error
output_df <- 
  output_df |>
  dplyr::mutate(error_perc = ((obs_dry_biomass_mg - dry_biomass_mg)/obs_dry_biomass_mg)*100,
                abs_error_perc = (abs(obs_dry_biomass_mg - dry_biomass_mg)/obs_dry_biomass_mg)*100) 


# get a habitat match variable
output_df[["habitat_match"]] <- 
  apply(output_df[, paste0(c("realm", "major_habitat_type", "ecoregion"), "_match") ], 1,
        function(x) sum(x, na.rm = TRUE) )

# absolute error

# check the distribution of these errors
hist(output_df$abs_error_perc)
summary(output_df$abs_error_perc)

# what about the log-transformed abs error perc
hist(log(output_df$abs_error_perc))

# error

# check the distribution
hist(output_df$error_perc)
summary(output_df$error_perc)

# how many data points do we have per taxon id
output_df |>
  dplyr::group_by(row) |>
  dplyr::summarise( n = n()) |>
  dplyr::ungroup() |>
  dplyr::summarise(mean = mean(n),
                   min = min(n),
                   max = max(n))

# write a function to min-max standardise
min_max <- function(x) {
  (x - min(x))/(max(x)-min(x))
}

# create a list with the relevant data
dat <- list(N = nrow(output_df),
            G = length(unique(output_df$row)),
            id = as.integer(as.factor(output_df$row)),
            td = min_max(output_df$tax_distance),
            bs = as.integer(output_df$body_size_range_match),
            hm = min_max(output_df$habitat_match),
            abs_error = (output_df$abs_error_perc) )
summary(dat)
str(dat)

# how many groups do we have now?
length(unique(dat$id))
length(dat$id)

# compile the stan growth rate model
m1 <- rstan::stan_model("companion_scripts/03_accuracy_analysis/04_model_error_variation.stan",
                        verbose = TRUE)
print(m1)

# sample the stan model: m1
m1_fit <- rstan::sampling(m1, data = dat, 
                          iter = 2500, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99),
                          seed = 54856)

# check the output
print(m1_fit)
m1_fit@model_pars

# extract the diagnostic parameters
diag <- summary(m1_fit)
par <- "abar"
diag$summary[grepl(par, row.names(diag$summary)), ]

# view some traceplots
traceplot(m1_fit, pars = c("abar[1]", "abar[2]", "abar[3]", "abar[4]"))

# extract the posterior distribution
m1_post <- extract(m1_fit)

# check how well the model fits the data
a_id <- m1_post$a_id
b1_td <- m1_post$b1_td
b2_bs <- m1_post$b2_bs
b3_hm <- m1_post$b3_hm
sigma <- m1_post$sigma

N <- length(dat$id)

# write a loop
y_list <- vector("list", length = length(sigma))
for (i in 1:length(sigma)) {
  
  mu <- sapply(1:N, function(x) {
    
    a_id[i,dat$id[x]] + b1_td[i,dat$id[x]] * dat$td[x] + b2_bs[i,dat$id[x]] * dat$bs[x] + b3_hm[i,dat$id[x]] * dat$hm[x]
    
  })
  
  y_list[[i]] <- exp(mu)
  
}

# bind into a matrix
y_df <- do.call("cbind", y_list)
dim(y_df)

# calculate the mean across samples
y <- apply(y_df, 1, mean)
y_PI_low <- apply(y_df, 1, function(x) {PI(x, 0.90)[1]})
y_PI_high <-apply(y_df, 1, function(x) {PI(x, 0.90)[2]}) 

# pull this into a data.frame
df_obs <- data.frame(obs_error = dat$abs_error,
                     est_error = y,
                     PI_low = y_PI_low,
                     PI_high = y_PI_high)

p1 <- 
  ggplot(data = df_obs,
       mapping = aes(x = obs_error, y = est_error)) +
  geom_point(shape = 16, alpha = 0.1, colour = wesanderson::wes_palette(name = "Darjeeling1", n = 1)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ylab("Predicted percentage prediction error (%)") +
  xlab("Observed percentage prediction error (%)") +
  theme_meta()

ggsave(filename = "figures/fig_S2.png", p1, dpi = 400,
       units = "cm", width = 10, height = 10)

# check how many mean values are predicted to be negative
sum(y<0)/length(y)
cor(y, dat$abs_error)

# calculate an r2 value
r2 <- 
  
  apply(y_df, 2, function(x) {
    
    r <- x - dat$abs_error
    r <- 1 - (var2(r)/var2(dat$abs_error))
    
    return(r)
    
  } )

# calculate the summary statistics of the r2 value
mean(r2)
range(r2)

# how to calculate the contrasts:
# https://vasishth.github.io/bayescogsci/book/ch-reg.html#thm:lognormal

# function to calculate the difference in median given values of covariates
median_diff <- function(cov1_name, cov1_val,
                        cov2_name, cov2_val,
                        do_name, do_val) {
  
  # covariate grid
  pred_grid1 <- expand.grid(v1 = cov1_val,
                            v2 = cov2_val)
  names(pred_grid1) <- c(cov1_name, cov2_name)
  
  # add the do_val
  pred_grid1[[do_name]] <- do_val[1]
  
  # add the do_val
  pred_grid2 <- pred_grid1
  pred_grid2[[do_name]] <- do_val[2]
  
  # median absolute error when td = 0
  m_diff_list <- vector("list", length = nrow(pred_grid1))
  for(i in 1:nrow(pred_grid1)) {
    
    m_ae1 <- with(pred_grid1[i,],
                  exp(a_id + b1_td*td1 + b2_bs*bs1 + b3_hm*hm1))
    
    m_ae2 <- with(pred_grid2[i,],
                  exp(a_id + b1_td*td1 + b2_bs*bs1 + b3_hm*hm1))
    
    m_diff <- sapply( (m_ae2 - m_ae1), function(x) x  )
    
    m_diff_list[[i]] <- m_diff
    
  }
  
  return(unlist(m_diff_list))
  
}

# calculate the effect of taxonomic distance
# difference between taxonomic distance of 4 versus 0
x <- median_diff(cov1_name = "bs1", cov1_val = c(0, 1), 
                 cov2_name = "hm1", cov2_val = seq(0, 1, 0.1), 
                 do_name = "td1", do_val = c(0, 1))
xm <- mean(x)
xpi <- PI(x, 0.90)

# calculate the effect of body size matching
# going from unmatching to matching
y <- median_diff(cov1_name = "td1", cov1_val = seq(0, 1, 0.1), 
                 cov2_name = "hm1", cov2_val = seq(0, 1, 0.1), 
                 do_name = "bs1", do_val = c(0, 1))
ym <- mean(y)
ypi <- PI(y, 0.90)

# calculate the effect of habitat matching
# going from 0 habitat match to 
z <- median_diff(cov1_name = "td1", cov1_val = seq(0, 1, 0.1), 
                 cov2_name = "bs1", cov2_val = c(0, 1), 
                 do_name = "hm1", do_val = c(0, 1))
zm <- mean(z)
zpi <- PI(z, 0.90)

# pull these results into a data.frame
coef_plot <- data.frame(effect = c("Taxonomic distance", 
                                   "Body size range match",
                                   "Habitat match"),
                        mean = c(xm, ym, zm),
                        PI_min = c(xpi[1], ypi[1], zpi[1]),
                        PI_max = c(xpi[2], ypi[2], zpi[2]))
coef_plot$effect <- factor(coef_plot$effect,
                           levels = c("Taxonomic distance", "Body size range match", "Habitat match"))

p2 <- 
  ggplot(data = coef_plot) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(mapping = aes(x = effect, y = mean), size = 3) +
  geom_errorbar(mapping = aes(x = effect, ymin = PI_min, ymax = PI_max),
                width = 0) +
  theme_meta() +
  ylab("Percentage prediction error (%)") +
  xlab(NULL) +
  scale_y_continuous(limits = c(-50, 85)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) +
  annotate(geom = "text", label = expression(r^{2}~" = 0.26"),
           x = 0.9, y = 80)

ggsave(filename = "figures/fig_6.png", p2, dpi = 400,
       units = "cm", width = 10, height = 11)
