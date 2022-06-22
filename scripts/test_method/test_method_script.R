
# Test the method for calculating species biomasses

library(here)
library(dplyr)
library(readr)
library(ggplot2)
source(here("scripts/use_database/01_search_equ_len_database_functions.R"))

# load the test data
t.dat <- readxl::read_xlsx("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/allometry_database/test_data.xlsx")
head(t.dat)

# make the life-stage column completely NA's
t.dat$Life_stage <- (rep(NA, nrow(t.dat)))
summary(t.dat)

# try the function
x <- get_taxa_mass(data.base = "itis",
                   max_tax_dist = 6,
                   data = t.dat %>% filter(Reference != "Eklof_2016"),
                   target.name.col = "Taxa",
                   life.stage.col = "Life_stage",
                   length.col = "Length_mm")
View(x)              

# join this data.frame to the actual mass data
y <- 
  full_join(x,
          t.dat %>%
            rename(target_name = Taxa,
                   target_life_stage = Life_stage,
                   size = Length_mm),
          by = c("target_name", "target_life_stage", "size")) %>%
  select(target_name, target_life_stage, size, dist_to_target,
         mass, Dry_weight_mg) %>%
  filter(!is.na(mass))

y$general_dw_Hebert <- exp(-4.814 + log(y$size) )

exp(-4.814 + log(27) )

y <- 
  y %>%
  tidyr::pivot_longer(cols = c("general_dw_Hebert", "mass"),
               names_to = "method",
               values_to = "DW")

ggplot(data = y,
       mapping = aes(x = log(Dry_weight_mg), y =log(DW), colour = method )) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_classic()
 
ggplot() +
  geom_point(data = y,
             mapping = aes(x = Dry_weight_mg, y =DW, colour = method )) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_classic()

ggplot() +
  geom_point(data = y %>% filter(method == "general_dw_Hebert"),
             mapping = aes(x = log(Dry_weight_mg), y = log(DW), colour = method )) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_classic()


# test the data by getting biomass conversion data for the Inselberg and Korranneberg data

ins <- read_csv(file = "C:/Users/james/Documents/github/predicting_trait_responses/data/biomass_conversions/aus_ins_bio.csv")

# use the method to get biomass data using default length data
ins$length_dat <- NA

ins.x <- 
  get_taxa_mass(data.base = "itis",
                max_tax_dist = 8,
                data = ins,
                target.name.col = "taxon",
                life.stage.col = "life_stage",
                length.col = "length_dat"
                )

View(ins.x)

ins.y <- 
  get_taxa_info(data.base = "itis",
                max_tax_dist = 6,
                target.name = ins$taxon,
                life.stage = ins$life_stage)

View(ins.y$equation_info)
View(ins.y$length_info)






