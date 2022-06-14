
# Test the method for calculating species biomasses

library(here)
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
                   data = t.dat,
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

library(ggplot2)

ggplot(data = y,
       mapping = aes(x = log(mass), log(Dry_weight_mg) )) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_classic()
 
ggplot(data = y %>% filter(mass < 2),
       mapping = aes(x = (mass), (Dry_weight_mg) )) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_classic()

cor(y$mass, y$Dry_weight_mg)              
cor(log(y$mass), log(y$Dry_weight_mg))              

y %>%
  filter(mass < 10, mass < 30, Dry_weight_mg > 10)

View(x %>% filter(target_name == "Gammarus"))


              
              