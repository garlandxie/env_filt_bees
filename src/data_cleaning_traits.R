# libraries --------------------------------------------------------------------
library(here) 
library(dplyr)
library(janitor)
library(stringr)
library(validate)
library(forcats)
library(readr)
library(stringr)

# import -----------------------------------------------------------------------
trait <- read_csv(here("data/original", "bee_trait_matrix.csv"))

# check packaging --------------------------------------------------------------

# data types, headers and footers
str(trait)
head(trait, n = 5)
tail(trait, n = 5)

# data validation
val <- trait %>%
  clean_names() %>%
  as.data.frame() %>%
  check_that(
    native_y_n %in% c("Y", "N"),
    emer_time_a > 0,
    !is.na(leaf_hair), 
    !is.na(leaf_cut), 
    !is.na(leaf_pulp), 
    !is.na(resin), 
    !is.na(mud), 
    !is.na(stone), 
    !is.na(none), 
    diet %in% c("Poly", "Oligo"), 
    volt_a %in% c(1, 2), 
    itd > 0) %>%
  summary()
  
# clean data -------------------------------------------------------------------

# convert nesting material into a binary variable (0/1)
nest_mat <- c("leaf_hair", 
              "leaf_cut", 
              "leaf_pulp", 
              "resin", 
              "mud", 
              "stone", 
              "none")  

# clean
trait_tidy <- trait %>%
  clean_names() %>%
  mutate(
    
    # fix some spelling mistakes 
    species = str_replace(species, pattern = "centuncularis", replacement = "centucularis") %>%
              str_replace(pattern = "atriventris",  replacement = "atriventis") %>%
              str_replace(pattern = "carinata", replacement = "crucifera"), 
  
    # change emergence time time levels to 
    # Spring: (Days 0 - 14)
    # Summer: (Days 14+ )
    emer_time2 = ifelse(emer_time <= 2, "Spring", "Summer") %>%
                 factor(),
    
    # change voltinism into a factor
    volt = factor(volt), 
    diet = factor(diet), 
  
    # create ID column to match with phylogenetic tree
    spp = paste(genus, species, sep = "_")) %>%
  
  # relabel the nesting material data into "Y" and "N"
  mutate_at(vars(nest_mat), ~replace(., is.na(.), "0") %>%
              str_replace("X", "1") %>%
              factor()) %>% 
  
  # change coding for native and exotic status
  # native status == 1
  # exotic status == 0 
  mutate(native_y_n = str_replace(native_y_n, "Y", "1") %>%
                      str_replace("N", "0") %>%
                      factor()
         ) %>%
  
# select relevant trait data (see below)
  select(spp, native_y_n, emer_time2, leaf_hair,
         leaf_cut, leaf_pulp, resin, mud, 
         stone, none, diet, volt, itd) 

# double-check
glimpse(trait_tidy)

# Save to disk -----------------------------------------------------------------

# comma-delimited file, for data repos
write.csv(trait_tidy, file = here("data/final", "trait_matrix.csv"))

# preserve data-types, save as RDS
saveRDS(trait_tidy, file = here("data/final", "trait_matrix.rds"))

