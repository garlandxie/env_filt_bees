# RLQ Analysis - Xie et al. 2018 
# Script is used to clean data for multiple databases before conducting an RLQ analysis
# Code developed by Garland Xie

# Load libraries ---------------------------------------------------------------
library(dplyr)
library(ade4)
library(readxl)
library(here)
library(janitor)
library(stringr)
library(readr)

# Import files -----------------------------------------------------------------

# community data matrix 
comm <- read_csv(here("data/original", "community_data_matrix.csv"))

# land-use metrics for each site (250m spatial scale)
met_250 <- read_csv(here("data/original", "land_use_metrics_250.csv"))

# land-use metrics for each site (500m spatial scale)
met_500 <- read_csv(here("data/original", "land_use_metrics_500.csv"))

# trait matrix of 31 species (excluded kleptoparasites)
trait <- read_csv(here("data/final", "trait_matrix.csv"))

# data - cleaning: comm matrix -------------------------------------------------

# Change taxon names to abbreviated versions 
# plotting purposes
taxon_change <- function (spList) {
  
  f <- spList %>%
    # separate genus name and species name onto a list
    strsplit(split = "_") %>%
    # create abbreviations of each species 
    # take first four letter of genus (e.g. Anthidium -> Anth)
    # take first three letters of species (e.g. manicatum -> man)
    lapply(function(x) paste(substr(x[1], 1, 3), 
                             substr(x[2], 1, 3), 
                             sep = "_"))
  # return vector of abbreviated taxon names
  return(f)
}  

# clean
comm_tidy <- comm %>%
  
  # alphabetical order
  select(order(colnames(.))) %>%
  
  # concise species code
  select_all(taxon_change) %>%
  rename(ID = X1_NA,
         Anth_man = Ant_man) %>%
  
  # remove kleptoparasites
  select(-c("Ste_ver", "Coe_say", "Coe_alt", "Coe_moe"))

# data - cleaning: trait matrix ------------------------------------------------

trait_clean <- trait %>%
  select(-"X1") %>%
  mutate(ID = taxon_change(spp) %>% 
           str_replace("Ant_man", "Anth_man")) %>%
  
  # change data types
  mutate(native_y_n     = factor(native_y_n),
         emer_time      = factor(emer_time),
         leaf_hair      = factor(leaf_hair),
         leaf_cut       = factor(leaf_cut),
         leaf_pulp      = factor(leaf_pulp),
         resin          = factor(resin),
         stone          = factor(stone),
         none           = factor(none),
         mud            = factor(mud),
         num_nest_mat   = factor(num_nest_mat),
         diet           = factor(diet),
         volt           = factor(volt)
         )

# double-check
str(trait_clean)

# data cleaning - environmental variables -------------------------------------------

# Re-order both environmental data matrices based on rows 
# Makes it easy to do value-matching 
met_250_clean <- met_250 %>%
  rename(ID = X1) %>%
  arrange(ID)

met_500_clean <- met_500 %>%
  rename(ID = X1) %>%
  arrange(ID)

# RLQ: 250 ---------------------------------------------------------------------

# correspondance analysis 
# community data matrix
dudiL_250 <- comm_tidy %>%
  inner_join(met_250_clean, by = "ID") %>%
  select(colnames(comm_tidy), -"ID") %>%
  dudi.coa(scannf = F)

# principal component analysis
# environmental variables
dudiR_250 <- met_250_clean %>%
  select("250_urban" = "total.area_250_urban",
         "250_grass" = "total.area_250_grass",
         "250_tree"  = "total.area_250_tree_canopy") %>%
  dudi.pca(row.w = dudiL_250$lw, scannf = F)

# hill-smith

dudiQ_250 <- trait_clean %>%
  select(-c(ID, spp)) %>%
  dudi.hillsmith(row.w = dudiL_250$cw, 
                 scannf = F)

# RLQ
RLQ_250 <- rlq(dudiR = dudiR_250, 
           dudiQ = dudiQ_250, 
           dudiL = dudiL_250)


