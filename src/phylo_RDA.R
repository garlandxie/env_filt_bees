# libraries -------------------------------------------------------------------
library(here)
library(vegan)
library(dplyr)
library(readr)
library(tibble)
library(faraway)  # for multicollinearity
library(ape)
library(phytools)

# import ----------------------------------------------------------------------

# file-paths
env_250_path <- here("data/original", "land_use_metrics_250.csv")
env_500_path <- here("data/original", "land_use_metrics_500.csv")
comm3A_path <- here("data/original", "community_data_matrix.csv")
tree_path <- here("data/final", "tree_rescale_ulm.new")

# read data
env_250_df <- read_csv(env_250_path)
env_500_df <- read_csv(env_500_path)
comm3A_df <- read_csv(comm3A_path)
tree_ulm <- read.tree(tree_path)

# check packaging --------------------------------------------------------------

str(env_250_df)
head(env_250_df, n = 5)
tail(env_250_df, n = 5)

str(env_500_df)
head(env_500_df, n = 5)
tail(env_500_df, n = 5)


# clean data: env 250 ----------------------------------------------------------

env_250_tidy <- env_250_df %>%
  
  # pick composite measures (e.g., %)
  select("ID" = X1,
         "percent_urban_250" = prop.landscape_250_urban,
         "percent_grass_250" = prop.landscape_250_grass,
         "percent_tree_250"  = prop.landscape_250_tree_canopy) %>%
  
  # subset sites by the community data matrix
  filter(ID %in% rownames(comm3A_df)) 
  


# data clean: env 500 ----------------------------------------------------------

env_500_tidy <- env_500_df %>%
  
  # pick composite measures (e.g., %)
  select("ID" = X1,
         "percent_urban_500" = prop.landscape_500_urban,
         "percent_grass_500" = prop.landscape_500_grass,
         "percent_tree_500"  = prop.landscape_500_tree_canopy) %>%
  
  # subset sites by the community data matrix
  filter(ID %in% rownames(comm3A_df)) 

# RDA: 250m --------------------------------------------------------------------

# check for co-linear environment variables
# % urban and % tree do covary
pairs(env_250_tidy[, -1])

# check for multicollinearity
# VIF values > 10 ==> should drop a variable
vif(env_250_tidy[, -1])

# 
comm_3A_tidy <- comm3A_df %>%
  rownames_to_column(var = "ID") %>%
  filter(ID %in% env_250_tidy$ID) %>%
  column_to_rownames(var = "ID") %>%
  mutate_if(., is.integer, as.numeric) %>%
  select(-c("Stelis_verticalis", 
            "Coelioxys_sayi", 
            "Coelioxys_moesta", 
            "Coelioxys_alternata")) %>%
  decostand(method = "hellinger")
  
# remove % tree from df
env_250_tidy2 <- env_250_tidy[, c("percent_urban_250", "percent_grass_250")]

# perform RDA
rda_250 <- rda(formula = comm_3A_tidy ~ percent_urban_250 + percent_grass_250, 
               data = env_250_tidy2)

# adjusted R-2
(R2adj <- RsquareAdj(rda_250)$adj.r.squared)

# global permutation test
anova(rda_250, permutations = how(nperm = 999))

# permutation test for each canonical axis
anova(rda_250, by = "axis", permutations = how(nperm = 999))

# get species scores that are constrained by env vars
(spp_scores <- rda_250$CCA$v)

# perform Blomberg's K test on species score of RDA axis 1
phylosig(tree = tree_ulm, 
         x = spp_scores[,1],
         method = "K",
         test = TRUE, 
         nsim = 999)

# RDA: 500 m -------------------------------------------------------------------

# check for co-linear environment variables
# % urban and % tree do covary
pairs(env_500_tidy[, -1])

# check for multicollinearity
# VIF values > 10 ==> should drop a variable
vif(env_500_tidy[, -1])

# 
comm_3A_tidy <- comm3A_df %>%
  rownames_to_column(var = "ID") %>%
  filter(ID %in% env_500_tidy$ID) %>%
  column_to_rownames(var = "ID") %>%
  mutate_if(., is.integer, as.numeric) %>%
  select(-c("Stelis_verticalis", 
            "Coelioxys_sayi", 
            "Coelioxys_moesta", 
            "Coelioxys_alternata")) %>%
  decostand(method = "hellinger")

# remove % tree from df
env_500_tidy2 <- env_500_tidy[, c("percent_urban_500", "percent_grass_500")]

# perform RDA
rda_500 <- rda(formula = comm_3A_tidy ~ percent_urban_500 + percent_grass_500, 
               data = env_500_tidy2)

# adjusted R-2
(R2adj <- RsquareAdj(rda_500)$adj.r.squared)

# global permutation test
anova(rda_500, permutations = how(nperm = 999))

# permutation test for each canonical axis
anova(rda_500, by = "axis", permutations = how(nperm = 999))

# get species scores that are constrained by env vars
(spp_scores_500 <- rda_500$CCA$v)

# perform Blomberg's K test on species score of RDA axis 1
phylosig(tree = tree_ulm, 
         x = spp_scores_500[,1],
         method = "K",
         test = TRUE, 
         nsim = 999)