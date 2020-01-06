# Partial Procrustes Analysis --------------------------------------------------

# Goal of script: to determine a relationship between:
# phylogenetic beta diversity (total, nest, turnover) and 
# environmental distance

# Author(s): Nicholas Sookhan, Garland Xie
# Institutional affiliation: University of Toronto

# libraries --------------------------------------------------------------------
library(vegan)   # for calculating partial procrustes analysis
library(here)    # for creating relative file-paths
library(picante) # for calculating principal coordinate analysis
library(tidyverse)

# import -----------------------------------------------------------------------

# relative file-paths
env_250_path <- here("data/original", "land_use_metrics_250.csv")
env_500_path <- here("data/original", "land_use_metrics_500.csv")
phylo_path <- here("data/working", "phylo_beta_matrices.rds")
dist_spa_path <- here("data/original", "spatial_distance_matrix.csv")

# load data
l_met_250 <- read_csv(env_250_path)
l_met_500 <- read_csv(env_500_path)
phylo_beta <- readRDS(phylo_path)
dist_spa <- read_csv(dist_spa_path)

# check packaging --------------------------------------------------------------

# environmental data
str(l_met_250)
head(l_met_250, n = 5)
tail(l_met_250, n = 5)

str(l_met_500)
head(l_met_500, n = 5)
tail(l_met_500, n = 5)

# spatial distance matrix
str(dist_spa)
head(dist_spa, head = 5)
tail(dist_spa, head = 5)

# phylo beta diversity matrix
str(as.matrix(phylo_beta$phylo.beta.sor)) 
head(as.matrix(phylo_beta$phylo.beta.sor))
tail(as.matrix(phylo_beta$phylo.beta.sor))

# phylo beta diversity matrix (spatial turn-over)
str(as.matrix(phylo_beta$phylo.beta.sim)) 
head(as.matrix(phylo_beta$phylo.beta.sim))
tail(as.matrix(phylo_beta$phylo.beta.sim))

# phylo beta diversity matrix (nestedness)
str(as.matrix(phylo_beta$phylo.beta.sne)) 
head(as.matrix(phylo_beta$phylo.beta.sne))
tail(as.matrix(phylo_beta$phylo.beta.sne))

# clean data: phylo beta diversity ---------------------------------------------

# spatial turnover
phylo_tu <- as.matrix(phylo_beta$phylo.beta.sim)

# nestedness
phylo_ne <- as.matrix(phylo_beta$phylo.beta.sne)

# total beta diversity 
phylo_tot <- as.matrix(phylo_beta$phylo.beta.sor)

# verify equal dimension sizes (rows and columns)
all(dim(phylo_ne), dim(phylo_tot), dim(phylo_tu))

# verify equal rownames
all(rownames(phylo_ne) == rownames(phylo_tot))
all(rownames(phylo_tu) == rownames(phylo_tot))
all(rownames(phylo_ne) == rownames(phylo_tu))

# verify symmetrical matrices
isSymmetric.matrix(phylo_ne)
isSymmetric.matrix(phylo_tu)
isSymmetric.matrix(phylo_tot)

# verify phylogenetic distances are equal or above zero
all(phylo_tu >= 0)
all(phylo_ne >= 0)
all(phylo_tot >= 0)

# clean data: environmental data -----------------------------------------------

# subset metric data to include sites with diversity data
met_250_clean <- l_met_250 %>%
   select(ID = "X1",
          grass_250_percent = "prop.landscape_250_grass",
          tree_250_percent  = "prop.landscape_250_tree_canopy",
          urban_250_percent = "prop.landscape_250_urban") %>%
   filter(ID %in% rownames(phylo_tu)) %>%
   column_to_rownames(var = "ID")

met_500_clean <- l_met_500 %>%
   select(ID = "X1",
          grass_500_percent = "prop.landscape_500_grass",
          tree_500_percent  = "prop.landscape_500_tree_canopy",
          urban_500_percent = "prop.landscape_500_urban") %>%
   filter(ID %in% rownames(phylo_tu)) %>%
   column_to_rownames(var = "ID")

# double check: subsets have an environmental gradient?
# high variation in % impervious cover
# use histograms for quick data viz 

met_250_clean %>%
   ggplot(aes(x = urban_250_percent)) + 
   geom_histogram(bins = 30) + 
   labs(x = "% Impervious Cover", 
        y = "", 
        title = "250m buffer") + 
   theme_minimal() 

met_500_clean %>%
   ggplot(aes(x = urban_500_percent)) + 
   geom_histogram(bins = 30) + 
   labs(x = "% Impervious Cover", 
        y = "", 
        title = "500m buffer") + 
   theme_minimal() 

# clean data: spatial distance  ------------------------------------------------

# subset to match sites in environmental distance matrix
dist_spa_250 <- dist_spa %>%
   select(ID = X1, 
          rownames(met_250_clean)) %>% 
   filter(ID %in% rownames(met_250_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()
 
dist_spa_500 <- dist_spa %>%
   select(ID = X1, 
          rownames(met_500_clean)) %>% 
   filter(ID %in% rownames(met_500_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# double checks
isSymmetric.matrix(dist_spa_250)
isSymmetric.matrix(dist_spa_500)

# PCA: environmental data (250m) -----------------------------------------------

# perform PCA on environmental variables
# standardize all variables
pc_env_250 <- prcomp(met_250_clean, scale = TRUE) 

# check for cumulative eigenvalues
summary(pc_env_250)
screeplot(pc_env_250)

# get scores for all three axes
scores_env_250 <- scores(pc_env_250, display = "sites", choice = 1:2)

# PCA - environmental distance (500m) ------------------------------------------

# perform PCA on environmental variables
# standardize all variables
pc_env_500 <- prcomp(met_500_clean, scale = TRUE) 

# check for cumulative eigenvalues
summary(pc_env_500)
screeplot(pc_env_500)

# get scores for both axes
scores_env_500 <- scores(pc_env_500, display = "sites", choice = 1:2)

# PCoA: spatial distance (250m) ------------------------------------------------

# conduct principal coordinate analysis
pc_spa_250 <- pcoa(sqrt(dist_spa_250))

# plot 
pc_spa_250$vectors %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point()

# get scores for both axes
scores_spa_250 <- pc_spa_250$vectors[, 1:2]

# PCoA: spatial distance (500m) ------------------------------------------------

# conduct principal coordinate analysis
pc_spa_500 <- pcoa(sqrt(dist_spa_500))

# plot 
pc_spa_500$vectors %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores for both axes
scores_spa_500 <- pc_spa_500$vectors[, 1:2]

# PCoA: phylo nestedness (250m) --------------------------------------------

# prep nestedness
phylo_ne_250 <- phylo_ne %>%
   as.data.frame() %>%
   select(rownames(met_250_clean)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(met_250_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# pcoa on nestedness
pc_phy_ne_250 <- pcoa(phylo_ne_250, correction = "lingoes")

# plot
pc_phy_ne_250$vectors.cor %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_phy_ne_250 <- pc_phy_ne_250$vectors.cor[, 1:2]

# PCoA: phylo nestedness (500m) --------------------------------------------

# prep nestedness
phylo_ne_500 <- phylo_ne %>%
   as.data.frame() %>%
   select(rownames(met_500_clean)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(met_500_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# pcoa on nestedness
pc_phy_ne_500 <- pcoa(phylo_ne_500, correction = "lingoes")

pc_phy_ne_500$vectors.cor %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_phy_ne_500 <- pc_phy_ne_500$vectors[, 1:2]

# PCoA: phylo turnover (250m) --------------------------------------------

phylo_tu_250 <- phylo_tu %>%
   as.data.frame() %>%
   select(rownames(met_250_clean)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(met_250_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# pcoa on nestedness
pc_phy_tu_250 <- pcoa(phylo_tu_250, correction = "lingoes")

# plot
pc_phy_tu_250$vectors.cor %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_phy_tu_250 <- pc_phy_tu_250$vectors.cor[, 1:2]

# PCoA: phylo turnover (500m) --------------------------------------------

phylo_tu_500 <- phylo_tu %>%
   as.data.frame() %>%
   select(rownames(met_500_clean)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(met_500_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# pcoa 
pc_phy_tu_500 <- pcoa(phylo_tu_500, correction = "lingoes")

# plot
pc_phy_tu_500$vectors.cor %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_phy_tu_500 <- pc_phy_tu_500$vectors.cor[, 1:2]

# PCoA: phylo total (250m) --------------------------------------------

phylo_tot_250 <- phylo_tot %>%
   as.data.frame() %>%
   select(rownames(met_250_clean)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(met_250_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# pcoa on nestedness
pc_phy_tot_250 <- pcoa(phylo_tot_250, correction = "lingoes")

# plot
pc_phy_tot_250$vectors %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_phy_tot_250 <- pc_phy_tot_250$vectors.cor[, 1:2]

# PCoA: phylo total (500m) --------------------------------------------

phylo_tot_500 <- phylo_tot %>%
   as.data.frame() %>%
   select(rownames(met_500_clean)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(met_500_clean)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# pcoa 
pc_phy_tot_500 <- pcoa(phylo_tot_500, correction = "lingoes")

# plot
pc_phy_tot_500$vectors.cor %>%
   as.data.frame() %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_phy_tot_500 <- pc_phy_tot_500$vectors.cor[, 1:2]

# Partial procrustes analysis: nestedness --------------------------------------

# control for space (250m buffer)
resid_spa_phy_ne_250 <- resid(lm(scores_phy_ne_250 ~ scores_spa_250))
resid_spa_env_250 <- resid(lm(scores_env_250 ~ scores_spa_250))

# control for space (500m buffer)
resid_spa_phy_ne_500 <- resid(lm(scores_phy_ne_500 ~ scores_spa_500))
resid_spa_env_500 <- resid(lm(scores_env_500 ~ scores_spa_500))

# run analysis 
parpro_ne_250 <- protest(X = resid_spa_env_250, Y = resid_spa_phy_ne_250)
parpro_ne_500 <- protest(X = resid_spa_env_500, Y = resid_spa_phy_ne_500)

# Partial procustes analysis: turnover -----------------------------------------

# control for space (250m buffer)
resid_spa_phy_tu_250<- resid(lm(scores_phy_tu_250 ~ scores_spa_250))
resid_spa_env_250 <- resid(lm(scores_env_250 ~ scores_spa_250))

# control for space (500m buffer)
resid_spa_phy_tu_500 <- resid(lm(scores_phy_tu_500 ~ scores_spa_500))
resid_spa_env_500 <- resid(lm(scores_env_500 ~ scores_spa_500))

# run analysis 
parpro_tu_250 <- protest(resid_spa_env_250, resid_spa_phy_tu_250)
parpro_tu_500 <- protest(resid_spa_env_500, resid_spa_phy_tu_500)

# Partial procustes analysis: total --------------------------------------------

# control for space (250m buffer)
resid_spa_phy_tot_250<- resid(lm(scores_phy_tot_250 ~ scores_spa_250))
resid_spa_env_250 <- resid(lm(scores_env_250 ~ scores_spa_250))

# control for space (500m buffer)
resid_spa_phy_tot_500 <- resid(lm(scores_phy_tot_500 ~ scores_spa_500))
resid_spa_env_500 <- resid(lm(scores_env_500 ~ scores_spa_500))

# run analysis 
parpro_tot_250 <- protest(resid_spa_env_250, resid_spa_phy_tot_250)
parpro_tot_500 <- protest(resid_spa_env_500, resid_spa_phy_tot_500)


   

