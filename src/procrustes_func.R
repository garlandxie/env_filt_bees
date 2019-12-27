# Calculate partial procruses analysis for functional beta diversity -----------
# Authors(s): Nicholas Sookhan, Garland Xie
# Institutional affiliation: University of Toronto

# libraries --------------------------------------------------------------------
library(vegan)     # for conducting partial procrustes analysis
library(here)      # for creating relative file-paths
library(readr)     # for reading csv files
library(dplyr)     # for manipulating data in R
library(tibble)    # for converting rownames to columnnames ( + vice versa)
library(ggplot2)   # for visualizing data
library(broom)     
library(picante)

# import -----------------------------------------------------------------------

# relative file-paths
env_250_path <- here("data/original", "land_use_metrics_250.csv")
env_500_path <- here("data/original", "land_use_metrics_500.csv")
func_beta_path <- here("data/working", "func_beta_matrices.rds")
dist_spa_path <- here("data/original", "spatial_distance_matrix.csv")

# import the data
env_250 <- read_csv(env_250_path)
env_500 <- read_csv(env_500_path)
func_beta <- readRDS(func_beta_path)
dist_spa <- read_csv(dist_spa_path)

# check packaging --------------------------------------------------------------

# environmental variables
glimpse(env_250)
head(env_250, n = 5)
tail(env_250, n = 5)

glimpse(env_500)
head(env_500, n = 5)
tail(env_500, n = 5)

# spatial distance
glimpse(dist_spa)
head(dist_spa, n = 5)
tail(dist_spa, n = 5)

# beta diversity (total)
dim(as.matrix(func_beta$funct.beta.sor))

# beta diversity (turnover)
dim(as.matrix(func_beta$funct.beta.sim))

# beta diversity (nestedness)
dim(as.matrix(func_beta$funct.beta.sne))

# clean data: functional beta diversity ----------------------------------------

# spatial turnover
func_tu <- as.matrix(func_beta$funct.beta.sim)

# nestedness
func_ne <- as.matrix(func_beta$funct.beta.sne)

# total beta diversity 
func_tot <- as.matrix(func_beta$funct.beta.sor)

# verify equal dimension sizes (rows and columns)
all(dim(func_ne), dim(func_tot), dim(func_tu))

# verify equal rownames
all(rownames(func_ne) == rownames(func_tot))
all(rownames(func_tu) == rownames(func_tot))
all(rownames(func_ne) == rownames(func_tu))

# verify symmetrical matrices
isSymmetric.matrix(func_ne)
isSymmetric.matrix(func_tu)
isSymmetric.matrix(func_tot)

# verify phylogenetic distances are equal or above zero
all(func_tu >= 0)
all(func_ne >= 0)
all(func_tot >= 0)

# clean data: environmental data -----------------------------------------------

# subset metric data to include sites with diversity data
# func beta diversity datasetes are restricted to sites with 2+ species 
# equal dimensions for func beta diversty data: use any matrix for subsetting 

env_250_tidy <- env_250  %>%
   select(ID = "X1",
          grass_250_percent = "prop.landscape_250_grass",
          tree_250_percent  = "prop.landscape_250_tree_canopy",
          urban_250_percent = "prop.landscape_250_urban") %>%
   filter(ID %in% rownames(func_tu)) %>%
   column_to_rownames(var = "ID")

env_500_tidy <- env_500 %>%
   select(ID = "X1",
          grass_500_percent = "prop.landscape_500_grass",
          tree_500_percent  = "prop.landscape_500_tree_canopy",
          urban_500_percent = "prop.landscape_500_urban") %>%
   filter(ID %in% rownames(func_tu)) %>%
   column_to_rownames(var = "ID")

# double check
glimpse(env_250_tidy)
glimpse(env_500_tidy)

# double check: subsets have an environmental gradient?
# high variation in % impervious cover
# use histograms for quick data viz 

env_250_tidy %>%
   ggplot(aes(x = urban_250_percent)) + 
   geom_histogram(bins = 30) + 
   labs(x = "% Impervious Cover", 
        y = "", 
        title = "250m buffer") + 
   theme_minimal() 

env_500_tidy %>%
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
          rownames(env_250_tidy)) %>% 
   filter(ID %in% rownames(env_250_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

dist_spa_500 <- dist_spa %>%
   select(ID = X1, 
          rownames(env_500_tidy)) %>% 
   filter(ID %in% rownames(env_500_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()
   
# double checks
isSymmetric.matrix(dist_spa_250)
isSymmetric.matrix(dist_spa_500)


# PCA: environmental data (250m) -----------------------------------------------

# perform PCA on environmental variables
# standardize all variables
pc_env_250 <- prcomp(env_250_tidy, scale = TRUE) 

# check for cumulative eigenvalues
summary(pc_env_250)
screeplot(pc_env_250)

# plot score loadings
biplot(pc_env_250)

# get scores for the first two PCA axes
scores_env_250 <- scores(pc_env_250, display = "sites", choice = 1:2)

# PCA: environmental distance (500m) ------------------------------------------

# perform PCA on environmental variables
# standardize all variables
pc_env_500 <- prcomp(env_500_tidy, scale = TRUE) 

# check for cumulative eigenvalues
summary(pc_env_500)
screeplot(pc_env_500)

# plot score loadings
biplot(pc_env_500)

# get scores for first two PCA axes
scores_env_500 <- scores(pc_env_500, display = "sites", choice = 1:2)

# PCoA: spatial distance (250m) ------------------------------------------------

# conduct principal coordinate analysis
# apply square root transformation on raw data to avoid negative eigenvalues
# Taken from Legendre (2018). Principal Coordinate Analysis. 
pc_spa_250 <- pcoa(sqrt(dist_spa_250))

# double check: no negative eigenvalues
plot(pc_spa_250$values$Eigenvalues)

# plot 
pc_spa_250$vectors %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores for both axes
scores_spa_250 <- pc_spa_250$vectors[, 1:2]

# PCoA: spatial distance (500m) ------------------------------------------------

# conduct principal coordinate analysis
# apply square root transformation on raw data to avoid negative eigenvalues
# Taken from Legendre (2018). Principal Coordinate Analysis. 
pc_spa_500 <- pcoa(sqrt(dist_spa_500))

# double check: no negative eigenvalues
plot(pc_spa_500$values$Eigenvalues)

# plot 
pc_spa_500$vectors %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores for both axes
scores_spa_500 <- pc_spa_500$vectors[, 1:2]

# PCoA: func nestedness (250m) --------------------------------------------

# prep nestedness
func_ne_250 <- func_ne %>%
   as.data.frame() %>%
   select(rownames(env_250_tidy)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(env_250_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# verify dimensions
dim(func_ne_250)

# pcoa on nestedness
pc_func_ne_250 <- pcoa(func_ne_250, correction = "lingoes")

# check for negative eigenvalues
all(pc_func_ne_250$values$Eigenvalues >= 0)

# check for cumulative eigenvalues
pc_func_ne_250$values$Rel_corr_eig

# plot
pc_func_ne_250$vectors %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_func_ne_250 <- pc_func_ne_250$vectors.cor[, 1:2]

# PCoA: func nestedness (500m) --------------------------------------------

# prep nestedness
func_ne_500 <- func_ne %>%
   as.data.frame() %>%
   select(rownames(env_500_tidy)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(env_500_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# pcoa on nestedness
pc_func_ne_500 <- pcoa(func_ne_500, correction = "lingoes")

# check for negative eigenvalues
all(pc_func_ne_500$values$Eigenvalues >= 0)

# check for cumulative eigenvalues
pc_func_ne_500$values$Cum_corr_eig

# plot
pc_func_ne_500$vectors %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_func_ne_500 <- pc_func_ne_500$vectors.cor[, 1:2]

# PCoA: func turnover (250m) --------------------------------------------

func_tu_250 <- func_tu %>%
   as.data.frame() %>%
   select(rownames(env_250_tidy)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(env_250_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# verify dimensions
dim(func_tu_250)

# pcoa on nestedness
pc_func_tu_250 <- pcoa(func_tu_250, correction = "lingoes")

# check for negative eigenvalues
all(pc_func_tu_250$values$Corr_eig >= 0)

# check for cumulative eigenvalues
pc_func_tu_250$values$Cum_corr_eig

pc_func_tu_250$vectors.cor %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_func_tu_250 <- pc_func_tu_250$vectors.cor[, 1:2]

# PCoA: func turnover (500m) --------------------------------------------

func_tu_500 <- func_tu %>%
   as.data.frame() %>%
   select(rownames(env_500_tidy)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(env_500_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# verify dimensions
dim(func_tu_500)

# pcoa on turnover (500m scale)
# apply a lingoez transformation to avoid negative eigenvalues
pc_func_tu_500 <- pcoa(func_tu_500, correction = "lingoes")

# check for negative eigenvalues
pc_func_tu_500$values$Corr_eig

# check for cumulative (lingoes-corrected) eigenvalues
pc_func_tu_500$vectors.cor

# quick plot
pc_func_tu_500$vectors.cor %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_func_tu_500 <- pc_func_tu_500$vectors.cor[, 1:2]

# PCoA: func total (250m) --------------------------------------------

func_tot_250 <- func_tot %>%
   as.data.frame() %>%
   select(rownames(env_250_tidy)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(env_250_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# verify dimensions
dim(func_tot_250)

# pcoa on total functional beta diversity 
# apply lingoes transformation to correct for negative eigenvalues
pc_func_tot_250 <- pcoa(func_tot_250, correction = "lingoes")

# double check lingoes-corrected eigenvalues
pc_func_tot_250$values$Corr_eig

# check for cumulative (lingoes-corrected) eigenvalues
pc_func_tot_250$values$Cum_corr_eig

# quick plot
pc_func_tot_250$vectors.cor %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_func_tot_250 <- pc_func_tot_250$vectors.cor[, 1:2]

# PCoA: func total (500m) --------------------------------------------

func_tot_500 <- func_tot %>%
   as.data.frame() %>%
   select(rownames(env_500_tidy)) %>%
   rownames_to_column(var = "ID") %>%
   filter(ID %in% rownames(env_500_tidy)) %>%
   column_to_rownames(var = "ID") %>%
   as.matrix()

# verify dimensions
dim(func_tot_500)

# pcoa on total functional beta diversity 
# apply lingoes transformation to correct for negative eigenvalues
pc_func_tot_500 <- pcoa(func_tot_500, correction = "lingoes")

# double check lingoes-corrected eigenvalues
pc_func_tot_500$values$Corr_eig

# check for cumulative (lingoes-corrected) eigenvalues
pc_func_tot_250$values$Cum_corr_eig

# quick plot
pc_func_tot_500$vectors.cor %>%
   as.data.frame() %>%
   select(Axis.1, Axis.2) %>%
   ggplot(aes(x = Axis.1, y = Axis.2)) + 
   geom_point() + 
   theme_minimal()

# get scores of all axes
scores_func_tot_500 <- pc_func_tot_500$vectors.cor[, 1:2]

# Partial procrustes analysis: nestedness --------------------------------------

# control for space (250m buffer)
resid_spa_phy_ne_250<- resid(lm(scores_func_ne_250 ~ scores_spa_250))
resid_spa_env_250 <- resid(lm(scores_env_250 ~ scores_spa_250))

# control for space (500m buffer)
resid_spa_func_ne_500 <- resid(lm(scores_func_ne_500 ~ scores_spa_500))
resid_spa_env_500 <- resid(lm(scores_env_500 ~ scores_spa_500))

# run analysis 
parpro_ne_250 <- protest(X = resid_spa_env_250, Y = resid_spa_phy_ne_250)
parpro_ne_500 <- protest(X = resid_spa_env_500, Y = resid_spa_func_ne_500)

# Partial procustes analysis: turnover -----------------------------------------

# control for space (250m buffer)
resid_spa_func_tu_250<- resid(lm(scores_func_tu_250 ~ scores_spa_250))
resid_spa_env_250 <- resid(lm(scores_env_250 ~ scores_spa_250))

# control for space (500m buffer)
resid_spa_func_tu_500 <- resid(lm(scores_func_tu_500 ~ scores_spa_500))
resid_spa_env_500 <- resid(lm(scores_env_500 ~ scores_spa_500))

# run analysis 
parpro_tu_250 <- protest(resid_spa_env_250, resid_spa_func_tu_250)
parpro_tu_500 <- protest(resid_spa_env_500, resid_spa_func_tu_500)

# Partial procustes analysis: total --------------------------------------------

# control for space (250m buffer)
resid_spa_func_tot_250<- resid(lm(scores_func_tot_250 ~ scores_spa_250))
resid_spa_env_250 <- resid(lm(scores_env_250 ~ scores_spa_250))

# control for space (500m buffer)
resid_spa_func_tot_500 <- resid(lm(scores_func_tot_500 ~ scores_spa_500))
resid_spa_env_500 <- resid(lm(scores_env_500 ~ scores_spa_500))

# run analysis 
parpro_tot_250 <- protest(resid_spa_env_250, resid_spa_func_tot_250)
parpro_tot_500 <- protest(resid_spa_env_500, resid_spa_func_tot_500)
