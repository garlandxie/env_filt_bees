# calculate functional beta diversity ------------------------------------------
# authors(s): Nicholas Sookhan, Garland Xie

# libraries --------------------------------------------------------------------
library(here)      # for creating relative file-paths
library(vegan)     # for analysing community matrices
library(betapart)  # for calculating turnover, total + nestedness matrices
library(dplyr)     # for manipulating data
library(tibble)    # for converting rownames to columns (+ vice verse)
library(cluster)   # for calculating gower's distance
library(ape)       # for running principal coordinate analysis

# import -----------------------------------------------------------------------

# relative file-paths
comm_path <- here("data/original", "community_data_matrix.csv")
trait_path <- here("data/final", "trait_matrix.rds")

# load data
comm <- read.csv(comm_path, row.names = 1)
trait <- readRDS(trait_path)

# check packaging --------------------------------------------------------------

# community data matrix
str(comm)
head(comm, n = 5)
tail(comm, n = 5)

# trait matrix
str(trait)
head(trait, n = 5)
tail(trait, n = 5)

# data cleaning: traits --------------------------------------------------------

# principal coordinate analysis
trait_tidy <- trait %>%
  column_to_rownames(var = "spp") 

# calculate gower's distance matrix
trait_dist <- as.matrix(daisy(trait_tidy, metric="gower", stand = T))

# avoid negative eigenvalues
trait_pcoa <- pcoa(trait_dist, correction = "cailliez")

# data cleaning: community data matrix -----------------------------------------

# to relative abundance
comm_rel <- decostand(comm, method = "pa") %>%
  
  # 
  rownames_to_column(var = "site") %>%
  
  # remove kleptoparasites
  select(site, rownames(trait_tidy)) %>%
  
  # functional beta.pair requires 3 or more species
  filter(rowSums(.[-1]) >= 3)  %>%
  
  column_to_rownames(var = "site")

# calculate functional beta-diversity ------------------------------------------

# take the first two PCoa axis
pcoa_axis2 <- trait_pcoa$vectors.cor[, 1:2]

# partition beta diversity into a list of 
# 1) total beta-diversity dissimilarity matrix
# 2) nestedness dissimilarity matrix
# 3) turnover dissimilarity matrix
func_beta <- functional.beta.pair(x = comm_rel, 
                                  traits = pcoa_axis2, 
                                  index.family = "sorensen")

# Save to disk -----------------------------------------------------------------

saveRDS(func_beta, here("data/working", "func_beta_matrices.rds"))