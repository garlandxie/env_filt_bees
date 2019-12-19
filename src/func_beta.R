#######################################################################################################
# calculate functional beta diversity (Dpw)
#######################################################################################################

# libraries --------------------------------------------------------------------
library(here)
library(vegan)
library(betapart)
library(dplyr)
library(tibble)
library(tidylog)
library(cluster)
library(ape)

# import -----------------------------------------------------------------------
comm <- read.csv(here("data/original", "community_data_matrix.csv"),
                      row.names = 1)

trait <- readRDS(here("data/final", "trait_matrix.rds"))

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

func_beta <- functional.beta.pair(x = comm_rel, 
                                  traits = pcoa_axis2, 
                                  index.family = "sorensen")

# Save to disk -----------------------------------------------------------------

saveRDS(func_beta, here("data/working", "func_beta_matrices.rds"))




