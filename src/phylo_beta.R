#######################################################################################################
# calculate phylogenetic beta diversity (Dpw)
#######################################################################################################

# libraries --------------------------------------------------------------------
library(ape)
library(picante)
library(vegan)
library(here)
library(dplyr)
library(betapart)

# import -----------------------------------------------------------------------
comm <- read.csv(here("data/original", "community_data_matrix.csv"), 
                 row.names = 1)

tree <- read.tree(here("data/original", "phylo_tree_ulm.new"))

# data cleaning - community data matrix ----------------------------------------

# vector of kleptoparasite names
kleptos <- c("Coelioxys_alternata", 
             "Coelioxys_moesta",       
             "Coelioxys_sayi",
             "Stelis_verticalis")
             
# presence or absence data 
comm_rel <- decostand(comm, method = "pa") %>%
  
  # change to column b/c of tidyverse stuff
  rownames_to_column(var = "site") %>%
  
  # remove kleptoparasites
  select(-kleptos) %>%
  
  # functional beta.pair requires 3 or more species
  # so make sure both matrices end up with the same sites
  filter(rowSums(.[-1]) >= 3)  %>%
  
  # change back to rownames
  column_to_rownames(var = "site")

# data cleaning: phylo tree ----------------------------------------------------

phylo_pruned <- drop.tip(tree, 
                         tip = kleptos)

# phylo beta diversity ---------------------------------------------------------

phylo_beta <- phylo.beta.pair(x = comm_rel, 
                              tree = phylo_pruned, 
                              index.family = "sorensen")

# save to disk -----------------------------------------------------------------

saveRDS(phylo_beta, file = here("data/working", "phylo_beta_matrices.rds"))
