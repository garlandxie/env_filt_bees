# calculate phylogenetic alpha diversity (ses.mpd) -----------------------------
# author(s): Garland Xie, Nicholas Sookhan

# libraries --------------------------------------------------------------------
library(ape)     # for reading in phylogenetic trees
library(picante) # for calculating ses.mpd values
library(vegan)   # for analysing community matrices
library(here)    # for creating relative file-paths
library(dplyr)   # for manipulating data 
library(ggplot2) # for visualising data 

# import -----------------------------------------------------------------------

# relative file-paths
comm_path <- here("data/original", "community_data_matrix.csv")
tree_path <- here("data/final", "tree_rescale_ulm.new")

# convert column ID into row-names 
# required for community ecology analyses
comm <- read.csv(comm_path, row.names = 1)

# phylogenetic tree
tree <- read.tree(tree_path)

# check packaging --------------------------------------------------------------
str(comm)
head(comm, n = 5)
tail(comm, n = 5)

# data cleaning ----------------------------------------------------------------

# convert raw to relative abundance for each species
comm_rel <- decostand(comm, method = "total")
  
# get pairwise phylogenetic distances
pw_phylo <- cophenetic(tree)

# double-checks: 
# (1) all pairwise distances are equal or above zero
# (2) pairwise distance matrix should be symmetric
all(pw_phylo >= 0) 
isSymmetric.matrix(pw_phylo)

# calculate ses.mpd ------------------------------------------------------------

# calculate ses.mpd according to a specific null model
ses_mpd_tx  <- ses.mpd(samp = comm_rel, 
                       dis = pw_phylo, 
                       runs = 999, 
                       abundance.weighted = TRUE)

# linear models: ses.mpd versus species richness -------------------------------

summary(lm(formula = mpd.obs.z ~ ntaxa, data = ses_mpd_tx))

# plots: ses.mpd vs species richness -------------------------------------------
  
(plot_tx <- ses_mpd_tx %>%
  filter(ntaxa != 0 & ntaxa != 1) %>%
  ggplot(aes(x = ntaxa, y = mpd.obs.z)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Species Richness",
       y = "SES.MPD") +
  theme_minimal())

# Visualize histogram of p-values ----------------------------------------------

(p_hist <- ses_mpd_tx %>%
  rename(ses_mpd = mpd.obs.z, 
         p_value = mpd.obs.p) %>%
  filter(ses_mpd < 0) %>%
  ggplot(aes(x = p_value)) +
  geom_histogram(bins = 30) + 
  geom_vline(xintercept = 0.05, linetype = "dashed") + 
  labs(y = "Frequency", 
       x = "p-values of under-dispersed ses.MPD") + 
  theme_minimal())

# Save the data! ---------------------------------------------------------------

# ses mpd csv file
write.csv(ses_mpd_tx, here("data/working", "ses_mpd.csv"))

# figure: histogram of ses.mpd values
ggsave(filename = here("output/figures", "fig1_crit1_ses_mpd.png"),
       plot = p_hist,
       units = "in",
       height = 4,
       width = 4)

