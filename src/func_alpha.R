# calculate functional alpha diversity (ses.mfd) -------------------------------
# author(s): Nicholas Sookhan and Garland Xie
# institutional affiliation: University of Toronto 

# libraries --------------------------------------------------------------------
library(StatMatch) # for calculating gower's distance
library(picante)   # for analyzing community matrices
library(here)      # for creating relative file-paths
library(tidyverse) 

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
tail(trait,  n = 5)

# data cleaning ----------------------------------------------------------------

# to relative abundance
comm_rel <- decostand(comm, method = "total")

# assign species as row names
trait_tidy <- trait %>%
  column_to_rownames(var = "spp") 

# calculate gower distance matrix
trait_dist <- gower.dist(trait_tidy)
colnames(trait_dist) <- rownames(trait_tidy)
rownames(trait_dist) <- rownames(trait_tidy)

# drop species from community data matrix that we do not have trait data for 
comm_tidy2 <- comm %>%
  select(rownames(trait_tidy))

# double checks
all(trait_dist >= 0) 
isSymmetric.matrix(trait_dist)

# calculate ses.mfd ------------------------------------------------------------

# use a taxa swap model 
ses_MFD_tx <- ses.mpd(samp = comm_tidy2, 
                      dis  = trait_dist, 
                      runs = 999, 
                      abundance.weighted = TRUE,
                      null.model = "taxa.labels")

# linear models ----------------------------------------------------------------

summary(lm(data = ses_MFD_tx, formula = mpd.obs.z ~ ntaxa))

# plots: ses.MFD versus richness -----------------------------------------------

(plot_tx <- ses_MFD_tx %>%
  filter(ntaxa != 0 & ntaxa != 1) %>%
  ggplot(aes(x = ntaxa, y = mpd.obs.z)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Species Richness",
       y = "SES.MFD") + 
  theme_minimal())

# Visualize histograms of p-values ---------------------------------------------

# rename columns to reflect mfd not mpd 
colnames(ses_MFD_tx) <- str_replace(colnames(ses_MFD_tx), "mpd", "mfd")

(p_hist <- ses_MFD_tx %>%
   rename(ses_mfd = mfd.obs.z, 
          p_value = mfd.obs.p) %>%
   filter(ses_mfd < 0) %>%
   ggplot(aes(x = p_value)) +
   geom_histogram(bins = 30) + 
   geom_vline(xintercept = 0.05, linetype = "dashed") + 
   labs(y = "Frequency", 
        x = "Randomized p-values of ses.MFD") + 
   theme_minimal())

# save to disk -----------------------------------------------------------------

# write to disk
write.csv(ses_MFD_tx, here("data/working", "ses_mfd.csv"))

# Save histogram of p-values from ses.MFD
ggsave(filename = here("output/figures", "fig1_crit1_ses_mfd.png"),
       plot = p_hist,
       units = "in",
       height = 4,
       width = 4)

