#######################################################################################################
# calculate phylogenetic alpha diversity (ses.mpd)
#######################################################################################################

# libraries --------------------------------------------------------------------
library(ape)
library(picante)
library(vegan)
library(here)
library(dplyr)
library(ggplot2)
library(purrr)
library(cowplot)
library(factoextra)

# import -----------------------------------------------------------------------

# convert column ID into row-names 
# required for community ecology analyses
comm <- read.csv(here("data/original", "community_data_matrix.csv"),
                 row.names = 1)

# phylogenetic tree
tree <- read.tree(here("data/final", "tree_rescale_ulm.new"))

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

# so many null models
ses_mpd <- partial(ses.mpd, 
                   samp = comm_rel, 
                   dis  = pw_phylo, 
                   runs = 999, 
                   abundance.weighted = TRUE)

# calculate ses.mpd according to a specific null model
ses_mpd_tx <- ses_mpd(null.model = "taxa.labels")
ses_mpd_R  <- ses_mpd(null.model = "richness")
ses_mpd_freq <- ses_mpd(null.model = "frequency")
ses_mpd_sp   <- ses_mpd(null.model = "sample.pool")
ses_mpd_pp   <- ses_mpd(null.model = "phylogeny.pool")
ses_mpd_is   <- ses_mpd(null.model = "independentswap")
ses_mpd_ts   <- ses_mpd(null.model = "trialswap")

# histogram of p-values from each null model -----------------------------------

hist(ses_mpd_freq$mpd.obs.p)
hist(ses_mpd_is$mpd.obs.p)
hist(ses_mpd_ts$mpd.obs.p)
hist(ses_mpd_pp$mpd.obs.p)
hist(ses_mpd_tx$mpd.obs.p)
hist(ses_mpd_R$mpd.obs.p)

# linear models: ses.mpd versus species richness -------------------------------

# linear models
lm_mpd <- partial(lm, formula = mpd.obs.z ~ ntaxa)

summary(lm_mpd(data = ses_mpd_tx))
summary(lm_mpd(data = ses_mpd_R))
summary(lm_mpd(data = ses_mpd_freq))
summary(lm_mpd(data = ses_mpd_sp))
summary(lm_mpd(data = ses_mpd_pp))
summary(lm_mpd(data = ses_mpd_is))
summary(lm_mpd(data = ses_mpd_ts))

# plots: ses.mpd vs species richness -------------------------------------------

mpd_sr <- function(df, mod_title = NULL) {
  
  if(is.data.frame(df)) {
    df %>%
      filter(ntaxa != 0 & ntaxa != 1) %>%
      ggplot(aes(x = ntaxa, y = mpd.obs.z)) + 
      geom_point() + 
      geom_smooth(method = "lm") + 
      labs(x = "Species Richness",
           y = "SES.MPD",
           title = as.character(mod_title)) + 
      theme_minimal()
  } else {
    print("df input is not a data-frame")
  }
}
  
plot_tx <- mpd_sr(ses_mpd_tx, mod_title = "taxa labels")
plot_R  <- mpd_sr(ses_mpd_R, mod_title = "richness")
plot_freq <- mpd_sr(ses_mpd_R, mod_title = "frequency")
plot_sp <- mpd_sr(ses_mpd_sp, mod_title = "speciess pool")
plot_pp <- mpd_sr(ses_mpd_pp, mod_title = "phylogeny pool")
plot_is <- mpd_sr(ses_mpd_is, mod_title = "ind swap")
plot_ts <- mpd_sr(ses_mpd_ts, mod_title = "trial swap")

(multi_pan <- plot_grid(plot_freq,
                       plot_is,
                       plot_ts,
                       plot_R,
                       plot_sp, 
                       plot_pp,
                       plot_tx,
                       ncol = 4))

# PCA --------------------------------------------------------------------------

# assuming all vectors have the same order of sites
pca_df <- data.frame(
  tx = ses_mpd_tx$mpd.obs.z,
  f  = ses_mpd_freq$mpd.obs.z, 
  is = ses_mpd_is$mpd.obs.z, 
  ts = ses_mpd_ts$mpd.obs.z, 
  pp = ses_mpd_pp$mpd.obs.z, 
  sp = ses_mpd_sp$mpd.obs.z, 
  r  = ses_mpd_R$mpd.obs.z) %>%
  filter(complete.cases(.))

# looks like I should keep the simplest null model then
# aka taxa labels
pca <- prcomp(pca_df, center = TRUE, scale = TRUE)
fviz_pca_var(pca, repel = TRUE)

# Save the file ----------------------------------------------------------------

write.csv(ses_mpd_tx,
          here("data/working", "d99_ses_mpd.csv"))

