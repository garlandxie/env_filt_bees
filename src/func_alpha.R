#######################################################################################################
# calculate functional alpha diversity (ses.mfd)
#######################################################################################################

# libraries --------------------------------------------------------------------
library(StatMatch)
library(picante)
library(here)
library(tidyverse)
library(cowplot)
library(factoextra)

# import -----------------------------------------------------------------------
comm <- read.csv(here("data/original", "community_data_matrix.csv"),
                 row.names = 1)

trait <- read.csv(here("data/final", "trait_matrix.csv"))

# check packaging --------------------------------------------------------------

# community data matrix
str(comm)
head(comm, n = 5)
tail(comm, n = 5)

# trait matrix
str(trait)
head(trait, n = 5)

# data cleaning ----------------------------------------------------------------

# to relative abundance
comm_rel <- decostand(comm, method = "total")

# assign appropriate data types to each trait
trait_tidy <- trait %>%
  mutate(emer_time = factor(emer_time, ordered = TRUE),
         num_nest_mat = factor(num_nest_mat, ordered = TRUE),
         volt = factor(volt)) %>%
  column_to_rownames(var = "spp") %>%
  select(-X)

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

# some prep
ses_MFD <- partial(ses.mpd, 
                   samp = comm_tidy2, 
                   dis  = trait_dist, 
                   runs = 999, 
                   abundance.weighted = TRUE)

# lots of null models!
ses_MFD_tx <- ses_MFD(null.model = "taxa.labels")
ses_MFD_R  <- ses_MFD(null.model = "richness")
ses_MFD_freq <- ses_MFD(null.model = "frequency")
ses_MFD_sp   <- ses_MFD(null.model = "sample.pool")
ses_MFD_pp   <- ses_MFD(null.model = "phylogeny.pool")
ses_MFD_is   <- ses_MFD(null.model = "independentswap")
ses_MFD_ts   <- ses_MFD(null.model = "trialswap")

# histograms -------------------------------------------------------------------

hist(ses_MFD_freq$mpd.obs.p)
hist(ses_MFD_is$mpd.obs.p)
hist(ses_MFD_ts$mpd.obs.p)
hist(ses_MFD_pp$mpd.obs.p)
hist(ses_MFD_tx$mpd.obs.p)
hist(ses_MFD_R$mpd.obs.p)

# linear models ----------------------------------------------------------------

# linear models
lm_mfd <- partial(lm, formula = mpd.obs.z ~ ntaxa)

summary(lm_mfd(data = ses_MFD_tx))
summary(lm_mfd(data = ses_MFD_R))
summary(lm_mfd(data = ses_MFD_freq))
summary(lm_mfd(data = ses_MFD_sp))
summary(lm_mfd(data = ses_MFD_pp))
summary(lm_mfd(data = ses_MFD_is))
summary(lm_mfd(data = ses_MFD_ts))

# plots: ses.MFD versus richness -----------------------------------------------

MFD_sr <- function(df, mod_title = NULL) {
  
  if(is.data.frame(df)) {
    df %>%
      filter(ntaxa != 0 & ntaxa != 1) %>%
      ggplot(aes(x = ntaxa, y = mpd.obs.z)) + 
      geom_point() + 
      geom_smooth(method = "lm") + 
      labs(x = "Species Richness",
           y = "SES.MFD",
           title = as.character(mod_title)) + 
      theme_minimal()
  } else {
    print("df input is not a data-frame")
  }
}

# individual plots
plot_tx <- MFD_sr(ses_MFD_tx, mod_title = "taxa labels")
plot_R  <- MFD_sr(ses_MFD_R, mod_title = "richness")
plot_freq <- MFD_sr(ses_MFD_R, mod_title = "frequency")
plot_sp <- MFD_sr(ses_MFD_sp, mod_title = "speciess pool")
plot_pp <- MFD_sr(ses_MFD_pp, mod_title = "phylogeny pool")
plot_is <- MFD_sr(ses_MFD_is, mod_title = "ind swap")
plot_ts <- MFD_sr(ses_MFD_ts, mod_title = "trial swap")

# multi-panel figure
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
  tx = ses_MFD_tx$mpd.obs.z,
  f  = ses_MFD_freq$mpd.obs.z, 
  is = ses_MFD_is$mpd.obs.z, 
  ts = ses_MFD_ts$mpd.obs.z, 
  pp = ses_MFD_pp$mpd.obs.z, 
  sp = ses_MFD_sp$mpd.obs.z, 
  r  = ses_MFD_R$mpd.obs.z) %>%
  filter(complete.cases(.))

# looks like I should keep the simplest null model then
# aka taxa labels
pca <- prcomp(pca_df, center = TRUE, scale = TRUE)
fviz_pca_var(pca, repel = TRUE)

# save to disk -----------------------------------------------------------------

# rename columns to reflect mfd not mpd 
colnames(ses_MFD_tx) <- str_replace(colnames(ses_MFD_tx), "mpd", "mfd")

# write to disk
write.csv(ses_MFD_tx, here("data/working", "d99_ses_mfd.csv"))

