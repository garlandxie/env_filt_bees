# RLQ Analysis - Xie et al. 2018 
# Script is used to clean data for multiple databases before conducting an RLQ analysis
# Code developed by Garland Xie

# Load libraries ---------------------------------------------------------------
library(dplyr)
library(ade4)
library(readxl)
library(here)
library(janitor)
library(stringr)
library(readr)
library(ggplot2)
library(ggrepel)

# Import files -----------------------------------------------------------------

# community data matrix 
comm <- read_csv(here("data/original", "community_data_matrix.csv"))

# land-use metrics for each site (250m spatial scale)
met_250 <- read_csv(here("data/original", "land_use_metrics_250.csv"))

# land-use metrics for each site (500m spatial scale)
met_500 <- read_csv(here("data/original", "land_use_metrics_500.csv"))

# trait matrix of 31 species (excluded kleptoparasites)
trait <- read_csv(here("data/final", "trait_matrix.csv"))

# Clean: comm matrix -------------------------------------------------

# Change taxon names to abbreviated versions 
# plotting purposes
taxon_change <- function (spList) {
  
  f <- spList %>%
    # separate genus name and species name onto a list
    strsplit(split = "_") %>%
    # create abbreviations of each species 
    # take first four letter of genus (e.g. Anthidium -> Anth)
    # take first three letters of species (e.g. manicatum -> man)
    lapply(function(x) paste(substr(x[1], 1, 3), 
                             substr(x[2], 1, 3), 
                             sep = "_"))
  # return vector of abbreviated taxon names
  return(f)
}  

# clean
comm_tidy <- comm %>%
  
  # alphabetical order
  select(order(colnames(.))) %>%
  
  # concise species code
  select_all(taxon_change) %>%
  rename(ID = X1_NA,
         Anth_man = Ant_man) %>%
  
  # remove kleptoparasites
  select(-c("Ste_ver", "Coe_say", "Coe_alt", "Coe_moe"))

# Clean: trait matrix ------------------------------------------------

trait_clean <- trait %>%
  select(-"X1") %>%
  mutate(ID = taxon_change(spp) %>% 
           str_replace("Ant_man", "Anth_man")) %>%
  
  # change data types
  mutate(native_y_n     = factor(native_y_n),
         emer_time      = factor(emer_time),
         leaf_hair      = factor(leaf_hair),
         leaf_cut       = factor(leaf_cut),
         leaf_pulp      = factor(leaf_pulp),
         resin          = factor(resin),
         stone          = factor(stone),
         none           = factor(none),
         mud            = factor(mud),
         num_nest_mat   = factor(num_nest_mat),
         diet           = factor(diet),
         volt           = factor(volt)
         )

# double-check
str(trait_clean)

# Clean: environmental vars ----------------------------------------------------

# Re-order both environmental data matrices based on rows 
# Makes it easy to do value-matching 
met_250_clean <- met_250 %>%
  rename(ID = X1) %>%
  arrange(ID)

met_500_clean <- met_500 %>%
  rename(ID = X1) %>%
  arrange(ID)

# RLQ: 250 ---------------------------------------------------------------------

# Procedure of RLQ follows closely to:
# Dray et al. 2014. Ecology. 
# https://doi.org/10.1890/13-0196.1
# Supplementary Info 1

# correspondance analysis 
# community data matrix
dudiL_250 <- comm_tidy %>%
  inner_join(met_250_clean, by = "ID") %>%
  select(colnames(comm_tidy)) %>%
  tibble::column_to_rownames(var = "ID") %>%
  dudi.coa(scannf = F)

# principal component analysis
# environmental variables
# weighted by site weights from previous correspondance analysis

dudiR_250 <- met_250_clean %>%
  tibble::column_to_rownames(var = "ID") %>%
  select("250_urban" = "total.area_250_urban",
         "250_grass" = "total.area_250_grass",
         "250_tree"  = "total.area_250_tree_canopy") %>%
  dudi.pca(row.w = dudiL_250$lw, scannf = F)

# hill-smith method: PCA of discrete/continuous data 
# reference: Hill and Smith. 1976. Taxon
# weighted by sites weight from previous correspondance analysis

dudiQ_250 <- trait_clean %>%
  select(-spp) %>%
  tibble::column_to_rownames(var = "ID") %>%
  dudi.hillsmith(row.w = dudiL_250$cw, 
                 scannf = F)

# RLQ: 250m spatial scale 
# obtain first two RLQ axes in a non-iteractive manner 
RLQ_250 <- rlq(dudiR = dudiR_250, 
               dudiQ = dudiQ_250, 
               dudiL = dudiL_250,
               scannf = FALSE, nf = 2)

# cumulative projected inertia (%)
summary(RLQ_250)

# RLQ: 500 ---------------------------------------------------------------------

# L table - species abundance
# correspondance analysis 
# community data matrix
dudiL_500 <- comm_tidy %>%
  inner_join(met_500_clean, by = "ID") %>%
  select(colnames(comm_tidy)) %>%
  tibble::column_to_rownames(var = "ID") %>%
  dudi.coa(scannf = F)

# Q table - traits
# hill-smith method: PCA of discrete/continuous data 
# weighted by sites weight from previous correspondance analysis
dudiQ_500 <- trait_clean %>%
  select(-spp) %>%
  tibble::column_to_rownames(var = "ID") %>%
  dudi.hillsmith(row.w = dudiL_500$cw, 
                 scannf = F)

# R table - env vars
# principal component analysis
# weighted by site weights from previous correspondance analysis
dudiR_500 <- met_500_clean %>%
  tibble::column_to_rownames(var = "ID") %>%
  select("500_urban" = "total.area_500_urban",
         "500_grass" = "total.area_500_grass",
         "500_tree"  = "total.area_500_tree_canopy") %>%
  dudi.pca(row.w = dudiL_500$lw, scannf = F)

# RLQ - 500
# obtain first two RLQ axes in an non-iteractive manner 
RLQ_500 <- rlq(dudiR = dudiR_500, 
               dudiQ = dudiQ_500, 
               dudiL = dudiL_500,
               scannf = FALSE, nf = 2)

# cumulative projected inertia (%)
summary(RLQ_500)

# RLQ plots: species -----------------------------------------------------------

lQ_rlq_250 <- 
 
  ggplot() +
  geom_vline(xintercept = 0, 
             lwd = 0.5, 
             col = "grey") + 
  geom_hline(yintercept = 0, 
             lwd = 0.5, 
             col = "grey") +
  geom_point(aes(x = AxcQ1, y = AxcQ2), 
             size = 2, 
             color = "grey", 
             alpha = 1, 
             data = RLQ_250$lQ) +
  geom_label_repel(aes(x = AxcQ1, 
                       y = AxcQ2, 
                       fill = (str_extract(rownames(RLQ_250$lQ), "^[^_]+")),  
                       label = rownames(RLQ_250$lQ)),
                   segment.alpha = 0.5,
                   point.padding = 0.5,
                   size = 3,
                   data = RLQ_250$lQ
                   ) + 
  scale_fill_discrete(name = "Genus") + 
  xlab("Axis 1") + 
  ylab("Axis 2") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black")
  ) 

lQ_rlq_500 <- 
  
  ggplot() +
  geom_vline(xintercept = 0, 
             lwd = 0.5, 
             col = "grey") + 
  geom_hline(yintercept = 0, 
             lwd = 0.5, 
             col = "grey") +
  geom_point(aes(x = AxcQ1, y = AxcQ2), 
             size = 2, 
             color = "grey", 
             alpha = 1, 
             data = RLQ_500$lQ) +
  geom_label_repel(aes(x = AxcQ1, 
                       y = AxcQ2, 
                       fill = (str_extract(rownames(RLQ_500$lQ), "^[^_]+")),  
                       label = rownames(RLQ_500$lQ)),
                   segment.alpha = 0.5,
                   point.padding = 0.5,
                   size = 3,
                   data = RLQ_500$lQ
  ) + 
  scale_fill_discrete(name = "Genus") + 
  xlab("Axis 1") + 
  ylab("Axis 2") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black")
  ) 
  
# RLQ plots: traits ------------------------------------------------------------

c1_rlq_500 <-
  
ggplot() + 
  geom_vline(xintercept = 0, 
             lwd = 0.5, 
             col = "grey") +  
  geom_hline(yintercept = 0, 
             lwd = 0.5,
             col = "grey") +
  geom_label_repel(aes(x = CS1, y = CS2, 
                       label = rownames(RLQ_500$c1)), 
                   point.padding = 0.5,
                   data = RLQ_500$c1) + 
  geom_segment(aes(x = 0, y = 0, xend = CS1, yend = CS2),
               arrow = arrow(length = unit(0.01, "npc")),
               colour = "blue",
               alpha = 0.25, 
               data = RLQ_500$c1) +  
  xlab("Axis 1") +
  ylab("Axis 2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
  ) 


c1_rlq_500 <- 
  
  ggplot() + 
  geom_vline(xintercept = 0, 
             lwd = 0.5, 
             col = "grey") +  
  geom_hline(yintercept = 0, 
             lwd = 0.5,
             col = "grey") +
  geom_label_repel(aes(x = CS1, y = CS2, 
                       label = rownames(RLQ_250$c1)), 
                   point.padding = 0.5,
                   data = RLQ_250$c1) + 
  geom_segment(aes(x = 0, y = 0, xend = CS1, yend = CS2),
               arrow = arrow(length = unit(0.01, "npc")),
               colour = "blue",
               alpha = 0.25, 
               data = RLQ_250$c1) +  
  xlab("Axis 1") +
  ylab("Axis 2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
  ) 
  