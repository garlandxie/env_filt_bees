# RLQ analysis -----------------------------------------------------------------
# Script is used to conduct an RLQ analysis ------------------------------------
# Code developed by Garland Xie
# Institutional affiliation: University of Toronto
# Contact info: garland.xie@mail.utoronto.ca

# Load libraries ---------------------------------------------------------------
library(ade4)
library(here)
library(janitor)
library(ggrepel)
library(tidyverse)

# Import files -----------------------------------------------------------------

# relative file-paths 
comm_path <- here("data/original", "community_data_matrix.csv")
met_250_path <- here("data/original", "land_use_metrics_250.csv")
met_500_path <- here("data/original", "land_use_metrics_500.csv")
trait_path <- here("data/final", "trait_matrix.rds")

# load data
comm <- read_csv(comm_path)
met_250 <- read_csv(met_250_path)
met_500 <- read_csv(met_500_path)
trait <- readRDS(trait_path)

# Clean: community matrix ------------------------------------------------------

# Change taxon names to abbreviated versions 
# for plotting purposes
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
  mutate(ID = spp %>% 
           as.character() %>%
           taxon_change() %>%
           str_replace("Ant_man", "Anth_man")
         ) %>%
  
  rename(LH = leaf_hair, 
         LC = leaf_cut, 
         LP = leaf_pulp)

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

# double check
str(met_250_clean)
str(met_500_clean)

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
  column_to_rownames(var = "ID") %>%
  dudi.coa(scannf = F)

# principal component analysis
# environmental variables
# weighted by site weights from previous correspondance analysis

dudiR_250 <- met_250_clean %>%
  column_to_rownames(var = "ID") %>%
  select("250_urban" = "total.area_250_urban",
         "250_grass" = "total.area_250_grass",
         "250_tree"  = "total.area_250_tree_canopy") %>%
  dudi.pca(row.w = dudiL_250$lw, scannf = F)

# hill-smith method: PCA of discrete/continuous data 
# reference: Hill and Smith. 1976. Taxon
# weighted by sites weight from previous correspondance analysis

dudiQ_250 <- trait_clean %>%
  select(-spp) %>%
  column_to_rownames(var = "ID") %>%
  dudi.hillsmith(row.w = dudiL_250$cw, scannf = F)

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


c1_rlq_250 <- 
  
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

# RLQ plots: env vars ----------------------------------------------------------

env_vars <- c("% Grass", "% Tree Canopy", "% Urban Area") 
rownames(RLQ_250$l1) <- env_vars
rownames(RLQ_500$l1) <- env_vars

l1_rlq_250 <- 
  ggplot() +
    geom_vline(xintercept = 0, lwd = 0.5, col = "grey") +  
    geom_hline(yintercept = 0, lwd = 0.5,col = "grey") +
    geom_label_repel(aes(x = RS1, 
                         y = RS2, 
                         label = rownames(RLQ_250$l1)), 
                     point.padding = 0.5,
                     data = RLQ_250$l1) +
    geom_segment(aes(x = 0, 
                     y = 0, 
                     xend = RS1, yend = RS2),
                 arrow = arrow(length = unit(0.01, "npc")),
                 colour = "red",
                 size = 1,
                 alpha = 0.25, 
                 data = RLQ_250$l1) +  
  xlab("Axis 1") +
  ylab("Axis 2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

l1_rlq_500 <- 
  ggplot() +
  geom_vline(xintercept = 0, lwd = 0.5, col = "grey") +  
  geom_hline(yintercept = 0, lwd = 0.5,col = "grey") +
  geom_label_repel(aes(x = RS1, 
                       y = RS2, 
                       label = rownames(RLQ_500$l1)), 
                   point.padding = 0.5,
                   data = RLQ_500$l1) +
  geom_segment(aes(x = 0, 
                   y = 0, 
                   xend = RS1, yend = RS2),
               arrow = arrow(length = unit(0.01, "npc")),
               colour = "red",
               size = 1,
               alpha = 0.25, 
               data = RLQ_500$l1) +  
  xlab("Axis 1") +
  ylab("Axis 2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) 

# RLQ plot: eigenvalues --------------------------------------------------------

eigs_250 <- data.frame(
  Axes = paste("Axis", 1:3, sep = " "),
  Eigenvalues =  RLQ_250$eig,
  Projected.Inertia = round(RLQ_250$eig/sum(RLQ_250$eig)*100, digits = 2)
  )

eigs_500 <- data.frame(
  Axes = paste("Axis", 1:3, sep = " "),
  Eigenvalues =  RLQ_500$eig,
  Projected.Inertia = round(RLQ_500$eig/sum(RLQ_500$eig)*100, digits = 2)
)

plot_eigs_250 <- 
  ggplot(data = eigs_250, aes(x = Axes, y = Projected.Inertia)) +
  geom_bar(colour = "black", stat = "identity") + 
  labs(x = "Axes", y = "Projected Inertia (%)") + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
  ) 

plot_eigs_500 <- 
  ggplot(data = eigs_500, aes(x = Axes, y = Projected.Inertia)) +
  geom_bar(colour = "black", stat = "identity") + 
  labs(x = "Axes", y = "Projected Inertia (%)") + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
  ) 

# misc plots -------------------------------------------------------------------

R_250_load <- RLQ_250$l1 %>%
  rownames_to_column(var = "class") %>%
  select(class, RS1) %>%
  mutate(class = factor(class) %>% fct_reorder(RS1)) %>%
  ggplot(aes(x = class, y = RS1, fill = class)) + 
  geom_col(show.legend = FALSE, alpha = 0.8) +
  labs(x = NULL,
       y = "Relative importance in principle component 1") + 
  coord_flip() 
  
R_500_load <- RLQ_500$l1 %>%
  rownames_to_column(var = "class") %>%
  select(class, RS1) %>%
  mutate(class = factor(class) %>% fct_reorder(RS1)) %>%
  ggplot(aes(x = class, y = RS1, fill = class)) + 
  geom_col(show.legend = FALSE, alpha = 0.8) +
  labs(x = NULL,
       y = "Relative importance in principle component (weighted by species abundance") + 
  coord_flip() 

traits_N <- c("stone.0","none.0","resin.0", "mud.0", 
              "LH.0", "LP.0", "LC.0")
 
RLQ_250_load <- RLQ_250$c1 %>%
  rownames_to_column(var = "traits")%>%
  select(traits, CS1) %>%
  filter(!traits %in% traits_N) %>%
  mutate(traits = as.character(traits), 
         traits = case_when(
           traits == "emer_.Summer" ~ "Emergence Time (Summer)", 
           traits == "emer_.Spring"~ "Emergence Time (Spring)", 
           traits == "diet.Poly" ~ "Diet Breadth (Polylectic)",
           traits == "diet.Oligo" ~ "Diet Breadth (Oligolectic)",
           traits == "volt.1" ~ "Univoltinism",
           traits == "volt.2" ~ "Multivolinism",
           traits == "stone.1" ~ "Nesting Material (Stone)",
           traits == "none.1" ~ "Nesting Material (None)",
           traits == "mud.1" ~ "Nesting Material (Mud)",
           traits == "LH.1" ~ "Nesting Material (Leaf Hair)",
           traits == "LC.1" ~ "Nesting Material (Leaf Cut)",
           traits == "LP.1" ~ "Nesting Material (Leaf Pulp)",
           traits == "resin.1" ~ "Nesting Material (Resin)",
           traits == "nativ.1"  ~ "Native Status",
           traits == "nativ.0" ~ "Exotic Status",
           traits == "itd" ~ "Intertegular Distance",
           TRUE ~ traits),
           traits = reorder(traits, CS1)) %>%
  ggplot(aes(x = traits, y = CS1)) + 
  geom_col(show.legend = FALSE, alpha = 0.8) +
  labs(x = NULL,
       y = "Relative importance in RLQ component 1") + 
  coord_flip() + 
  theme_minimal()

RLQ_500_load <- RLQ_500$c1 %>%
  rownames_to_column(var = "traits")%>%
  select(traits, CS1) %>%
  filter(!traits %in% traits_N) %>%
  mutate(traits = as.character(traits), 
         traits = case_when(
           traits == "emer_.Summer" ~ "Emergence Time (Summer)", 
           traits == "emer_.Spring"~ "Emergence Time (Spring)", 
           traits == "diet.Poly" ~ "Diet Breadth (Polylectic)",
           traits == "diet.Oligo" ~ "Diet Breadth (Oligolectic)",
           traits == "volt.1" ~ "Univoltinism",
           traits == "volt.2" ~ "Multivolinism",
           traits == "stone.1" ~ "Nesting Material (Stone)",
           traits == "none.1" ~ "Nesting Material (None)",
           traits == "mud.1" ~ "Nesting Material (Mud)",
           traits == "LH.1" ~ "Nesting Material (Leaf Hair)",
           traits == "LC.1" ~ "Nesting Material (Leaf Cut)",
           traits == "LP.1" ~ "Nesting Material (Leaf Pulp)",
           traits == "resin.1" ~ "Nesting Material (Resin)",
           traits == "nativ.1"  ~ "Native Status",
           traits == "nativ.0" ~ "Exotic Status",
           traits == "itd" ~ "Intertegular Distance",
           TRUE ~ traits),
         traits = reorder(traits, CS1)) %>%
  ggplot(aes(x = traits, y = CS1)) + 
  geom_col(show.legend = FALSE, alpha = 0.8) +
  labs(x = NULL,
       y = "Relative importance in RLQ component 1") + 
  coord_flip() + 
  theme_minimal()

L_250_load <- RLQ_250$mQ %>%
  rownames_to_column(var = "species") %>%
  select(species, NorS1) %>%
  mutate(species = reorder(species, NorS1),
         genus = str_split(species, pattern = "_") %>%
                sapply("[", 1)) %>%
  ggplot(aes(x = species, y = NorS1, fill = genus)) + 
  geom_col(show.legend = TRUE, alpha = 0.8) +
  labs(x = NULL,
       y = "Relative importance in CA component 1") + 
  coord_flip()

L_500_load <- RLQ_500$mQ %>%
  rownames_to_column(var = "species") %>%
  select(species, NorS1) %>%
  mutate(species = factor(species) %>% fct_reorder(NorS1),
         genus = str_split(species, pattern = "_") %>%
           sapply("[", 1)) %>%
  ggplot(aes(x = species, y = NorS1, fill = genus)) + 
  geom_col(show.legend = TRUE, alpha = 0.8) +
  labs(x = NULL,
       y = "Relative importance in CA component") + 
  coord_flip() + 
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(color = "white"))

