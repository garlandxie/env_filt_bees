# Partial Procrustes Analysis --------------------------------------------------

# Goal of script: to determine a relationship between:
# phylogenetic beta diversity (total, nest, turnover) and 
# environmental distance

# Author(s): Nicholas Sookhan, Garland Xie
# Institutional affiliation: University of Toronto

# libraries --------------------------------------------------------------------
library(vegan)
library(picante)
library(here)
library(readr)
library(dplyr)

# import -----------------------------------------------------------------------

# relative file-paths
env_250_path <- here("data/original", "land_use_metrics_250.csv")
env_500_path <- here("data/original", "land_use_metrics_500.csv")
phylo_path <- here("data/working", "phylo_beta_matrices.rds")
dist_spa_path <- here("data/original", "spatial_distance_matrix.csv")

# load data
l_met_250 <- read_csv(env_250_path)
l_met_500 <- read_csv(env_500_path)
phylo_beta <- readRDS(phylo_path)
dist_spa <- read_csv(dist_spa_path)

# check packaging --------------------------------------------------------------

# environmental data
str(l_met_250)
head(l_met_250, n = 5)
tail(l_met_250, n = 5)

str(l_met_500)
head(l_met_500, n = 5)
tail(l_met_500, n = 5)

# spatial distance matrix
str(dist_spa)
head(dist_spa, head = 5)
tail(dist_spa, head = 5)

# phylo beta diversity matrix
str(as.matrix(phylo_beta$phylo.beta.sor)) 
head(as.matrix(phylo_beta$phylo.beta.sor))
tail(as.matrix(phylo_beta$phylo.beta.sor))

# phylo beta diversity matrix (spatial turn-over)
str(as.matrix(phylo_beta$phylo.beta.sim)) 
head(as.matrix(phylo_beta$phylo.beta.sim))
tail(as.matrix(phylo_beta$phylo.beta.sim))

# phylo beta diversity matrix (nestedness)
str(as.matrix(phylo_beta$phylo.beta.sne)) 
head(as.matrix(phylo_beta$phylo.beta.sne))
tail(as.matrix(phylo_beta$phylo.beta.sne))

# clean data: phylo beta diversity ---------------------------------------------

# spatial turnover
phylo_tu <- as.matrix(phylo_beta$phylo.beta.sim)

# nestedness
phylo_ne <- as.matrix(phylo_beta$phylo.beta.sne)

# total beta diversity 
phylo_tot <- as.matrix(phylo_beta$phylo.beta.sor)

# verify equal dimension sizes (rows and columns)
all(dim(phylo_ne), dim(phylo_tot), dim(phylo_tu))

# verify equal rownames
all(rownames(phylo_ne) == rownames(phylo_tot))
all(rownames(phylo_tu) == rownames(phylo_tot))
all(rownames(phylo_ne) == rownames(phylo_tu))

# verify symmetrical matrices
isSymmetric.matrix(phylo_ne)
isSymmetric.matrix(phylo_tu)
isSymmetric.matrix(phylo_tot)

# verify phylogenetic distances are equal or above zero
all(phylo_tu >= 0)
all(phylo_ne >= 0)
all(phylo_tot >= 0)

# environmental distance -------------------------------------------------------

# subset diversity data to include sites using either the 250m or 500m buffer
dpw_phylo_250 <- dpw_phylo[rownames(dpw_phylo) %in% rownames(l_met_250), rownames(dpw_phylo) %in% rownames(l_met_250)]
dpw_phylo_500 <- dpw_phylo[rownames(dpw_phylo) %in% rownames(l_met_500), rownames(dpw_phylo) %in% rownames(l_met_500)]

# clean data: environmental data -----------------------------------------------

# subset metric data to include sites with diversity data
l_met_250 <- 
   
   
   l_met_250[rownames(l_met_250) %in% rownames(dpw_phylo_250),]
l_met_500 <- l_met_500[rownames(l_met_500) %in% rownames(dpw_phylo_500),]

# subset metric data to include columns used in analysis
l_met_250 <- l_met_250[,c("prop.landscape_250_grass", "prop.landscape_250_tree_canopy", "prop.landscape_250_urban", "edge.density_250_total")]
l_met_500 <- l_met_500[,c("prop.landscape_500_grass", "prop.landscape_500_tree_canopy", "prop.landscape_500_urban", "edge.density_500_total")]

# clean data: spatial distance  ------------------------------------------------

# subset spatial distance matrix to include sites used in analysis
dist_spa_250 <- dist_spa[rownames(dist_spa)%in%rownames(l_met_250),colnames(dist_spa)%in%rownames(l_met_250)]
dist_spa_500 <- dist_spa[rownames(dist_spa)%in%rownames(l_met_500),colnames(dist_spa)%in%rownames(l_met_500)]

#######################################################################################################
### prepare spatial, environmental and diversity data for partial procrustes ###

# PCA of spatial distance matrix and take scores of first 3 PCs for the 250m and the 500m buffer
pc_env_250 <- scores( prcomp(dist(l_met_250),scale=T),display = "sites", choices=1:3)
pc_env_500 <- scores( prcomp(dist(l_met_500),scale=T),display = "sites", choices=1:3)

# PCA of spatial distance matrix and take scores of first 3 PCs for the 250m and the 500m buffer
pc_spa_250 <- scores( prcomp(dist_spa_250,scale=T), choices=1:3)
pc_spa_500 <- scores( prcomp(dist_spa_500,scale=T), choices=1:3)

# PCA of phylogenetic beta diversity and take scores of first 3 PCs for the 250m and the 500m buffer
pc_dpw_250 <- scores( prcomp(dpw_phylo_250,scale=T), choices=1:3)
pc_dpw_500 <- scores( prcomp(dpw_phylo_500,scale=T), choices=1:3)

#######################################################################################################
### partial procrustes analysis ###
# control for space (250m buffer)
resid_spa_dpw_250 <- resid(lm(pc_dpw_250~pc_spa_250))
resid_spa_env_250 <- resid(lm(pc_env_250~pc_spa_250))
# control for space (500m buffer)
resid_spa_dpw_500 <- resid(lm(pc_dpw_500~pc_spa_500))
resid_spa_env_500 <- resid(lm(pc_env_500~pc_spa_500))

# run analysis for both the 250m and the 500m buffer
parpro_250 <- protest(resid_spa_env_250,resid_spa_dpw_250)
parpro_500 <- protest(resid_spa_env_500,resid_spa_dpw_500)



#######################################################################################################
### summarize
parpro_250
parpro_500

#######################################################################################################
### plot
setwd(paste0(wd,"/d99/results/procrustes"))

pdf('d99_procrustes_phylo.pdf',width=8,height=6)
par(mfrow=c(1,2),mar=c(4,4,1,1))
# 250
plot(as.dist(dpw_phylo_250)~dist(l_met_250),
     col="grey",pch=16,cex=0.7,
     xlab="Environmental Distance (250m Buffer)", ylab=expression("Phylogenetic Beta Diversity"*" ("*italic(D)*scriptstyle(pw)*")"), cex.lab=0.8 )
mtext("a", side=3, line=-2, adj=0.05, cex=1.5)
r2 <- paste0(signif(parpro_250$t0, 2),",")
p <- paste0(signif(parpro_250$signif, 2))
mtext(bquote(R^2 == .(r2) ~ p ~ '>' ~ 0.05),
      cex=0.8, side = 3, line = -2, adj = 0.90)
# 500
plot(as.dist(dpw_phylo_500)~dist(l_met_500),
     col="grey",pch=16,cex=0.7,
     xlab="Environmental Distance (500m Buffer)", ylab=expression("Phylogenetic Beta Diversity"*" ("*italic(D)*scriptstyle(pw)*")"), cex.lab=0.8 )
mtext("b", side=3, line=-2, adj=0.05, cex=1.5)
r2 <- paste0(signif(parpro_500$t0, 2),",")
p <- paste0(signif(parpro_500$signif, 2))
mtext(bquote(R^2 == .(r2) ~ p ~ '>' ~ 0.05),
      cex=0.8, side = 3, line = -2, adj = 0.90)
dev.off()