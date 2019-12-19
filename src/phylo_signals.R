# libraries --------------------------------------------------------------------
library(here) 
library(dplyr)
library(janitor)
library(stringr)
library(ape)
library(geiger)
library(caper)

# import -----------------------------------------------------------------------

# trait matrix
trait <- read.csv(here("data/final", "trait_matrix.csv"))

# phylo tree (with soft polytomies)
tree <- read.tree(here("data/original", "phylo_tree_ulm.new"))

# check packaging --------------------------------------------------------------

# raw trait matrix
glimpse(trait)
head(trait, n = 5)
tail(trait, n = 5)

# basic information on the phylogeny
plot(tree)

is.rooted(tree) 
is.ultrametric(tree) 

tree$edge
tree$Nnode
tree$node.label
tree$edge.length
cbind(tree$edge, tree$edge.length)
tree$tip.label

# tree scaling -----------------------------------------------------------------

# imported tree is not ultrametric - requires tree scaling
# estimating divergence times using penalized maximum likelihood 
# REF: Sanderson. 2002. Molecular Biology and Evolution. 19(1): 101-109

# outlined steps below are closely adapted from: 
# Cadotte and Davies. 2015. Phylogenies in Ecology 
# Chapter 2

# first, we need to obtain three set of parameters from a chronos object class
# (1) penalized log-likelihood
# (2) rate transformation for each edge in the phylogeny
# (3) lambda values 

# create a gradient of lambda values
l <- c(0:10000)

# specific an empty object to hold plogLik values
LL.out <- NULL

# run a for loop to contain logLik values
for (i in 1:length(l)) {
  LL.out[i] <- attributes(chronos(tree,lambda=l[i]))$ploglik
}

# select maximum plogLik value and corresponding lambda
# lambda value was 0
lambda.max <- l[LL.out==max(LL.out)]

# double check: plogLik vs log(lambda)
plot(log(l+1), LL.out,
     type = "l", lwd = 3, col = "gray",
     xlab = expression(paste("Log(",lambda, ")")),
     ylab = "Log likelihood") 

# rescale tree using lambda.max value
tree_rescale <- chronos(tree, lambda = lambda.max)

# more double-checks
plot(tree_rescale)
is.ultrametric(tree_rescale)

# save tree
write.tree(tree_rescale, here("data/final", "tree_rescale_ulm.new"))

# read in tree as a phylo object class
tree_rescale_ulm <- read.tree(here("data/final", "tree_rescale_ulm.new"))

# remove node labels on phylogenetic tree
tree_nodes_null <- tree_rescale_ulm 
tree_nodes_null$node.label <- NULL

# phylogenetic signals: ITD ----------------------------------------------------
  
# prep
ITD <- trait_clean$itd
names(ITD) <- trait_clean$spp
tree_di2multi <- multi2di(tree_rescale_ulm) # binary tree

# trait evolution models 
lambda.itd <- fitContinuous(tree_di2multi, ITD, model = "lambda")
star.tree.itd <- rescale(tree_di2multi, "lambda", 0)
lambda.star.itd <- fitContinuous(star.tree.itd, ITD, model = "BM")

# likelihood ratio test
# p < 0.05: infer phylogenetic structure for the given focal trait
LR.itd <-2*(lambda.itd$opt$lnL - lambda.star.itd$opt$lnL)
pchisq(LR.itd, df = 1, lower.tail = F)

# phylogenetic signals: native status ------------------------------------------

# fit model using phylogenetic D statistic
# statistical inference: 0 = BM model, 1 = no signal
# result: follows a BM model
phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = native_y_n)

# phylogenetic signal: leaf hair -----------------------------------------------

# fit model using phylogenetic D statistic
# statistical inference: 0 = BM model, 1 = no signal
# result: follows a BM model
phylo.d(data = as.data.frame(trait),
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = leaf_hair)

# phylogenetic signal: leaf cut -------------------------------------------------

# fit model using phylogenetic D statistic
# statistical inference: 0 = BM model, 1 = no signal
# result: follows a BM model
phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = leaf_cut)

# phylogenetic signal: leaf pulp -----------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = leaf_pulp)

# phylogenetic signal: resin ---------------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = resin)

# phylogenetic signal: mud -----------------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = mud)

# phylogenetic signal: stone ---------------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = stone)

# phylogenetic signal: none ----------------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = none)

# phylogenetic signal: diet breadth --------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = diet)

# phylogenetic signal: voltinism -----------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = volt)

# phylogenetic signal: emergence time -------------------------------------------

phylo.d(data = as.data.frame(trait), 
        phy = tree_nodes_null, 
        names.col = spp, 
        binvar = emer_time2)
