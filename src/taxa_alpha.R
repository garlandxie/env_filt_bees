# Calculate standardized effect size of species richness -----------------------
# Authors of script: Nicholas Sookhan, Garland Xie
# Institutional affiliations: University of Toronto Scarborough
# Contact info: garland.xie@mail.utoronto.ca

# libraries --------------------------------------------------------------------
library(vegan)
library(picante)
library(EcoSimR)
library(here)

# load data --------------------------------------------------------------------

# file-paths
comm_path <- here("data/original", "community_data_matrix.csv")
dist_path <- here("data/working", "dist_spa.csv")

# import 
comm <- read.csv(comm_path, row.names = 1)
dist_spa <- read.csv(dist_path, row.names = 1)

# functions --------------------------------------------------------------------

# gotelli simulation 2 while maintaining species abundances
sim2 <- function(comm){
  
  # create empty data matrix to store randomly generated community
  rand_comm <- matrix(NA,
                      nrow = dim(comm)[1],
                      ncol = dim(comm)[2],
                      dimnames= list(rownames(comm),colnames(comm)))
  
  for (i in 1:ncol(comm)){
    
    # get abundance of species i
    abun <- specnumber(comm, MARGIN = 2)[i]
    
    # sort abundances of species i randomly across sites 
    rand <- list()
    for (j in 1:abun){
      rand[[j]] <- sample(c(1, rep(0,nrow(comm)-1)), 
                          size=nrow(comm), 
                          replace = F)
    }
    
    rand <- do.call("cbind",rand)
    rand <- apply(rand,1,sum)
    
    #assign randomly generated community to rand_comm
    rand_comm[,i] <- rand
  }
  return(rand_comm)
}

# ses for richness - adapt mpd from picante
ses.richness <- function (samp,null.model = "frequency", runs = 999) {
  
  sr.obs <- specnumber(samp)
  
  sr.rand <-  switch(null.model,
                      frequency = t(replicate(runs, specnumber(randomizeMatrix(samp, "frequency")))),
                      richness = t(replicate(runs, specnumber(randomizeMatrix(samp, "richness")))),
                      independentswap = t(replicate(runs, specnumber( randomizeMatrix(samp, "independentswap")))),
                                        sim2 = t(replicate(runs, specnumber( sim2(samp))))
                      
  )
  
  sr.rand.mean <- apply(X = sr.rand, 
                        MARGIN = 2,
                        FUN = mean, 
                        na.rm = TRUE)
  
  sr.rand.sd <- apply(X = sr.rand, 
                      MARGIN = 2, 
                      FUN = sd, 
                      na.rm = TRUE)
  
  sr.obs.z <- (sr.obs - sr.rand.mean)/sr.rand.sd
  
  sr.obs.rank <- apply(X = rbind(sr.obs, sr.rand),
                       MARGIN = 2, 
                       FUN = rank)[1, ]
  
  sr.obs.rank <- ifelse(is.na(sr.rand.mean), NA, sr.obs.rank)
  
  sr.obs.p <- sr.obs.rank/(runs+1)
  
  # two-tailed
  sr.obs.p <- apply(matrix(sr.obs.p), 2,
                     function(sr.obs.p){
                       ifelse(sr.obs.p<0.5, 
                              sr.obs.p*2,
                              (1-sr.obs.p)*2)
                       }
  )
  
  data.frame(sr.obs, 
             sr.rand.mean, 
             sr.rand.sd, 
             sr.obs.rank, 
             sr.obs.z, 
             sr.obs.p, 
             runs = runs,
             row.names = row.names(samp))
}

# calculate metrics ------------------------------------------------------------
ses_richness <- ses.richness(comm,"sim2")

# save and clean ---------------------------------------------------------------

# write to disk
file_path <- here("data/working", "ses_richness.csv")
write.csv(ses_richness, file_path)

