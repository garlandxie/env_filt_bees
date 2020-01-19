# libraries -------------------------------------------------------------------
library(here)  # for creating relative file-paths
library(dplyr)
library(vegan)
library(ggplot2)
# import ----------------------------------------------------------------------

# file-paths
comm_path <- here("data/original", "community_data_matrix.csv")

# read csv
comm_df <- read.csv(comm_path, row.names = 1)

# check packaging --------------------------------------------------------------

str(comm_df)
head(comm_df, n = 5)
tail(comm_df, n = 5)

# clean data -------------------------------------------------------------------

# abundance cut-off: plots with more than 3 individuals
comm_3Ab <- comm_df[rowSums(comm_df) >= 3, ]

# abundance cut-off: plots with more than 5 individuals
comm_5Ab <- comm_df[rowSums(comm_df) >= 5, ]

# abundance cut-off: plots with more than 7 individuals
comm_7Ab <- comm_df[rowSums(comm_df) >= 7, ]

# plot 

plot_comm_3Ab <- data.frame(
  SR = rowSums(decostand(comm_3Ab, method = "pa")),
  Ab = rowSums(comm_3Ab)) %>%
  ggplot(aes(x = Ab, y = SR)) + 
  geom_point()

plot_comm_5Ab <- data.frame(
  SR = rowSums(decostand(comm_5Ab, method = "pa")),
  Ab = rowSums(comm_5Ab)) %>%
  ggplot(aes(x = Ab, y = SR)) + 
  geom_point()

plot_comm_5Ab <- data.frame(
  SR = rowSums(decostand(comm_7Ab, method = "pa")),
  Ab = rowSums(comm_7Ab)) %>%
  ggplot(aes(x = Ab, y = SR)) + 
  geom_point()
  
# save to disk -----------------------------------------------------------------

saveRDS(comm_3Ab, file = here("data/working", "comm_3Ab.RDS"))
saveRDS(comm_5Ab, file = here("data/working", "comm_5Ab.RDS"))
saveRDS(comm_7Ab, file = here("data/working", "comm_7Ab.RDS"))
