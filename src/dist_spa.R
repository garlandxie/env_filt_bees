# libraries --------------------------------------------------------------------
library(sp)
library(here)
library(validate)
library(dplyr)
library(reshape2)
library(ggplot2)

# load data --------------------------------------------------------------------

# relative file-paths
site_path <- here("data/original", "site_info_coords.csv")

# read csv: avoid strings as factors
site <- read.csv(site_path, 
                 stringsAsFactors = FALSE,
                 row.names = 1)

# check packaging --------------------------------------------------------------

# check for data-types, header and footer
str(site)
head(site, n = 5)
tail(site, n = 5)

# validate data
# note: all sites are in the GTA = narrow coordinate range
site %>%
  check_that(
    is.character(site),
    is.character(name),
    wgs84_lat >= 43, 
    wgs84_long <= -79,
    nad27_lat >= 43,
    nad27_long <= -79,
    mtm3deg_nor > 0,
    mtm3deg_eas > 0) %>%
  summary()

# visualize coordinates --------------------------------------------------------

# site distributions looks like Toronto, ON
site %>%
  ggplot(aes(x = wgs84_long, y = wgs84_lat)) + 
  geom_point() + 
  labs(x = "longitude",
       y = "latitude") + 
  theme_minimal()

# geospatial analysis # --------------------------------------------------------

# new df for spatial projections
site_sp <- site

# coerce site database to SpatialPointsDataFrame obj
coordinates(site_sp) = ~mtm3deg_nor + mtm3deg_eas
proj4string(site_sp) <- CRS("+proj=tmerc +lat_0=0 +lon_0=-79.5 +k=0.9999 +x_0=304800 +y_0=0 +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat")

# calculate spatial distance matrix --------------------------------------------

# create an empty n x n distance matrix
# n = number of sites
# row and column names are site names

dist_spa <- matrix(0, 
                   nrow = nrow(site_sp),
                   ncol = nrow(site_sp),
                   dimnames = list(row.names(site_sp), row.names(site_sp)) 
                   )

# apply spatial distances within the appropriate cells
for (i in 1:nrow(site_sp)) { 
  dist_spa[,i] <- spDistsN1(site_sp[1:nrow(site),1] ,site_sp[i,1])
}

# double checks
all(dist_spa >= 0)
isSymmetric.matrix(dist_spa)

# save to disk -----------------------------------------------------------------

dist_path <- here("data/working", "dist_spa.csv")
write.csv(dist_spa,"dist_spa.csv")

