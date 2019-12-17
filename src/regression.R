# libraries --------------------------------------------------------------------
library(here)
library(dplyr)
library(readr)
library(validate)
library(sp)
library(tibble)
library(spdep)
library(ggplot2)

# import -----------------------------------------------------------------------

# relative file-paths
site_path <- here("data/original", "site_info_coords.csv")
mpd_path  <- here("data/working", "ses_mpd.csv") 
mfd_path  <- here("data/working", "ses_mfd.csv")
rich_path <- here("data/working", "ses_richness.csv")
met_250_path <- here("data/original", "land_use_metrics_250.csv")
met_500_path <- here("data/original", "land_use_metrics_500.csv")

# import
site <- read_csv(site_path)
ses_mpd <- read_csv(mpd_path)
ses_mfd <- read_csv(mfd_path)
ses_richness <- read_csv(rich_path)
l_met_250 <- read_csv(met_250_path)
l_met_500 <- read_csv(met_500_path)

# check packaging --------------------------------------------------------------

# site
str(site)
head(site, n = 5)
tail(site, n = 5)

# ses mpd
str(ses_mpd)
head(ses_mpd, n = 5)
tail(ses_mpd, n = 5)

# ses mfd
str(ses_mfd)
head(ses_mfd, n = 5)
tail(ses_mfd, n = 5)

# ses richness
str(ses_richness)
head(ses_richness, n = 5)
tail(ses_richness, n = 5)

# met 250
str(l_met_250)
head(l_met_250, n = 5)
tail(l_met_250, n = 5)

# met 500
str(l_met_500)
head(l_met_500, n = 5)
tail(l_met_500, n = 5)

# data cleaning ----------------------------------------------------------------

df <- site %>%
  full_join(l_met_250, by = "X1") %>%
  full_join(l_met_500, by = "X1") %>%
  full_join(ses_mpd, by = "X1") %>%
  full_join(ses_mfd, by = "X1") %>%
  full_join(ses_richness, by = "X1") %>%
  select(ID = "X1",
         mtm3deg_nor,
         mtm3deg_eas,
         grass_250_percent = "prop.landscape_250_grass",
         tree_250_percent   = "prop.landscape_250_tree_canopy",
         urban_250_percent  = "prop.landscape_250_urban",
         grass_500_percent  = "prop.landscape_500_grass",
         tree_500_percent   = "prop.landscape_500_tree_canopy",
         urban_500_percent  = "prop.landscape_500_urban",
         ses_mfd = "mfd.obs.z",
         ses_mpd = "mpd.obs.z",
         ses_R = "sr.obs.z")

# linear models: model fitting -------------------------------------------------

# ses richness
lm_250_sr <- df %>%
  filter(!is.na(urban_250_percent)) %>%
  lm(ses_R ~ urban_250_percent, data = .) 

lm_500_sr <- df %>%
  filter(!is.na(urban_500_percent)) %>%
  lm(ses_R ~ urban_500_percent, data = .)

# ses mpd
lm_250_mpd <- df %>%
  filter(!is.na(urban_250_percent)) %>%
  lm(ses_mpd ~ urban_250_percent, data = .)

lm_500_mpd <- df %>%
  filter(!is.na(urban_500_percent)) %>%
  lm(ses_mpd ~ urban_500_percent, data = .)

# ses mfd
lm_250_mfd <- df %>%
  filter(!is.na(urban_250_percent)) %>%
  lm(ses_mfd ~ urban_250_percent, data = .)

lm_500_mfd <- df %>%
  filter(!is.na(urban_500_percent)) %>%
  lm(ses_mfd ~ urban_500_percent, data = .)

# linear models: model summary -------------------------------------------------

# sr
summary(lm_250_sr) 
summary(lm_500_sr) 

# mpd
summary(lm_250_mpd) 
summary(lm_500_mpd)

# mfd
summary(lm_250_mfd) 
summary(lm_500_mfd)

# linear models: model adequacy ------------------------------------------------

# sr
plot(lm_250_sr) 
plot(lm_500_sr)

hist(resid(lm_250_sr))
hist(resid(lm_500_sr))

# mpd
plot(lm_250_mpd) 
plot(lm_500_mpd)

hist(resid(lm_250_mpd))
hist(resid(lm_500_mpd))

# mfd
plot(lm_250_mfd) 
plot(lm_500_mfd)

hist(resid(lm_250_mfd))
hist(resid(lm_500_mfd))

# linear models: spatial autocorrelation ---------------------------------------

# 250 m spatial scale
coords_250 <- df %>%
  filter(!is.na(urban_250_percent)) %>%
  select(ID, mtm3deg_nor, mtm3deg_eas) %>%
  column_to_rownames(var = "ID")

coordinates(coords_250) = ~mtm3deg_eas + mtm3deg_nor
proj4string(coords_250) <- CRS("+proj=tmerc +lat_0=0 +lon_0=-79.5 \n 
                               +k=0.9999 +x_0=304800 +y_0=0 +datum=NAD27 \n 
                               +units=m +no_defs +ellps=clrk66 \n
                               +nadgrids=@conus,@alaska,@ntv2_0.gsb, \n
                               @ntv1_can.dat")

# 500m spatial scale 
coords_500 <- df %>%
  filter(!is.na(urban_500_percent)) %>%
  select(ID, mtm3deg_nor, mtm3deg_eas) %>%
  column_to_rownames(var = "ID")

coordinates(coords_500) = ~mtm3deg_eas + mtm3deg_nor
proj4string(coords_500) <- CRS("+proj=tmerc +lat_0=0 +lon_0=-79.5 \n 
                               +k=0.9999 +x_0=304800 +y_0=0 +datum=NAD27 \n 
                               +units=m +no_defs +ellps=clrk66 \n
                               +nadgrids=@conus,@alaska,@ntv2_0.gsb, \n
                               @ntv1_can.dat")

# moran i test for spatial autocorrelation in residuals

lm_morantest_250 <- partial(
  lm.morantest, 
  listw = nb2listw(knn2nb(knearneigh(coords_250, 8)),style = "W"),
)

lm_morantest_500 <- partial(
  lm.morantest, 
  listw = nb2listw(knn2nb(knearneigh(coords_500, 8)),style = "W"),
)

# sr
lm_morantest_250(lm_250_sr)
lm_morantest_500(model = lm_500_sr)

# mpd
lm_morantest_250(model = lm_250_mpd)
lm_morantest_500(model = lm_500_mpd)

# mfd
lm_morantest_250(lm_250_mfd)
lm_morantest_500(lm_500_mfd)

# plots ------------------------------------------------------------------------

plot_lm_R_250 <- df %>%
  filter(!is.na(urban_250_percent) & !is.na(ses_R)) %>%
  ggplot(aes(x = urban_250_percent, y = ses_R)) + 
  geom_jitter() + 
  labs(y = "ses.Richness", 
       x = "% Impervious Cover (250m scale)") + 
  theme_minimal()

plot_lm_R_500 <- df %>%
  filter(!is.na(urban_500_percent) & !is.na(ses_R)) %>%
  ggplot(aes(x = urban_500_percent, y = ses_R)) + 
  geom_jitter() + 
  labs(y = "ses.Richness", 
       x = "% Impervious Cover (500m scale)") + 
  theme_minimal()

plot_lm_mpd_250 <- df %>%
  filter(!is.na(urban_250_percent) & !is.na(ses_mpd)) %>%
  ggplot(aes(x = urban_250_percent, y = ses_mpd)) + 
  geom_smooth(method = "lm") + 
  geom_jitter() + 
  labs(y = "ses.MPD", 
       x = "% Impervious Cover (250m scale)") + 
  theme_minimal()

plot_lm_mpd_500 <- df %>%
  filter(!is.na(urban_500_percent) & !is.na(ses_mpd)) %>%
  ggplot(aes(x = urban_500_percent, y = ses_mpd)) + 
  geom_smooth(method = "lm") + 
  geom_jitter() + 
  labs(y = "ses.MPD", 
       x = "% Impervious Cover (500m scale)") + 
  theme_minimal()

plot_lm_mfd_250 <- df %>%
  filter(!is.na(urban_250_percent) & !is.na(ses_mfd)) %>%
  ggplot(aes(x = urban_250_percent, y = ses_mfd)) + 
  geom_jitter() + 
  labs(y = "ses.MFD", 
       x = "% Impervious Cover (250m scale)") + 
  theme_minimal()

plot_lm_mfd_500 <- df %>%
  filter(!is.na(urban_500_percent) & !is.na(ses_mfd)) %>%
  ggplot(aes(x = urban_500_percent, y = ses_mfd)) + 
  geom_smooth(method = "lm") + 
  geom_jitter() + 
  labs(y = "ses.MFD", 
       x = "% Impervious Cover (500m scale)") + 
  theme_minimal()