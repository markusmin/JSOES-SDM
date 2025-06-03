# Code to generate SDMs for species caught in the JSOES trawl

## load libraries
library(tidyverse)
library(here)
library(readxl)
library(janitor)
library(sf)
library(sdmTMB)
library(TMB)
library(rnaturalearth)
library(sdmTMBextra)
library(fmesher)
library(lubridate)
library(geosphere)
library(geos)

#### load and reformat trawl data ####
jsoes <- clean_names(read_excel(here::here("Data", "Markus_Min_Trawl_CTD_Chl_Nuts_Thorson_Scheuerell_5.2.24_FINAL.xlsx")))

species_column_names <- colnames(jsoes)[32:ncol(jsoes)]

# change nmi to km to keep units consistent
jsoes %>% 
  mutate(km_from_shore = nmi_from_shore * 1.852) -> jsoes

jsoes %>% 
  pivot_longer(., cols = all_of(species_column_names), values_to = "n_per_km", names_to = "species") %>% 
  mutate(n = n_per_km*trawl_dist_km) -> jsoes_long

# for now, let's drop any data that's missing lat/longs
jsoes_long <- filter(jsoes_long, !(is.na(mid_lat)))

# Drop any samples that don't have temperature data (we could interpolate both of these later)
jsoes_long <- filter(jsoes_long, !(is.na(x3m_temp_c)))



# Add UTM coordinates to get an equal distance projection and get lat/longs
utm_crs <- get_crs(jsoes_long, c("mid_long", "mid_lat")) # check that

jsoes_long <- add_utm_columns(
  jsoes_long,
  ll_names = c("mid_long", "mid_lat"),
  ll_crs = 4326,
  utm_names = c("Lon.km", "Lat.km"),
  utm_crs = utm_crs,
  units = c("km")
)


# Change Lon.km and Lat.km columns to X and Y to avoid problems later
jsoes_long %>% 
  dplyr::rename(X = Lon.km, Y = Lat.km) -> jsoes_long

#### load and reformat spatial data ####
usa_spdf <- st_read("/Users/markusmin/Documents/ESA_RF_2021/map_files/USA_adm0.shp")
# load BC
CAN_spdf <- st_read("/Users/markusmin/Documents/ESA_RF_2021/map_files/canada/lpr_000b16a_e.shp")
BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
BC_proj <- st_transform(BC_spdf, crs = 4326)


# crop them to our desired area
US_west_coast <- sf::st_crop(usa_spdf,
                             c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))

BC_coast <- sf::st_crop(BC_proj,
                        c(xmin = -126, ymin = 44, xmax = -123, ymax = 48.5))



# convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
UTM_zone_10_crs <- 32610

US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)

# make this projection into kilometers
US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)



#### load and reformat SST data ####
SST1 <- read.csv(here::here("Data", "noaacrwsstDaily_9bdb_17ce_0aa0.csv"), skip = 1)
SST2 <- read.csv(here::here("Data", "noaacrwsstDaily_36e3_7cb3_1109.csv"), skip = 1)
SST3 <- read.csv(here::here("Data", "noaacrwsstDaily_519c_f955_0fea.csv"), skip = 1)
SST4 <- read.csv(here::here("Data", "noaacrwsstDaily_727c_5c10_e6a0.csv"), skip = 1)
SST5 <- read.csv(here::here("Data", "noaacrwsstDaily_1176_fc5f_53d2.csv"), skip = 1)
SST6 <- read.csv(here::here("Data", "noaacrwsstDaily_ae0c_9197_083f.csv"), skip = 1)
SST7 <- read.csv(here::here("Data", "noaacrwsstDaily_fe4f_9c9c_442e.csv"), skip = 1)

SST1 %>% 
  bind_rows(., SST2, SST3, SST4, SST5, SST6, SST7) -> SST

SST %>% 
  dplyr::rename(time = UTC, latitude = degrees_north, longitude = degrees_east, SST = degree_C) -> SST
SST %>% 
  mutate(time = ymd_hms(time)) %>% 
  mutate(day_of_month = mday(time)) %>% 
  mutate(year = year(time)) %>% 
  mutate(julian = yday(time))-> SST


# Looks good to me! Now just need to resolve CRS conflicts

# Convert SST to UTM zone 10
# st_as_sf(SST, coords = c("latitude", "longitude"), crs = "4326") -> SST_sf
st_as_sf(SST, coords = c("longitude", "latitude"), crs = 4326) -> SST_sf
US_west_coast

sf::st_transform(SST_sf, crs = UTM_zone_10_crs) -> SST_sf_proj

# SST_plots <- survey_area_basemap +
#   geom_sf(data = SST_sf_proj, aes(color = SST)) +
#   scale_color_viridis_c(option = "A")+
#   facet_wrap(~year, nrow = 3) +
#   theme(legend.position = c(0.95, 0.1))

#### create base map for visualizing data ####
survey_area_basemap_km <- ggplot(US_west_coast_proj_km) +
  geom_sf() +
  geom_sf(data = BC_coast_proj_km) +
  # coord_sf(ylim = c(4922.052, 5342.052), xlim = c(334.8638, 404.8638)) +
  # coord_sf(ylim = c(44,48.5),  xlim = c(-126, -123)) +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks = element_blank(),
        axis.text = element_blank())
# temporary fix for lat/long on axes to just get rid of them - start here to actually fix: https://forum.posit.co/t/converting-axes-to-lat-lon/27181

survey_area_basemap_km

survey_area_basemap <- ggplot(US_west_coast) +
  geom_sf() +
  geom_sf(data = BC_coast) +
  coord_sf(ylim = c(44,48.5),  xlim = c(-126, -123)) +
  scale_x_continuous(breaks = c(124,125,126)) +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

#### Create a survey prediction grid ####

# Create a survey grid from scratch
# Make sure that it aligns with our SST

# only keep one year, since that's all we need to make a grid
SST_sf_proj %>% 
  filter(year == 2000) -> SST_sf_proj_2000

st_as_sf(SST_sf_proj_2000$geometry/1000, crs = UTM_zone_10_crs) -> SST_sf_proj_km_2000

SST_sf_proj_km_2000$dist_shore <- as.numeric(st_distance(SST_sf_proj_km_2000, US_west_coast_proj_km))


# 100 km offshore is the max, then 20 km north/south of limits of survey region
# The vast majority of stations are within 60 km of shore, but there are just a couple that are further out in the Columbia River Plume up to 83 km from shore
# We also want to exclude values that are within 0.5 km of shore, because these
# don't have measureable SST values (they're basically on shore)
hist(jsoes_long$km_from_shore)

SST_sf_proj_km_2000 %>% 
  filter(dist_shore <= 100 & dist_shore >= 0.5) -> grid_within_100km_shore


# let's inspect the grid that we created
ggplot(grid_within_100km_shore) +
  geom_sf()

# looks okay, but there are still sections that we need to crop out below

min(jsoes_long$Y) -> min_Y
max(jsoes_long$Y) -> max_Y


survey_domain <- st_crop(grid_within_100km_shore, xmin = 0, xmax = 100000, ymin=min_Y-20, ymax=max_Y+20)

st_difference(survey_domain, US_west_coast_proj_km) -> survey_domain

# inspect survey domain - still needs things cropped
ggplot(survey_domain) +
  geom_sf()

# manually remove sections for strait, hood canal (and other inland waters), Grays Harbor, Willapa Bay, Columbia River estuary
strait <- st_as_sfc(st_bbox(c(xmin=-124.65, xmax=-122, ymin=47.9375, ymax=49), crs = "WGS84"))
strait_proj <- sf::st_transform(strait, crs = UTM_zone_10_crs)
strait_proj_km <- st_as_sf(strait_proj/1000, crs = UTM_zone_10_crs)

inland_waters <- st_as_sfc(st_bbox(c(xmin=-123.5, xmax=-120, ymin=44, ymax=49.5), crs = "WGS84"))
inland_waters_proj <- sf::st_transform(inland_waters, crs = UTM_zone_10_crs)
inland_waters_proj_km <- st_as_sf(inland_waters_proj/1000, crs = UTM_zone_10_crs)

grays_harbor <- st_as_sfc(st_bbox(c(xmin=-124.15, xmax=-123.6, ymin=46.83, ymax=47.09), crs = "WGS84"))
grays_harbor_proj <- sf::st_transform(grays_harbor, crs = UTM_zone_10_crs)
grays_harbor_proj_km <- st_as_sf(grays_harbor_proj/1000, crs = UTM_zone_10_crs)

willapa_bay <- st_as_sfc(st_bbox(c(xmin=-124.06, xmax=-123.5, ymin=46.34, ymax=46.8), crs = "WGS84"))
willapa_bay_proj <- sf::st_transform(willapa_bay, crs = UTM_zone_10_crs)
willapa_bay_proj_km <- st_as_sf(willapa_bay_proj/1000, crs = UTM_zone_10_crs)

estuary <- st_as_sfc(st_bbox(c(xmin=-124.02, xmax=-123.5, ymin=46.12, ymax=46.35), crs = "WGS84"))
estuary_proj <- sf::st_transform(estuary, crs = UTM_zone_10_crs)
estuary_proj_km <- st_as_sf(estuary_proj/1000, crs = UTM_zone_10_crs)

survey_domain %>% 
  st_difference(., strait_proj_km) %>% 
  st_difference(., inland_waters_proj_km) %>% 
  st_difference(., grays_harbor_proj_km) %>% 
  st_difference(., willapa_bay_proj_km) %>% 
  st_difference(., estuary_proj_km) -> survey_domain

# Create a concave hull around points
st_concave_hull(st_union(survey_domain), ratio = 0.1) -> survey_domain_polygon

ggplot(survey_domain_polygon) +
  geom_sf()

ggplot(survey_domain) +
  geom_sf()

# let's inspect our survey domain - looks pretty good!
survey_area_basemap_km +
  geom_sf(data = survey_domain_polygon, fill = "blue", alpha = 0.2)

#### Create a prediction grid ####
# Take SST grid and intersect with survey domain
SST_sf_proj %>% 
  mutate(geometry = geometry/1000) %>% 
  st_set_crs(UTM_zone_10_crs) %>% 
  st_intersection(survey_domain_polygon) %>% 
  dplyr::select(SST, year, geometry)-> SST_survey_domain

# add in distance from shore
SST_survey_domain %>% 
  mutate(sf_dist_shore = as.numeric(st_distance(geometry, US_west_coast_proj_km))) %>% 
  # I'm not sure why we have to filter again here, but make sure to exclude points that are 
  # > 100 km from shore and < 0.5 km from shore - maybe some rounding errors?
  filter(sf_dist_shore <= 100 & sf_dist_shore >= 0.5) -> survey_domain_cov

# plot it - looks good!
ggplot(survey_domain_cov) + geom_sf()

# Now let's reformat it to be in the form of a data frame that we can use to predict
as.data.frame(st_coordinates(survey_domain_cov)) -> survey_domain_cov_coords
survey_predict_grid <- data.frame(X = survey_domain_cov_coords$X,
                                  Y = survey_domain_cov_coords$Y,
                                  SST = survey_domain_cov$SST,
                                  year = survey_domain_cov$year,
                                  dist_shore = as.numeric(survey_domain_cov$sf_dist_shore))

# rescale prediction grid
survey_predict_grid %>% 
  mutate(SST_scaled = as.numeric(scale(SST)),
         dist_shore_scaled = as.numeric(scale(dist_shore))) -> survey_predict_grid


#### Create the mesh for SDM ####

# get the sampling coordinates
jsoes_long %>% 
  distinct(year, station_code, .keep_all = TRUE) %>% 
  dplyr::select(-species) -> jsoes_samples

# The object jsoes_samples contains the each sample collected by the survey,
# along with associated lat/lon information (plotted below)
ggplot(jsoes_samples, aes(x = X, y = Y)) +
  geom_point()

# From this object, we then can create a mesh that will be used by the SPDE method.

# Create mesh for SPDE
inla_mesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(jsoes_samples$X, jsoes_samples$Y), # coordinates
  cutoff = 10 # minimum triangle edge length
)

plot(inla_mesh)

spde <- fm_fem( inla_mesh, 
                order = 2 )

# create projection matrix from vertices to samples
A_is = fm_evaluator(inla_mesh, loc = as.matrix(data.frame(X = jsoes_samples$X, Y = jsoes_samples$Y)))$proj$A

# create projection matrix from vertices to prediction grid
# the grid dimensions are the same each year, so just select one year and drop fields besides geometry
survey_domain_cov %>% 
  filter(year == 1998) %>% 
  dplyr::select(geometry) -> survey_domain_cov_grid

A_gs = fm_evaluator( inla_mesh, loc=st_coordinates(survey_domain_cov_grid))$proj$A
# A_gs = fm_evaluator( inla_mesh, loc=st_coordinates(st_centroid(survey_grid_df)))$proj$A

# Make the sdmTMB mesh
# for consistency, use the inla_mesh for all
mesh_sdmTMB <- make_mesh(jsoes_samples, xy_cols = c("X", "Y"), mesh = inla_mesh)

# confirm that the meshes used for bespoke TMB and sdmTMB are the same
plot(inla_mesh)
plot(mesh_sdmTMB)

# yes, mesh is the same

png(here::here("figures", "SDMs", "jsoes_mesh.png"), width=4, height=6, res=200, units="in")
plot(mesh_sdmTMB)
dev.off()

#### SDM model selection process ####

# We will fit SDMs for Interior Spring chinook
csyis <- subset(jsoes_long, species == "chinook_salmon_yearling_interior_sp")

#### Model 1: Null model ####

compile(here::here("analysis", "jsoes_null_sdm.cpp"))
dyn.load( dynlib(here::here("analysis", "jsoes_null_sdm")))

Data = list( "D_i"= csyis$n_per_km, "t_i" = csyis$year - min(csyis$year), 
             "n_t" = length(unique(csyis$year)),
             "A_is"=A_is, "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta_t"=rep(0, length(unique(csyis$year))), "ln_tau"=0, "ln_kappa"=0, 
               "ln_phi" = 0, "finv_power" = 0, "omega_s"=rnorm(nrow(spde$c0)))
null_Obj = MakeADFun( data=Data, parameters=Params, random= c("omega_s"),
                      DLL = "jsoes_null_sdm")

# Optimize
null_Opt = nlminb( start=null_Obj$par, obj=null_Obj$fn, grad=null_Obj$gr )
null_Opt$SD = sdreport( null_Obj, bias.correct=TRUE )
null_report = null_Obj$report()


# Fit sdmTMB with effect of year and effect of space, no interaction (model 1, null model)
jsoes_null_sdmTMB <- sdmTMB(n_per_km ~ 0 + as.factor(year),
                            data = csyis,
                            mesh = mesh_sdmTMB,
                            time = "year",
                            spatial = "on",
                            spatiotemporal = "off",
                            family = tweedie(link = "log"))

#### Model 2: Spatiotemporal effect ####

compile(here::here("analysis", "jsoes_st_sdm.cpp"))
dyn.load( dynlib(here::here("analysis", "jsoes_st_sdm")))

Data = list( "D_i"= csyis$n_per_km, "t_i" = csyis$year - min(csyis$year), 
             "n_t" = length(unique(csyis$year)),
             "A_is"=A_is, "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta_t"=rep(0, length(unique(csyis$year))), 
               "ln_tau_omega"=0, "ln_tau_epsilon"=0,
               "ln_kappa"=0, "ln_phi" = 0, 
               "finv_power" = 0, "omega_s"=rnorm(nrow(spde$c0)),
               "epsilon_st"=matrix(0, nrow=nrow(spde$c0), ncol=Data$n_t))

st_Obj = MakeADFun( data=Data, parameters=Params, random= c("omega_s", "epsilon_st"),
                    DLL = "jsoes_st_sdm")

# Optimize
st_Opt = nlminb( start=st_Obj$par, obj=st_Obj$fn, grad=st_Obj$gr )
# Opt$SD = sdreport( Obj, bias.correct=TRUE )
st_report = st_Obj$report()
# Fit sdmTMB with effect of year and effect of space and spatio-temporal effect
jsoes_st_sdmTMB <- sdmTMB(n_per_km ~ 0 + as.factor(year),
                          data = csyis,
                          mesh = mesh_sdmTMB,
                          time = "year",
                          spatial = "on",
                          spatiotemporal = "iid",
                          family = tweedie(link = "log"))

# check that they're the same
st_report$beta_t
jsoes_st_sdmTMB

#### Model 3: Spatiotemporal effect as AR1 ####

compile(here::here("analysis", "jsoes_stac_sdm.cpp"))
dyn.load( dynlib(here::here("analysis", "jsoes_stac_sdm")))

Data = list( "D_i"= csyis$n_per_km, "t_i" = csyis$year - min(csyis$year), 
             "n_t" = length(unique(csyis$year)),
             "A_is"=A_is, "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta_t"=rep(0, length(unique(csyis$year))), 
               "ln_tau_omega"=0, "ln_tau_epsilon"=0,
               "ln_kappa"=0, "logit_rhoE" = 0,  "ln_phi" = 0, 
               "finv_power" = 0, "omega_s"=rnorm(nrow(spde$c0)),
               "epsilon_st"=matrix(0, nrow=nrow(spde$c0), ncol=Data$n_t))

stac_Obj = MakeADFun( data=Data, parameters=Params, random= c("omega_s", "epsilon_st"),
                      DLL = "jsoes_stac_sdm")

# Optimize
stac_Opt = nlminb( start=stac_Obj$par, obj=stac_Obj$fn, grad=stac_Obj$gr )
# Opt$SD = sdreport( Obj, bias.correct=TRUE )
stac_report = stac_Obj$report()

# Fit sdmTMB with effect of year and effect of space and spatio-temporal effect
jsoes_stac_sdmTMB <- sdmTMB(n_per_km ~ 0 + as.factor(year),
                            data = csyis,
                            mesh = mesh_sdmTMB,
                            time = "year",
                            spatial = "on",
                            spatiotemporal = "ar1",
                            family = tweedie(link = "log"))


# compare with bespoke TMB implementation
jsoes_stac_sdmTMB
stac_report$beta_t
# exactly the same

#### Model 4: intercept as temporal effect, spatiotemporal as IID (Same as model 2, but now with AR1 intercept) ####

compile(here::here("analysis",  "jsoes_tempac_sdm.cpp"))
dyn.load( dynlib(here::here("analysis", "jsoes_tempac_sdm")))

Data = list( "D_i"= csyis$n_per_km, "t_i" = csyis$year - min(csyis$year), 
             "n_t" = length(unique(csyis$year)),
             "A_is"=A_is, "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta_t"=rep(0, length(unique(csyis$year))), 
               "ln_tau_omega"=0, "ln_tau_epsilon"=0,
               "ln_kappa"=0, "logit_rhoB" = 0,  "ln_phi" = 0, 
               "finv_power" = 0, "ln_sigmaB" = 0, "omega_s"=rnorm(nrow(spde$c0)),
               "epsilon_st"=matrix(0, nrow=nrow(spde$c0), ncol=Data$n_t))

tempac_Obj = MakeADFun( data=Data, parameters=Params, random= c("omega_s", "epsilon_st", "beta_t"),
                        DLL = "jsoes_tempac_sdm")

# Optimize
tempac_Opt = nlminb( start=tempac_Obj$par, obj=tempac_Obj$fn, grad=tempac_Obj$gr )
# Opt$SD = sdreport( tempac_Obj, bias.correct=TRUE )
tempac_report = tempac_Obj$report()

# Fit sdmTMB with effect of year and effect of space and spatio-temporal effect
# model 3:
# jsoes_stac_sdmTMB <- sdmTMB(n_per_km ~ 0 + as.factor(year),
#                             data = csyis,
#                             mesh = mesh_sdmTMB,
#                             time = "year",
#                             spatial = "on",
#                             spatiotemporal = "ar1",
#                             family = tweedie(link = "log"))

# model 4:
# jsoes_tempac_sdmTMB <- sdmTMB(n_per_km ~ 0 + as.factor(year),
#                               data = csyis,
#                               mesh = mesh_sdmTMB,
#                               time = "year",
#                               time_varying = ~ 0 + as.factor(year),
#                               time_varying_type = "ar1",
#                               spatial = "on",
#                               spatiotemporal = "iid",
#                               family = tweedie(link = "log"),
#                               silent = FALSE)

# This doesn't work - it runs the optimizer but then gives this error message:
# Error in solve.default(h, g) : system is computationally singular: reciprocal condition number = 1.84191e-29

# what if we only include year in the time-varying part but leave it out of the main formula?
# jsoes_tempac_sdmTMB <- sdmTMB(n_per_km ~ 0,
#                               data = csyis,
#                               mesh = mesh_sdmTMB,
#                               time = "year",
#                               time_varying = ~ 0 + as.factor(year),
#                               time_varying_type = "ar1",
#                               spatial = "on",
#                               spatiotemporal = "iid",
#                               family = tweedie(link = "log"),
#                               silent = FALSE)
# this also doesn't work - Error in solve.default(h, g) : 
# system is computationally singular: reciprocal condition number = 2.18518e-20

# what if we add an intercept?
# jsoes_tempac_sdmTMB <- sdmTMB(n_per_km ~ 1,
#                               data = csyis,
#                               mesh = mesh_sdmTMB,
#                               time = "year",
#                               time_varying = ~ 0 + as.factor(year),
#                               time_varying_type = "ar1",
#                               spatial = "on",
#                               spatiotemporal = "iid",
#                               family = tweedie(link = "log"),
#                               silent = FALSE)
# nope - Error in solve.default(h, g) : 
# system is computationally singular: reciprocal condition number = 1.79181e-22

# What if we change spatiotemporal to ar1 as well?
# jsoes_tempac_sdmTMB <- sdmTMB(n_per_km ~ 1,
#                               data = csyis,
#                               mesh = mesh_sdmTMB,
#                               time = "year",
#                               time_varying = ~ 0 + as.factor(year),
#                               time_varying_type = "ar1",
#                               spatial = "on",
#                               spatiotemporal = "ar1",
#                               family = tweedie(link = "log"),
#                               silent = FALSE)
# Error in solve.default(h, g) : 
# system is computationally singular: reciprocal condition number = 3.60052e-19

# "time_varying takes a one-sided formula. ~ 1 implies a time-varying intercept."
# Okay, let's try getting rid of as.factor(year)
jsoes_tempar1_sdmTMB <- sdmTMB(n_per_km ~ 0,
                              data = csyis,
                              mesh = mesh_sdmTMB,
                              time = "year",
                              time_varying = ~ 1,
                              time_varying_type = "ar1",
                              spatial = "on",
                              spatiotemporal = "iid",
                              family = tweedie(link = "log"),
                              silent = FALSE)
# this works - but does it give us what we think it should?
# Well, we now have intercept terms. Perhaps before it was singular
# because it was trying to estimate a time-varying intercept for year, and also
# a separate effect from the year column? That would lead to perfect collinearity...

# Ok - let's try and get the TMB code for an AR1 process the same
compile(here::here("analysis", "jsoes_tempar1_sdm.cpp"))
dyn.load( dynlib(here::here("analysis", "jsoes_tempar1_sdm")))

Data = list( "D_i"= csyis$n_per_km, "t_i" = csyis$year - min(csyis$year), 
             "n_t" = length(unique(csyis$year)),
             "A_is"=A_is, "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta_t"=rep(0, length(unique(csyis$year))), 
               "ln_tau_omega"=0, "ln_tau_epsilon"=0,
               "ln_kappa"=0, "logit_rhoB" = 0,  "ln_phi" = 0, 
               "finv_power" = 0, "ln_sigmaB" = 0, "omega_s"=rnorm(nrow(spde$c0)),
               "epsilon_st"=matrix(0, nrow=nrow(spde$c0), ncol=Data$n_t))

tempar1_Obj = MakeADFun( data=Data, parameters=Params, random= c("omega_s", "epsilon_st", "beta_t"),
                        DLL = "jsoes_tempar1_sdm")

# Optimize
tempar1_Opt = nlminb( start=tempar1_Obj$par, obj=tempar1_Obj$fn, grad=tempar1_Obj$gr )
# Opt$SD = sdreport( tempar1_Obj, bias.correct=TRUE )
tempar1_report = tempar1_Obj$report()

# compare with sdmTMB
jsoes_tempar1_sdmTMB
tempar1_report$beta_t

# they're the same! wooooooo
# the only thing I don't understand is the term for the variance in the AR1 process.
# not sure why it's different here vs. in Jim's code
# We are going to trust the sdmTMB formulation for an AR1, so keep the additional
# correction factor for the variance

#### Model 5: Intercept as temporal (AR1), spatiotemporal as AR1 ####

compile(here::here("analysis", "jsoes_star1_tempar1_sdm.cpp"))
dyn.load( dynlib(here::here("analysis", "jsoes_star1_tempar1_sdm")))

Data = list( "D_i"= csyis$n_per_km, "t_i" = csyis$year - min(csyis$year), 
             "n_t" = length(unique(csyis$year)),
             "A_is"=A_is, "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta_t"=rep(0, length(unique(csyis$year))), 
               "ln_tau_omega"=0, "ln_tau_epsilon"=0,
               "ln_kappa"=0, "logit_rhoB" = 0, "logit_rhoE" = 0, "ln_phi" = 0, 
               "finv_power" = 0, "ln_sigmaB" = 0, "omega_s"=rnorm(nrow(spde$c0)),
               "epsilon_st"=matrix(0, nrow=nrow(spde$c0), ncol=Data$n_t))

star1_tempar1_Obj = MakeADFun( data=Data, parameters=Params, random= c("omega_s", "epsilon_st", "beta_t"),
                             DLL = "jsoes_star1_tempar1_sdm")

# Optimize
star1_tempar1_Opt = nlminb( start=star1_tempar1_Obj$par, obj=star1_tempar1_Obj$fn, grad=star1_tempar1_Obj$gr )
# Opt$SD = sdreport( star1_tempar1_Obj, bias.correct=TRUE )
star1_tempar1_report = star1_tempar1_Obj$report()

# Fit sdmTMB model
jsoes_star1_tempar1_sdmTMB <- sdmTMB(n_per_km ~ 0,
                                   data = csyis,
                                   mesh = mesh_sdmTMB,
                                   time = "year",
                                   time_varying = ~ 1,
                                   time_varying_type = "ar1",
                                   spatial = "on",
                                   spatiotemporal = "ar1",
                                   family = tweedie(link = "log"),
                                   silent = FALSE)

# compare with bespoke sdmTMB implementation
jsoes_star1_tempar1_sdmTMB
star1_tempar1_report$beta_t
# once again slightly different

#### Models without covariates: cross-validation ####

# Ok - so with cross-validation, because it's a random process, it'll vary
# quite a bit between runs. If you change the seed, you'll change the result.
# In our case, the change is quite drastic.

# Set up for parallelization
library(future)
plan(multisession, workers = 4)

# Model 1: Null model

# Let's compare variability vs. number of folds
jsoes_null_cv_comp <- data.frame(model_run = 1:45,
                                 n_folds = rep(c(4, 7, 10), each = 15),
                                 sum_loglik = rep(NA, 45))

set.seed(123)
for (i in 1:nrow(jsoes_null_cv_comp)){
  # set.seed(jsoes_null_cv_comp$model_run[i] * jsoes_null_cv_comp$n_folds[i])
  jsoes_null_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year),
                                    data = csyis,
                                    mesh = mesh_sdmTMB,
                                    time = "year",
                                    spatial = "on",
                                    spatiotemporal = "off",
                                    family = tweedie(link = "log"),
                                    k_folds = jsoes_null_cv_comp$n_folds[i])
  jsoes_null_cv_comp$model_run[i] <- i
  jsoes_null_cv_comp$sum_loglik[i] <- jsoes_null_sdmTMB_cv$sum_loglik
  
}

ggplot(jsoes_null_cv_comp, aes(x = model_run, y = sum_loglik, color = as.factor(n_folds))) +
  geom_point()

# sometimes this doesn't really work at all?
# drop those points and try again
ggplot(subset(jsoes_null_cv_comp, sum_loglik > -10000), aes(x = model_run, y = sum_loglik, color = as.factor(n_folds))) +
  geom_point()
# so yes, as number of folds goes up, it generally becomes more tightly clustered

set.seed(124)
jsoes_null_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year),
                                  data = csyis,
                                  mesh = mesh_sdmTMB,
                                  time = "year",
                                  spatial = "on",
                                  spatiotemporal = "off",
                                  family = tweedie(link = "log"),
                                  k_folds = 4)

jsoes_null_sdmTMB_cv$sum_loglik

# Model 2: Spatiotemporal, IID
jsoes_st_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year),
                                data = csyis,
                                mesh = mesh_sdmTMB,
                                time = "year",
                                spatial = "on",
                                spatiotemporal = "iid",
                                family = tweedie(link = "log"),
                                k_folds = 4)


# Model 3: Spatiotemporal, AR1
jsoes_stac_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year),
                                  data = csyis,
                                  mesh = mesh_sdmTMB,
                                  time = "year",
                                  spatial = "on",
                                  spatiotemporal = "ar1",
                                  family = tweedie(link = "log"),
                                  k_folds = 4)

# Model 4: Intercept as temporal (AR1), spatiotemporal as IID
jsoes_tempar1_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0,
                               data = csyis,
                               mesh = mesh_sdmTMB,
                               time = "year",
                               time_varying = ~ 1,
                               time_varying_type = "ar1",
                               spatial = "on",
                               spatiotemporal = "iid",
                               family = tweedie(link = "log"),
                               k_folds = 4)


# Model 5: Intercept as temporal (AR1), spatiotemporal as AR1
jsoes_star1_tempar1_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0,
                                     data = csyis,
                                     mesh = mesh_sdmTMB,
                                     time = "year",
                                     time_varying = ~ 1,
                                     time_varying_type = "ar1",
                                     spatial = "on",
                                     spatiotemporal = "ar1",
                                     family = tweedie(link = "log"),
                                     k_folds = 4)


# Compare models
AIC(jsoes_null_sdmTMB)
AIC(jsoes_st_sdmTMB)
AIC(jsoes_stac_sdmTMB)
AIC(jsoes_tempar1_sdmTMB)
AIC(jsoes_star1_tempar1_sdmTMB)




jsoes_null_sdmTMB_cv$sum_loglik
jsoes_st_sdmTMB_cv$sum_loglik
jsoes_stac_sdmTMB_cv$sum_loglik
jsoes_tempar1_sdmTMB_cv$sum_loglik
jsoes_star1_tempar1_sdmTMB_cv$sum_loglik

#### Repeated k-fold cross-validation ####

# As we see in the example above, there is quite a bit of random variation in
# the cross-validation performance of the model depending on how the dataset
# is split randomly.

# for the number of folds, there are not hard and fast rules. We will use 10
dim(csyis)
# we have 1287 rows here. That is a reasonably large dataset and using a default
# value like 10 folds seems reasonable.

# We will run the cross-validation 100 times for each model, with 10 folds for each model
cv_repeats <- 10
n_folds <- 10
n_models <- 5

# create a dataframe to store our results
jsoes_csyis_cv_comp <- data.frame(model_run = rep(NA, cv_repeats*n_models),
                                  model_name = rep(NA, cv_repeats*n_models),
                                 sum_loglik = rep(NA, cv_repeats*n_models))


for (i in 1:cv_repeats){
  # Model 1: Null model
  jsoes_null_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year),
                                    data = csyis,
                                    mesh = mesh_sdmTMB,
                                    time = "year",
                                    spatial = "on",
                                    spatiotemporal = "off",
                                    family = tweedie(link = "log"),
                                    k_folds = n_folds)
  
  jsoes_csyis_cv_comp$model_run[(i-1)*5+1] <- (i-1)*5+1
  jsoes_csyis_cv_comp$model_name[(i-1)*5+1] <- "1-null"
  jsoes_csyis_cv_comp$sum_loglik[(i-1)*5+1] <- jsoes_null_sdmTMB_cv$sum_loglik
  
  # Model 2: Spatiotemporal, IID
  jsoes_st_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year),
                                  data = csyis,
                                  mesh = mesh_sdmTMB,
                                  time = "year",
                                  spatial = "on",
                                  spatiotemporal = "iid",
                                  family = tweedie(link = "log"),
                                  k_folds = n_folds)
  
  jsoes_csyis_cv_comp$model_run[(i-1)*5+2] <- (i-1)*5+2
  jsoes_csyis_cv_comp$model_name[(i-1)*5+2] <- "2-st_iid"
  jsoes_csyis_cv_comp$sum_loglik[(i-1)*5+2] <- jsoes_st_sdmTMB_cv$sum_loglik
  
  
  # Model 3: Spatiotemporal, AR1
  jsoes_stac_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year),
                                    data = csyis,
                                    mesh = mesh_sdmTMB,
                                    time = "year",
                                    spatial = "on",
                                    spatiotemporal = "ar1",
                                    family = tweedie(link = "log"),
                                    k_folds = n_folds)
  
  jsoes_csyis_cv_comp$model_run[(i-1)*5+3] <- (i-1)*5+3
  jsoes_csyis_cv_comp$model_name[(i-1)*5+3] <- "3-st_ar1"
  jsoes_csyis_cv_comp$sum_loglik[(i-1)*5+3] <- jsoes_stac_sdmTMB_cv$sum_loglik
  
  # Model 4: Intercept as temporal (AR1), spatiotemporal as IID
  jsoes_tempar1_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0,
                                       data = csyis,
                                       mesh = mesh_sdmTMB,
                                       time = "year",
                                       time_varying = ~ 1,
                                       time_varying_type = "ar1",
                                       spatial = "on",
                                       spatiotemporal = "iid",
                                       family = tweedie(link = "log"),
                                       k_folds = n_folds)
  
  jsoes_csyis_cv_comp$model_run[(i-1)*5+4] <- (i-1)*5+4
  jsoes_csyis_cv_comp$model_name[(i-1)*5+4] <- "4_temp_ar1_st_iid"
  jsoes_csyis_cv_comp$sum_loglik[(i-1)*5+4] <- jsoes_tempar1_sdmTMB_cv$sum_loglik
  
  
  # Model 5: Intercept as temporal (AR1), spatiotemporal as AR1
  jsoes_star1_tempar1_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0,
                                             data = csyis,
                                             mesh = mesh_sdmTMB,
                                             time = "year",
                                             time_varying = ~ 1,
                                             time_varying_type = "ar1",
                                             spatial = "on",
                                             spatiotemporal = "ar1",
                                             family = tweedie(link = "log"),
                                             k_folds = n_folds)
  
  jsoes_csyis_cv_comp$model_run[(i-1)*5+5] <- (i-1)*5+5
  jsoes_csyis_cv_comp$model_name[(i-1)*5+5] <- "5_temp_ar1_st_ar1"
  jsoes_csyis_cv_comp$sum_loglik[(i-1)*5+5] <- jsoes_star1_tempar1_sdmTMB_cv$sum_loglik
  
  
  
}

ggplot(subset(jsoes_csyis_cv_comp, sum_loglik > -10000), aes(x = model_name, y = sum_loglik)) +
  geom_boxplot()

# Ok, let's explore the distribution of observations, because this might be why we
# have very different performance (and causes CV runs to be absurd)

# distribution has some very large outliers - one tow that was full of Chinook
ggplot(csyis, aes(x = n_per_km)) +
  geom_histogram()



#### Explore covariate data ####

ggplot(jsoes, aes(x = as.factor(nmi_from_shore), y = chinook_salmon_yearling_interior_sp)) +
  geom_boxplot()

ggplot(jsoes, aes(x = stn_depth_m, y = chinook_salmon_yearling_interior_sp)) +
  geom_point()

ggplot(jsoes, aes(x = x3m_temp_c, y = chinook_salmon_yearling_interior_sp)) +
  geom_point()

# look at correlation between nmi_from_shore and stn_depth_m
ggplot(jsoes, aes(x = nmi_from_shore, y = stn_depth_m)) +
  geom_point()
cor(jsoes$nmi_from_shore, jsoes$stn_depth_m)
# yeah pretty correlated, should probably only choose one.
# I am going to say that for now, I think that we'll include distance from shore,
# since for a pelagic species that seems to make more ecological sense

#We also need to know the date range to use for SST for interpolation.

jsoes %>% 
  mutate(julian_day = yday(sample_date)) -> jsoes

hist(jsoes$julian_day)
# for simplicity (and to keep our dataset a reasonable size), let's take the mean day
round(mean(jsoes$julian_day))
# corresponds to june 24 (on non-leap years)

#### Calculate distance to shore ####

dplyr::select(csyis, station, km_from_shore, X, Y) -> csyis_spatial
st_as_sf(csyis_spatial, coords = c("X", "Y"), crs = UTM_zone_10_crs) -> csyis_sf
st_distance(csyis_sf, US_west_coast_proj_km) -> csyis_dist_shore
csyis$sf_dist_shore <- as.numeric(csyis_dist_shore)

# confirm that these are reasonable, based on points from survey
# they're similar - but not identical
hist(as.numeric(csyis$sf_dist_shore)-csyis_spatial$km_from_shore)
# and they're pretty biased - the distances calculated using sf tend to be smaller
# Let's just use the sf distance from shore calculated for both prediction grid and
# the samples themselves


# rescale data
csyis %>% 
  mutate(SST_scaled = as.numeric(scale(x3m_temp_c)),
         dist_shore_scaled = as.numeric(scale(sf_dist_shore))) -> csyis


#### Fit model with covariates ####

# for distance - just pull from prediction grid, for a single year
distance_projection_vector <- as.numeric(filter(survey_predict_grid, year == 1998)$dist_shore_scaled)

# for SST - take prediction grid, add SST with locations as rows and years as column
survey_predict_grid %>% 
  dplyr::select(-dist_shore) %>% 
  pivot_wider(names_from = year, values_from = SST_scaled) %>% 
  dplyr::select(-c(X, Y)) %>% 
  as.matrix() -> temperature_projection_matrix


# sanity check that all our covariates and projection matrices are the same dimensions
length(distance_projection_vector)
dim(temperature_projection_matrix)
dim(A_gs)
# very nice!

# look at projection histograms
hist(distance_projection_vector)
hist(temperature_projection_matrix)

# look at data histograms
hist(csyis$n_per_km)
hist(csyis$SST_scaled)
hist(csyis$dist_shore_scaled)

summary(csyis$SST_scaled)
summary(csyis$dist_shore_scaled)

compile(here::here("analysis", "jsoes_st_cov_sdm.cpp"))
dyn.load( dynlib(here::here("analysis", "jsoes_st_cov_sdm")))

Data = list( "D_i"= csyis$n_per_km, "t_i" = csyis$year - min(csyis$year), 
             "temp_i" = csyis$SST_scaled, "dist_i" = csyis$dist_shore_scaled,
             "dist_g" = distance_projection_vector, "temp_gt" = temperature_projection_matrix, 
             "n_t" = length(unique(csyis$year)),
             "A_is"=A_is, "A_gs"=A_gs, "M0"=spde$c0, "M1"=spde$g1, "M2"=spde$g2 )
Params = list( "beta_t"=rep(0, length(unique(csyis$year))), 
               "beta_temp" = 0, "beta_dist" = 0,
               "ln_tau_omega"=0, "ln_tau_epsilon"=0,
               "ln_kappa"=0, "ln_phi" = 0, 
               "finv_power" = 0, "omega_s"=rnorm(nrow(spde$c0)),
               "epsilon_st"=matrix(0, nrow=nrow(spde$c0), ncol=Data$n_t))

st_cov_Obj = MakeADFun( data=Data, parameters=Params, random= c("omega_s", "epsilon_st"),
                        DLL = "jsoes_st_cov_sdm")


# Optimize
st_cov_Opt = nlminb( start=st_cov_Obj$par, obj=st_cov_Obj$fn, grad=st_cov_Obj$gr )
# Opt$SD = sdreport( Obj, bias.correct=TRUE )
st_cov_report = st_cov_Obj$report()

#### Fit sdmTMB with effect of year and effect of space and spatio-temporal effect
csyis_stac_cov_sdmTMB <- sdmTMB(n_per_km ~ 0 + as.factor(year) + SST_scaled + dist_shore_scaled,
                                data = csyis,
                                mesh = mesh_sdmTMB,
                                time = "year",
                                spatial = "on",
                                spatiotemporal = "ar1",
                                family = tweedie(link = "log"))

# Katie' model
load("/Users/markusmin/Documents/NCCMovementModel/final/sdmTMB/mayjune_bspde.RData")
fit <- sdmTMB(formula = density ~ 0 +
                as.factor(year) + #month +
                scale_depth + I(scale_depth^2) +
                scale_yday + scale_yday2,
              spatial="on",
              spatiotemporal = "iid",
              time = "month_year",
              mesh = bspde,
              family = tweedie(link = "log"),
              data = dat)

sanity(csyis_stac_cov_sdmTMB)
AIC(csyis_stac_cov_sdmTMB)


# check for alignment
csyis_stac_cov_sdmTMB
st_cov_report$beta_t
st_cov_report$beta_dist
st_cov_report$beta_temp
# the same! Now fit for other species

# Cross validation for comparison with other models
set.seed(556)
# csyis_stac_cov_sdmTMB_cv <- sdmTMB_cv(n_per_km ~ 0 + as.factor(year) + SST_scaled + dist_shore_scaled,
#                                       data = csyis,
#                                       mesh = mesh_sdmTMB,
#                                       time = "year",
#                                       spatial = "on",
#                                       spatiotemporal = "ar1",
#                                       family = tweedie(link = "log"),
#                                       k_folds = 4)
# 
# csyis_stac_cov_sdmTMB_cv$sum_loglik


#### Predict using the model with covariates ####

survey_predict_grid
# check for no missing values in prediction grid
which(is.na(survey_predict_grid))


csyis_stac_cov_sdmTMB_pred <- predict(csyis_stac_cov_sdmTMB, newdata = survey_predict_grid)

# plot the predictions
plot_density_predictions <- function(predictions){
  # I'm sure I could dial the geom_tile argument in a bit better based on the actual grid size of 0.05 lat x 0.05 long, but oh well... according to NOAA these tiles should be approx 4 km wide by 6 km tall
  density_predictions_map <- survey_area_basemap_km +
    geom_tile(data = predictions, aes(x = X, y = Y, fill = exp(est)),
              width = 7, height = 7) +
    scale_fill_viridis_c( trans = "sqrt",
                          # trim extreme high values to make spatial variation more visible
                          na.value = "yellow", limits = c(0, quantile(exp(predictions$est), 0.995)),
                          name = "Predicted density\n(n per km)") +
    facet_wrap(~year, nrow = 3) +
    theme(axis.text = element_blank(),
          legend.position = c(0.95, 0.1))
  
  return(density_predictions_map)
  
}

csyis_prediction_map <- plot_density_predictions(predictions = csyis_stac_cov_sdmTMB_pred)

# estimate uncertainty
sim <- predict(csyis_stac_cov_sdmTMB, newdata = survey_predict_grid, nsim = 100)

# generate CV estimates for each year
# start with the first year
sim_year <- sim[survey_predict_grid$year == unique(survey_predict_grid$year)[1], ]
pred_cv_df <- csyis_stac_cov_sdmTMB_pred[csyis_stac_cov_sdmTMB_pred$year == unique(survey_predict_grid$year)[1], ]
pred_cv_df$lwr <- apply(exp(sim_year), 1, quantile, probs = 0.025)
pred_cv_df$upr <- apply(exp(sim_year), 1, quantile, probs = 0.975)
pred_cv_df$sd <- round(apply(exp(sim_year), 1, function(x) sd(x)), 2)
pred_cv_df$cv <- round(apply(exp(sim_year), 1, function(x) sd(x) / mean(x)), 2)

# then, loop through the remaining years and bind them
for (i in 2:length(unique(survey_predict_grid$year))){
  sim_year <- sim[survey_predict_grid$year == unique(survey_predict_grid$year)[i], ]
  pred_year <- csyis_stac_cov_sdmTMB_pred[csyis_stac_cov_sdmTMB_pred$year == unique(survey_predict_grid$year)[i], ]
  pred_year$lwr <- apply(exp(sim_year), 1, quantile, probs = 0.025)
  pred_year$upr <- apply(exp(sim_year), 1, quantile, probs = 0.975)
  pred_year$sd <- round(apply(exp(sim_year), 1, function(x) sd(x)), 2)
  pred_year$cv <- round(apply(exp(sim_year), 1, function(x) sd(x) / mean(x)), 2)
  
  pred_cv_df %>%  
    bind_rows(., pred_year) -> pred_cv_df
}

# Plot all of the years with a facet wrap
ggplot(pred_cv_df, aes(X, Y, fill = cv)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_viridis_c() +
  facet_wrap(~year, nrow = 3)

sim_last <- sim[survey_predict_grid$year == max(survey_predict_grid$year), ] # just plot last year
pred_last <- csyis_stac_cov_sdmTMB_pred[csyis_stac_cov_sdmTMB_pred$year == max(survey_predict_grid$year), ]
pred_last$lwr <- apply(exp(sim_last), 1, quantile, probs = 0.025)
pred_last$upr <- apply(exp(sim_last), 1, quantile, probs = 0.975)
pred_last$sd <- round(apply(exp(sim_last), 1, function(x) sd(x)), 2)
pred_last$cv <- round(apply(exp(sim_last), 1, function(x) sd(x) / mean(x)), 2)

csyis_prediction_cv_map <- ggplot(pred_last, aes(X, Y, fill = cv)) +
  geom_tile(width = 7, height = 7) +
  scale_fill_viridis_c()
# plot the standard error of the density predictions
ggsave(here::here("figures", "csyis_prediction_cv_map.png"), plot = csyis_prediction_cv_map, height = 12, width = 17)




ggsave(here::here("figures", "csyis_prediction_map.png"), plot = csyis_prediction_map, height = 12, width = 17)



#### Replicate Katie's indices ####

# load some objects from Katie's project:
# read in our grid 
new_df <- read.csv("/Users/markusmin/Documents/NCCMovementModel/final/sdmTMB/new_grid.csv")

# subset so our abundance is only over a comparable area to IPM
curr_newdf <- subset(new_df, Lat.km >= 5113 & Lat.km <= 5354)

# now, take our prediction grid and ensure that it covers the same area as Katie's grid
range(survey_predict_grid$Y)
my_newdf <- subset(survey_predict_grid, Y >= 5113 & Y <= 5354)

# here, use our model instead of Katie's and change the covariates that we're using to predict
pred.index = predict(csyis_stac_cov_sdmTMB, newdata = my_newdf, return_tmb_object = TRUE)

# area_vec <- curr_newdf$area_km2
hist(curr_newdf$area_km2)
# so these are the areas of the cells that she's predicting to; they're mostly 2x2 km, but some
# are cropped (I assume because they intersect with the coast)


# we don't have this column, but in our calculations we assumed our cells were each 24 km^2;
# we can pass a scalar instead of a vector if we assume they're all the same size
my_index_total = get_index(pred.index, bias_correct = TRUE, area = 24)

# let's check the area calculations
# the index total area is the number of cells times the area
nrow(subset(pred.index$data, year == 1998))*24

# the area of our survey is approximately a rectangle containing the boundaries

(max(pred.index$data$X) - min(pred.index$data$X)) * (max(pred.index$data$Y) - min(pred.index$data$Y))
# you'd expect this to be bigger. Passes the sanity check



# read in Katie's index_total
katie_index_total <- read.csv("/Users/markusmin/Documents/NCCMovementModel/final/sdmTMB/index_total.csv")

# keep only june for comparison
katie_index_total_june <- subset(katie_index_total, month == 6)

abundance_index_comp <- ggplot(my_index_total, aes(x = year, y = est)) +
  geom_point(color = "blue") + 
  geom_errorbar(aes(ymin = lwr, ymax = upr), color = "blue") +
  geom_point(data = katie_index_total_june, color = "red") +
  geom_errorbar(data = katie_index_total_june, color = "red", aes(ymin = lwr, ymax = upr)) +
  ggtitle("Abundance index - Katie (red), Markus (blue)")

ggsave(here::here("figures", "migration_model", "abundance_index_comparison.png"), plot = abundance_index_comp)
write.csv(my_index_total, here::here("figures", "markus_interior_spring_chinook_yearlings_index.csv"))

### Compare center of gravity 

# our calculations
my_newdf <- subset(survey_predict_grid, Y >= 5113 & Y <= 5354)

pred.index = predict(csyis_stac_cov_sdmTMB, newdata = my_newdf, return_tmb_object = TRUE)
my_cog <- get_cog(pred.index, format = "wide")



# Katie's code
# curr_new_df <- subset(new_df, Lat.km <= 5354 & Lat.km >= 5113)
# 
# pred.index = predict(fit, newdata = curr_new_df, return_tmb_object = TRUE)
# cog <- get_cog(pred.index, format = "wide")
# 
# 
# 
# write.csv(cog, "sdmTMB/cog.csv", row.names=FALSE)




### load Katie's calculations
katie_cog <- read.csv("/Users/markusmin/Documents/NCCMovementModel/final/sdmTMB/cog.csv")

df <- read.csv(here::here("/Users/markusmin/Documents/NCCMovementModel/final/data/df_JSOES_Chinook_salmon_yearling_Interior_Sp.csv")) %>%
  mutate(scale_lat = scale(as.numeric(Lat.km))[ , 1]) %>%
  mutate(scale_yday = scale(as.numeric(yday))[,1]) %>%
  mutate(scale_depth2 = scale_depth^2) %>%
  mutate(scale_yday2 = scale_yday^2)

# figuring out a time variable that's month_year
dat <- df %>% 
  dplyr::rename(density = 'numPkm2') %>%
  mutate(month_year_full = paste(month, year))
dat <- transform(dat, month_year=match(month_year_full, unique(month_year_full))) %>%
  mutate(month = match(month, month.name)) 

# keeping track of month/year combos 
month_year_combos <- data.frame(combos = c(unique(dat$month_year_full))) %>%
  mutate(month = str_trim(substring(combos, 1, 4))) %>%
  mutate(month = match(month, month.name)) %>%
  mutate(year = str_trim(substring(combos,5)))

katie_cog <- katie_cog %>%
  mutate(month_year_full = month_year_combos$combos[month_year]) %>%
  mutate(month = month_year_combos$month[month_year]) %>% 
  mutate(year = as.numeric(gsub("June |May ", "", month_year_full))) # my addition for plotting comparisons


katie_cog_df <- data.frame(month_year = katie_cog$month_year,
                     month = katie_cog$month,
                     year = month_year_combos$year[katie_cog$month_year],
                     sdm_katie_cog = katie_cog$est_y,
                     sdm_se = katie_cog$se_y,
                     ipm_cog = rep(0, length(katie_cog$month_year)),
                     ipm_se = rep(0, length(katie_cog$month_year)))

### compare our cog

cog_comp <- ggplot(my_cog, aes(x = year, y = est_y)) +
  geom_point(color = "blue") +
  geom_errorbar(aes(ymin = lwr_y, ymax = upr_y), color = "blue") +
  geom_point(data = subset(katie_cog, month == 6), color = "red") +
  geom_errorbar(data = subset(katie_cog, month == 6), aes(ymin = lwr_y, ymax = upr_y), color = "red") +
  ggtitle("Center of gravity - Katie (red), Markus (blue)") +
  ylab("Y-coordinate (km north of equator)")

ggsave(here::here("figures", "migration_model", "center_of_gravity_comparison.png"), plot = cog_comp)


### compare ratio north

# Take the whole survey and predict across it
range(survey_predict_grid$Y)

# here, use our model instead of Katie's and change the covariates that we're using to predict
pred.full = predict(csyis_stac_cov_sdmTMB, newdata = survey_predict_grid)

# Get the abundance index by year - this is not bias corrected but is fairly close
# to the output from sdmTMB::get_index()

# for full survey area
annual_abundance_index <- data.frame(year = unique(pred.full$year),
                                     est = rep(NA, length(unique(pred.full$year))))

for (i in 1:length(unique(pred.full$year))){
  pred_annual <- subset(pred.full, year == unique(pred.full$year)[i])
  annual_abundance_index$est[i] <- sum(exp(pred_annual$est) * 24)
}


# for northern point
midpoint <- 5180

pred.full.north <- subset(pred.full, Y >= midpoint)

annual_abundance_index_north <- data.frame(year = unique(pred.full.north$year),
                                     est = rep(NA, length(unique(pred.full.north$year))))

for (i in 1:length(unique(pred.full.north$year))){
  pred_annual_north <- subset(pred.full.north, year == unique(pred.full.north$year)[i])
  annual_abundance_index_north$est[i] <- sum(exp(pred_annual_north$est) * 24)
}

percent_north_df <- data.frame(year = unique(pred.full.north$year),
                               full = annual_abundance_index$est,
                               north = annual_abundance_index_north$est,
                               percent = annual_abundance_index_north$est/annual_abundance_index$est)


density_predictions_map <- survey_area_basemap_km +
  geom_tile(data = pred.full, aes(x = X, y = Y, fill = exp(est)),
            width = 7, height = 7) +
  geom_hline(yintercept = midpoint, lty = 2, color = "yellow") +
  # geom_hline(yintercept = 5354, lty = 2, color = "orange") +
  scale_fill_viridis_c( trans = "sqrt",
                        # trim extreme high values to make spatial variation more visible
                        na.value = "yellow", limits = c(0, quantile(exp(pred.full$est), 0.995)),
                        name = "Predicted density\n(n per km)") +
  facet_wrap(~year, nrow = 3) +
  theme(axis.text = element_blank(),
        legend.position = c(0.95, 0.1))

ggsave(here::here("figures", "migration_model", "percent_north_map.png"), plot = density_predictions_map, height = 12, width = 16)

precent_north_plot <- ggplot(percent_north_df, aes(x = year, y = percent)) +
  geom_line() +
  ggtitle("Percent North of midpoint line")

ggsave(here::here("figures", "migration_model", "percent_north_plot.png"), plot = precent_north_plot, height = 12, width = 16)
