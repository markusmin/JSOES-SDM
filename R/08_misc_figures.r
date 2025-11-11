### 08_misc_figures

## load libraries
library(tidyverse)
library(here)
library(TMB)
library(ggforce)
library(sdmTMB)
library(kableExtra)
library(viridis)

## make meshes

source(here::here("R", "CCES_make_mesh.R"))
source(here::here("R", "hake_survey_make_mesh.R"))
source(here::here("R", "PRS_PWCC_make_mesh.R"))
source(here::here("R", "JSOES_seabirds_make_mesh.R"))
source(here::here("R", "JSOES_make_mesh.R"))

## load data

anchovy_data <- subset(cces_long, species == "Northern Anchovy")





#### Survey coverage comparison figure ####

csyif %>% 
  # filter(!(duplicated(station)))
  filter(year == 2017) %>% 
  mutate(survey = "JSOES") -> jsoes_stations

birds_long %>% 
  filter(species == "common_murre") %>% 
  filter(year == 2017) %>% 
  mutate(survey = "JSOES Seabirds")  -> bird_stations

anchovy_data %>% 
  filter(year == 2017) %>% 
  filter(lat >= 44) %>% 
  mutate(survey = "CCES") %>% 
  mutate(month = month(datetime)) -> cces_stations

rf %>% 
  filter(year == 2011) %>% 
  filter(Y >= 4880)%>% 
  mutate(survey = "PRS") -> prs_stations

hake %>% 
  filter(year == 2013) %>% 
  filter(Lat >= 44) %>% 
  mutate(survey = "Hake") %>% 
  mutate(Date = ymd(Date)) %>% 
  mutate(month = month(Date)) -> hake_stations

dplyr::select(jsoes_stations, survey, X, Y) %>% 
  # bind_rows(dplyr::select(bird_stations, survey, X, Y)) %>% 
  bind_rows(dplyr::select(cces_stations, survey, X, Y)) %>% 
  bind_rows(dplyr::select(prs_stations, survey, X, Y)) %>% 
  bind_rows(dplyr::select(hake_stations, survey, X, Y)) -> survey_data

survey_data$survey <- factor(survey_data$survey, levels = c("JSOES", "PRS", "CCES", "Hake"))

survey_data_no_CCES <- subset(survey_data, survey != "CCES")

# survey_area_basemap_km +
#   geom_point(data = jsoes_stations, aes(x = X, y = Y))
# 
# survey_area_basemap_km +
#   geom_point(data = bird_stations, aes(x = X, y = Y))
# 
# survey_area_basemap_km +
#   geom_point(data = cces_stations, aes(x = X, y = Y))
# 
# survey_area_basemap_km +
#   geom_point(data = prs_stations, aes(x = X, y = Y))
# 
# survey_area_basemap_km +
#   geom_point(data = hake_stations, aes(x = X, y = Y))

survey_comp_plot <- survey_area_basemap_km +
  geom_point(data = survey_data_no_CCES, aes(x = X, y = Y), size = 1) +
  facet_wrap(~survey, nrow = 1) +
  theme(strip.text.x = element_text(size = 20),
        strip.background =element_rect(fill="white"))

ggsave(here::here("figures", "presentation_figures", "survey_comp_plot.png"), survey_comp_plot,  
       height = 6, width = 12)


# Compare just hake and JSOES
survey_area_basemap_km +
  geom_point(data = jsoes_stations, aes(x = X, y = Y), size = 1, color = "red") +
  geom_point(data = hake_stations, aes(x = X, y = Y), size = 1, color = "blue") +
  theme(strip.text.x = element_text(size = 20),
        strip.background =element_rect(fill="white"))


## Plot comparison

# make a fake RF mesh just for visualization purposes
rf_new_bnd <- INLA::inla.nonconvex.hull(cbind(subset(rf_samples_study_domain, sf_dist_shore < 45*1.852 & Y > 4885)$X, subset(rf_samples_study_domain, sf_dist_shore < 45*1.852 & Y > 4885)$Y), convex = -0.1)

inla_mesh_cutoff15_jsoes_domain_fake <- fmesher::fm_mesh_2d_inla(
  loc = cbind(rf_jsoes_all_points$X, rf_jsoes_all_points$Y),
  cutoff = 15,
  boundary = rf_new_bnd
)

# transform mesh to sf
fm_as_sfc(inla_mesh_cutoff15_jsoes_domain_fake) %>% 
  st_set_crs(st_crs(survey_domain_jsoes_cov_rf_years_grid)) -> inla_mesh_cutoff15_jsoes_domain_fake_sf

PRS_v_JSOES_samples <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA, color = "red") +
  geom_sf(data = inla_mesh_cutoff15_jsoes_domain_fake_sf, fill = NA, color = "blue") +
  geom_point(data = jsoes_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3) +
  geom_point(data = subset(rf, Y > 4885), aes(x = X, y = Y), color = "white",fill = "blue", shape = 24, size = 3)


ggsave(here::here("figures", "presentation_figures", "PRS_v_JSOES_samples.png"), PRS_v_JSOES_samples,  height = 8, width = 4)

# new versions of figures where they're separate

JSOES_samples_mesh_map <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_jsoes_domain_fake_sf, fill = NA, color = NA) +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA, color = "red") +
  geom_point(data = jsoes_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)
  # geom_point(data = subset(rf, Y > 4885), aes(x = X, y = Y), color = "white",fill = "blue", shape = 24, size = 3)

ggsave(here::here("figures", "presentation_figures", "JSOES_samples_mesh_map.png"), JSOES_samples_mesh_map,  height = 8, width = 4)

rf_samples_mesh_map <- survey_area_basemap_km +
  geom_sf(data = inla_mesh_cutoff15_sf, fill = NA, color = NA) +
  geom_sf(data = inla_mesh_cutoff15_jsoes_domain_fake_sf, fill = NA, color = "blue") +
  # geom_point(data = jsoes_samples, aes(x = X, y = Y), color = "white",fill = "red", shape = 24, size = 3)
geom_point(data = subset(rf, Y > 4885), aes(x = X, y = Y), color = "white",fill = "blue", shape = 24, size = 3)

ggsave(here::here("figures", "presentation_figures", "rf_samples_mesh_map.png"), rf_samples_mesh_map,  height = 8, width = 4)





