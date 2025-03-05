#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

#### Prepare data for inputs to shiny ####
# Load packages
library(readxl)
library(here)
library(janitor)
library(tidyverse)
library(kableExtra)
library(sf)
library(gstat)
library(sp)
library(ggthemes)
library(ape)

#### Load and reformat trawl data ####
bongo <- clean_names(read_excel("Markus_Min_Plankton Density for Trophic Summary Query_10.11.24.xlsx"))
# Drop any samples that don't have spatial data
filter(bongo, !(is.na(dec_long))) -> bongo
# this drops one sample: 062904NH30 (NH30 in June 2004)

# turn this into an sf object
st_as_sf(bongo, coords = c("dec_long", "dec_lat"), crs = 4326) -> bongo_sf

# change CRS to UTM zone 10 (to work in meters)
UTM_zone_10_crs <- 32610
bongo_sf_proj <- sf::st_transform(bongo_sf, crs = UTM_zone_10_crs)

# make this projection into kilometers to help with interpretability
bongo_sf_proj_km <- st_as_sf(bongo_sf_proj$geometry/1000, crs = UTM_zone_10_crs)

# extract geometry
as.data.frame(st_coordinates(bongo_sf_proj_km)) -> bongo_km
# add this back to bongo (X and Y now represent eastings and northings in km)
bind_cols(bongo, bongo_km) -> bongo

trophic_groupings <- clean_names(read_excel("Markus_Min_Trophic Groupings_10.2.24.xlsx"))

# change genus species to be consistent across these two files, and usual 
# latin name and common name capitalization
bongo$genus_species <- str_to_sentence(bongo$genus_species)
trophic_groupings$genus_species <- str_to_sentence(trophic_groupings$genus_species)
trophic_groupings$common_name <- str_to_title(trophic_groupings$common_name)

trophic_groupings %>% 
  dplyr::select(genus_species, common_name) %>% 
  filter(!(duplicated(genus_species))) -> common_name_df

# get all unique tows
bongo %>% 
  mutate(sampleID = paste0(station, "_", sample_date)) %>% 
  relocate(sampleID) %>% 
  dplyr::select(-c(genus_species, life_history_stage, id_code, count, sum_of_density_number_m3)) %>% 
  filter(!(duplicated(sampleID))) -> bongo_tows

# complete bongo data, so that we have zeros for each taxon in tows where they weren't caught
bongo %>% 
  mutate(sampleID = paste0(station, "_", sample_date)) %>% 
  relocate(sampleID) %>% 
  # drop tow information; we will add this back later
  dplyr::select(-c(id_code, cruise_number, month, year, station_code, net_type, sample_date,
                   sample_time, dec_lat, dec_long, station_depth_m, transect_name, station_number, 
                   station, n_cr_s, day_night, volume_filtered_m3)) %>% 
  complete(sampleID, 
           nesting(genus_species, life_history_stage), 
           fill = list(count = 0, sum_of_density_number_m3 = 0)) -> bongo_catch

bongo_catch %>% 
  left_join(bongo_tows, by = "sampleID") -> bongo_catch





# add common name
bongo_catch %>% left_join(common_name_df, by = "genus_species") -> bongo_catch

# create a version of this where we collapse all life history stages together
bongo_catch %>% 
  dplyr::select(sampleID, genus_species, count, sum_of_density_number_m3) %>% 
  group_by(sampleID, genus_species) %>% 
  summarise_if(is.numeric, sum) -> bongo_catch_LH_collapse

bongo_catch_LH_collapse %>% 
  left_join(bongo_tows, by = "sampleID") %>% 
  ungroup() -> bongo_catch_LH_collapse



bongo_species <- unique(bongo_catch_LH_collapse$genus_species)

#### Define UI ####
ui <- fluidPage(
  
  # Application title
  titlePanel("JSOES Bongo, annual abundance time series"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectInput("species", "Choose a species:",
                  choices = bongo_species)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("ts_plot")
    )
  )
)

#### Define server logic ####
server <- function(input, output) {
  #### filter data and generate plot ####
  
  output$ts_plot <- renderPlot({
    dataset <- subset(bongo_catch_LH_collapse, genus_species == input$species)
    dataset %>% 
        mutate(log_density = log(sum_of_density_number_m3 + 1)) %>% 
        group_by(year) %>% 
        summarise(mean_density = mean(log_density),
                  sd_density = sd(log_density)) %>% 
        # mutate(CI_95_upper = mean_density + 1.96*sd_density) %>% 
        # mutate(CI_95_lower = mean_density - 1.96*sd_density) %>% 
        mutate(taxon = input$species) -> ts
    
    # generate a plot for one species based on input$species from ui.R
      ggplot(ts, aes(x = year, y = mean_density)) +
        geom_point() +
        geom_line() +
        ylab("Mean log (density/m3 + 1)")
    # scale_color_tableau(palette = "Tableau 10")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

# Deploy the app
# library(rsconnect)
# rsconnect::deployApp('shiny/jsoes_bongo_annual_time_series')
