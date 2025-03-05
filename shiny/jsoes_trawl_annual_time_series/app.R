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
jsoes <- clean_names(read_excel("Markus_Min_Trawl_CTD_Chl_Nuts_Thorson_Scheuerell_5.2.24_FINAL.xlsx"))

species_column_names <- colnames(jsoes)[32:ncol(jsoes)]

# change nmi to km to keep units consistent
jsoes %>% 
  mutate(km_from_shore = nmi_from_shore * 1.852) -> jsoes

jsoes %>% 
  pivot_longer(., cols = all_of(species_column_names), values_to = "n_per_km", names_to = "species") %>% 
  mutate(n = n_per_km*trawl_dist_km) -> jsoes_long

# Drop any samples that don't have spatial data
filter(jsoes_long, !(is.na(mid_long))) -> jsoes_long
# this drops one two samples: 061999CR35 and 062416QR14

# turn this into an sf object
st_as_sf(jsoes_long, coords = c("mid_long", "mid_lat"), crs = 4326) -> jsoes_long_sf

# change CRS to UTM zone 10 (to work in meters)
UTM_zone_10_crs <- 32610
jsoes_long_sf_proj <- sf::st_transform(jsoes_long_sf, crs = UTM_zone_10_crs)

# make this projection into kilometers to help with interpretability
jsoes_long_sf_proj_km <- st_as_sf(jsoes_long_sf_proj$geometry/1000, crs = UTM_zone_10_crs)

# extract geometry
as.data.frame(st_coordinates(jsoes_long_sf_proj_km)) -> jsoes_long_km
# add this back to jsoes_long (X and Y now represent eastings and northings in km)
bind_cols(jsoes_long, jsoes_long_km) -> jsoes_long

# Handle repeat tows: If a tow was repeated, keep only the repeat
jsoes_long %>% 
  mutate(tow_year = paste0(station, "_", year)) -> jsoes_long
repeat_tows <- unique(subset(jsoes_long, `repeat` == TRUE)$tow_year)
jsoes_long %>% 
  filter(!(tow_year %in% repeat_tows & `repeat` == FALSE)) -> jsoes_long


# get all unique tows
jsoes_long %>% 
  dplyr::select(-c(species, n_per_km, n)) %>% 
  filter(!(duplicated(station_code))) -> jsoes_tows



trawl_species <- unique(jsoes_long$species)

#### Define UI ####
ui <- fluidPage(

    # Application title
    titlePanel("JSOES Trawl, annual abundance time series"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          selectInput("species", "Choose a species:",
                      choices = trawl_species)
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
      dataset <- subset(jsoes_long, species == input$species)
      dataset %>% 
        filter(year != 1998)  %>% 
        mutate(log_density = log(n_per_km + 1)) %>% 
        group_by(year) %>% 
        summarise(mean_density = mean(log_density),
                  sd_density = sd(log_density)) %>% 
        # mutate(CI_95_upper = mean_density + 1.96*sd_density) %>% 
        # mutate(CI_95_lower = mean_density - 1.96*sd_density) %>% 
        mutate(taxon = input$choices) -> ts
      
        # generate a plot for one species based on input$choices from ui.R
      ggplot(ts, aes(x = year, y = mean_density)) +
        geom_point() +
        geom_line() +
        ylab("Mean log (n per km + 1)")
        # scale_color_tableau(palette = "Tableau 10")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

# Deploy the app
# library(rsconnect)
# rsconnect::deployApp('shiny/jsoes_trawl_annual_time_series')
