# plot_distribution - function to plot the density of a taxon

plot_distribution_PRS_PWCC <- function(data, taxon_name){
  # create facet_wrap plot for distribution across all years
  data %>% 
    mutate(encounter = ifelse(total == 0, "zero", "non-zero")) -> data
  
  # Drop 1998 (no water jellies recorded in that year)
  data %>% 
    filter(year != 1998) -> data
  
  survey_area_basemap_km +
    geom_point(data = data, aes(x = X, y = Y, size = total, color = encounter),
               alpha = 0.5) +
    scale_color_manual(values = c("zero" = "#fc9272", "non-zero" = "#2ca25f")) +
    facet_wrap(~year, nrow = 2) +
    ggtitle(taxon_name) +
    theme(legend.position = "right",
          # legend.position = c(0.925, 0.20),
          legend.key.height = unit(0.35, "cm"),
          legend.key.width = unit(0.25, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(0.01, 'cm')) -> species_distribution_plot
  
  return(species_distribution_plot)
}
