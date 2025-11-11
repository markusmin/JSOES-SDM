# 01.5 PIT tag data visualizations

# Description: This script takes the output of 01_PIT_tag_data.Rmd and creates visualizations
# for the supplemental figures for Min et al. 2026.

library(tidyverse)
library(here)

# load the SAR data - the output of 01_PIT_tag_data.Rmd
SAR_data <- read.csv(here::here("model_inputs", "chinook_det_hist.csv"))

# inspect distribution of outmigration timing

SAR_data %>% 
  mutate(BON_juv_det_date = substr(BON_juv_det_time, 1, 10)) %>% 
  mutate(BON_juv_det_date = ymd(BON_juv_det_date)) %>% 
  mutate(outmigration_date = yday(BON_juv_det_date)) -> SAR_data

# plot distribution of outmigration timing
outmigration_timing_plot <- ggplot(SAR_data, aes(x = outmigration_date)) +
  geom_histogram() +
  facet_wrap(~group + run_name)

ggsave(here::here("figures", "outmigration_timing_plot_by_group.png"), outmigration_timing_plot,
       height = 12, width = 12)


# inspect distribution of years with adult returns
subset(SAR_data, adult_det == 1) -> SAR_adult_det

adult_returns_plot <- ggplot(SAR_adult_det, aes(x = run_year)) +
  geom_histogram(breaks = seq(1998, 2025, 1)) +
  facet_wrap(~group + run_name)

ggsave(here::here("figures", "adult_returns_plot.png"), adult_returns_plot,
       height = 12, width = 12)

SAR_adult_det %>% 
  group_by(group, run_name, run_year) %>% 
  summarise(N = sum(adult_det)) -> adult_det_by_year

# inspect distribution of juvenile run years

adult_returns_plot <- ggplot(SAR_data, aes(x = run_year)) +
  geom_histogram(breaks = seq(1998, 2025, 1)) +
  facet_wrap(~group + run_name)

ggsave(here::here("figures", "annual_juveniles_BON_det_plot.png"), adult_returns_plot,
       height = 12, width = 12)


#### Create visualizations for our two target groups: Upper Columbia and Snake River Fall Chinook ####

### show outmigration timing, with June 15 clearly demarcated

SRF_outmigration_plot <- ggplot(subset(SAR_data, run_name == "Fall" & group == "Snake River"), aes(x = outmigration_date)) +
  geom_histogram() +
  geom_vline(xintercept = 167, lty = 2, color = "red") +
  annotate(geom = "text", x = 168, y = 150000, label = "June 15", hjust = 0, color = "red") +
  xlab("Outmigration Date (Julian Day)") +
  ylab("Count") +
  ggtitle("Snake River Fall Chinook")

ggsave(here::here("figures", "PIT_tag_data", "SRF_outmigration_timing_plot.png"), SRF_outmigration_plot,
       height = 8, width = 8)

UCF_outmigration_plot <- ggplot(subset(SAR_data, run_name == "Fall" & group == "Upper Columbia"), aes(x = outmigration_date)) +
  geom_histogram() +
  geom_vline(xintercept = 167, lty = 2, color = "red") +
  annotate(geom = "text", x = 168, y = 5000, label = "June 15", hjust = 0, color = "red") +
  xlab("Outmigration Date (Julian Day)") +
  ylab("Count") +
  ggtitle("Upper Columbia Fall Chinook")

ggsave(here::here("figures", "PIT_tag_data", "UCF_outmigration_timing_plot.png"), UCF_outmigration_plot,
       height = 8, width = 8)

UCS_outmigration_plot <- ggplot(subset(SAR_data, run_name == "Summer" & group == "Upper Columbia"), aes(x = outmigration_date)) +
  geom_histogram() +
  geom_vline(xintercept = 167, lty = 2, color = "red") +
  annotate(geom = "text", x = 168, y = 23000, label = "June 15", hjust = 0, color = "red") +
  xlab("Outmigration Date (Julian Day)") +
  ylab("Count") +
  ggtitle("Upper Columbia Summer Chinook")

ggsave(here::here("figures", "PIT_tag_data", "UCS_outmigration_timing_plot.png"), UCS_outmigration_plot,
       height = 8, width = 8)

### show data availability by year - juvenile and adult counts
detection_colors <- c("Adult" = "black",
                      "Juvenile" = "black")

detection_fills <- c("Adult" = "black",
                      "Juvenile" = "white")

SAR_data %>% 
  group_by(group, run_name, run_year) %>% 
  summarise(N_adult = sum(adult_det),
            N_juv = sum(juv_det)) %>% 
  pivot_longer(cols = c(N_adult, N_juv)) %>% 
  mutate(name = ifelse(name == "N_adult", "Adult", "Juvenile")) %>% 
  dplyr::rename(Detections = name) -> SAR_det_by_year

UCS_detections_plot <- ggplot(subset(SAR_det_by_year, run_name == "Summer" & group == "Upper Columbia"), 
                              aes(x = run_year, y = value, color = Detections, fill = Detections)) +
  scale_fill_manual(values = detection_fills) +
  scale_color_manual(values = detection_colors) +
  geom_bar(stat = "identity", position = "dodge2") +
  scale_x_continuous(expand = c(0,0), limits = c(1997, 2026)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 18000)) +
  ylab("Detections") +
  xlab("Run Year") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

ggsave(here::here("figures", "PIT_tag_data", "UCS_detections_plot.png"), UCS_detections_plot,
       height = 8, width = 12)

SRF_detections_plot <- ggplot(subset(SAR_det_by_year, run_name == "Fall" & group == "Snake River"), 
                              aes(x = run_year, y = value, color = Detections, fill = Detections)) +
  scale_fill_manual(values = detection_fills) +
  scale_color_manual(values = detection_colors) +
  geom_bar(stat = "identity", position = "dodge2") +
  scale_x_continuous(expand = c(0,0), limits = c(1997, 2026)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 85000)) +
  ylab("Detections") +
  xlab("Run Year") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

ggsave(here::here("figures", "PIT_tag_data", "SRF_detections_plot.png"), SRF_detections_plot,
       height = 8, width = 12)








