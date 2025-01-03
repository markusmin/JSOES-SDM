# Update website

rmarkdown::render(here::here("docs", "index.Rmd"))
rmarkdown::render(here::here("docs", "project_description.Rmd"))
rmarkdown::render(here::here("docs", "bongo_exploratory.Rmd"))