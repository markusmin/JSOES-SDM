# Update website

rmarkdown::render(here::here("docs", "index.Rmd"))
rmarkdown::render(here::here("docs", "project_description.Rmd"))
rmarkdown::render(here::here("docs", "bongo_exploratory.Rmd"))
rmarkdown::render(here::here("docs", "jsoes_trawl_exploratory.Rmd"))
rmarkdown::render(here::here("docs", "jsoes_survey_description.Rmd"))
rmarkdown::render(here::here("docs", "CCES_survey_description.Rmd"))
rmarkdown::render(here::here("docs", "PRS_survey_description.Rmd"))
rmarkdown::render(here::here("docs", "hake_survey_description.Rmd"))