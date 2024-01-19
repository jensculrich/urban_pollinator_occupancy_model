# general idea and code apadted from:
# https://aammd.github.io/posts/2023-11-08-tree-cooling/#the-dag

library(brms)
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(cmdstanr)
library(ggdag)

dagified <- dagify(
  occurrence ~ herbaceous_flower_abundance,
  occurrence ~ herbaceous_flower_diversity,
  occurrence ~ soil,
  occurrence ~ vegetation_height,
  occurrence ~ woody_plant_abundance,
  occurrence ~ woody_plant_diversity,
  herbaceous_flower_abundance ~ mowing_enhancement,
  herbaceous_flower_diversity ~ mowing_enhancement,
  soil ~ mowing_enhancement,
  vegetation_height ~ mowing_enhancement
  labels = c(
    "occurrence" = "occurrence",
    "herbaceous_flower_abundance" = "herbaceous\n flower abundance",
    "herbaceous_flower_diversity" = "herbaceous\n flower diversity",
    "soil" = "soil",
    "vegetation_height" = "vegetation\n height",
    "woody_flower_abundance" = "woody\n flower abundance",
    "woody_flower_diversity" = "woody\n flower diversity",
    "mowing_enhancement" = "no-mow\n enhancement"
  ),
  exposure = 'past_land_use',
  outcome = 'cooling',
  coords = list(x = c(occurrence = 0, herbaceous_flower_abundance = -1, herbaceous_flower_diversity = 0, 
                      woody_flower_abundance = 1, woody_flower_diversity = 2,
                      vegetation_height = 1, soil = 0, mowing_enhancement = 0),
                y = c(occurrence = 3, herbaceous_flower_abundance = 2, herbaceous_flower_diversity = 2, 
                      woody_flower_abundance = 2, woody_flower_diversity = 2,
                      vegetation_height = 1, soil = 1, mowing_enhancement = 0))) %>%
  tidy_dagitty() %>%
  mutate(status = case_when(name == "occurrence" ~ 'outcome',
                            name == "mowing_enhancement" ~ 'exposure',
                            .default = 'NA'))

ggplot(dagified, aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_dag() + 
  geom_dag_point(aes(color = status)) +
  geom_dag_label_repel(aes(label = label, fill = status),
                       color = "white", fontface = "bold") +
  geom_dag_edges() + 
  scale_fill_manual(values = c('darkseagreen', 'grey', 'lightskyblue')) + 
  scale_colour_manual(values = c('darkseagreen', 'grey', 'lightskyblue')) + 
  theme(legend.position = 'none')