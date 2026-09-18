library(tidyverse)

x <- read_csv("../../data/little_rancheria/caribou.csv") |>
  filter(id %in% sample(unique(id), 5)) |>
  select(id, timestamp, longitude, latitude, season, migration, elevation, phase)
 
write_csv(x, "demo.csv")
