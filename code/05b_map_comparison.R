# load required packages
rm(list = ls())
library(Racmacs)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)

source("./functions/map_functions.R")
# set seed
set.seed(100)

match_names <- read.csv("data/metadata/wang_et_al_name_matching.csv", row.names = "old_name")

map <- read.acmap("data/maps/map-merged-adj.ace")

human_only_map <- read.acmap("data/maps/roessler_netzl_et_al2023.ace")

# compare to human map, NHP map and wang map. need to make wang map

wang_data <- read_xlsx("data/titer_data/wang_et_al_data.xlsx", sheet = "Data for Figure 3", range = "C51:FM80") %>%
  as.data.frame()
rownames(wang_data) <- wang_data[,1]
wang_data <- wang_data[,-1]
wang_data[is.na(wang_data)] <- "*"

wang_map <- make.acmap(wang_data, number_of_dimensions = 2, number_of_optimizations = 2000)
agNames(wang_map) <- match_names[agNames(wang_map), ]
wang_map <- realignMap(wang_map, map)


nhp_map <- read.acmap("data/maps/roessler_netzl_et_al2025_map.ace")

replace_names <- agNames(nhp_map) %in% rownames(match_names)
agNames(nhp_map)[replace_names] <- match_names[agNames(nhp_map)[replace_names],]


# save procrustes
lims_no_zoom <- Racmacs:::mapPlotLims(map, sera = FALSE)

xlim_no_zoom <- round(lims_no_zoom$xlim)
ylim_no_zoom <- round(lims_no_zoom$ylim)
ylim_no_zoom[1] <- ylim_no_zoom[1]-1

png("figures/labelled_map/map_procrustes_to_published.png", 10,5, units = 'in', res=300, pointsize = 12)
par(mfrow = c(1, 3))
plot(procrustesMap(map, human_only_map, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom,
     margins = NULL)

plot(procrustesMap(map, wang_map, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom,
     margins = NULL)

plot(procrustesMap(map, nhp_map, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom,
     margins = NULL)

dev.off()


