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


map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure.ace")
map_merged <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-merged2.ace")


# reduce jn1 reactivity by avg >LOD titer in BNT/BNT
jn1_table <- logtiterTable(map)["JN.1", srGroups(map) == "BNT/BNT"]
jn1_table <- jn1_table[jn1_table > log2(1.6)]
jn1_react <- mean(jn1_table)


ag_reacts <- rep(0, length(agNames(map)))
ag_reacts[agNames(map) == "JN.1"] <- -jn1_react
ag_reacts[agNames(map) == "P.1.1"] <- -1

map_adjust <- optimizeAgReactivity(map, fixed_ag_reactivities = ag_reacts)
map_adjust <- optimize_and_realign_map(map_adjust, map, optimization_numbers = 1000, option_list = list(ignore_disconnected = TRUE, dim_annealing = TRUE))


merged_adjust <- optimizeAgReactivity(map_merged, fixed_ag_reactivities = ag_reacts)
merged_adjust <- optimize_and_realign_map(merged_adjust, map, optimization_numbers = 1000, option_list = list(ignore_disconnected = TRUE, dim_annealing = TRUE))

# now optimize BA.1 reactivity
ag_reacts[grep("BA.1", agNames(map))] <- NA

mapba1_optim <- optimizeAgReactivity(map, fixed_ag_reactivities = ag_reacts)
mapba1_optim <- optimize_and_realign_map(mapba1_optim, map, optimization_numbers = 1000, option_list = list(ignore_disconnected = TRUE, dim_annealing = TRUE))


mergedba1_optim <- optimizeAgReactivity(map_merged, fixed_ag_reactivities = ag_reacts)
mergedba1_optim <- optimize_and_realign_map(mergedba1_optim, map, optimization_numbers = 1000, option_list = list(ignore_disconnected = TRUE, dim_annealing = TRUE))

stop()

#---------- Plot maps
lims_no_zoom <- Racmacs:::mapPlotLims(mergedba1_optim, sera = TRUE)

xlim_no_zoom <- round(lims_no_zoom$xlim)
ylim_no_zoom <- round(lims_no_zoom$ylim)

png("figures/labelled_map/sfig_hamster_maps_merged_options.png", 9,9, units = 'in', res=300, pointsize = 12)
layout(matrix(c(1:9), ncol = 3, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0.1, 0.1, 0))
plot(map, plot_stress = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)

plot(map_adjust, plot_stress = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)

plot(mapba1_optim, plot_stress = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "C", cex = 1.2)

plot(map_merged, plot_stress = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "D", cex = 1.2)

plot(merged_adjust, plot_stress = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "E", cex = 1.2)

plot(mergedba1_optim, plot_stress = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "F", cex = 1.2)

# Procrustes from hamster to human
plot(procrustesMap(map_merged, map, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "G", cex = 1.2)

plot(procrustesMap(merged_adjust, map_adjust, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "H", cex = 1.2)

plot(procrustesMap(mergedba1_optim, mapba1_optim, sera = FALSE), xlim = xlim_no_zoom, ylim = ylim_no_zoom)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "I", cex = 1.2)

dev.off()
layout(matrix(c(1), ncol = 1, byrow = T))


save.acmap(mergedba1_optim, "data/maps/map-merged-adj.ace")
