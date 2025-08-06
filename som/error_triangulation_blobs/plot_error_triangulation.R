rm(list = ls())
library(Racmacs)
library(cowplot)
source("functions/map_functions.R")
set.seed(100)

neut <- read.acmap('data/maps/map-merged-adj.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

# subsetMap to only positioned sera
map_positioned <- remove_na_coords(neut)
 

png("revision/error_triangulation_blobs/error_triangulation.png", width = 12, height = 6, units = 'in', res=300, pointsize = 18)
par(mfrow = c(1,2))
plot(map_positioned, show_error = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE, margins = NULL)
  # title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "D", cex = 1.2)

plot(triangulationBlobs(relaxMap(map_positioned), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE, margins = NULL)
  # title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "E", cex = 1.2)

dev.off()


