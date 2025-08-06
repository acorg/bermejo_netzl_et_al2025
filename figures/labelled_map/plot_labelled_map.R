# labelled map 
rm(list = ls())
library(Racmacs)
source("functions/map_functions.R")

# read in map
map <- read.acmap('data/maps/map-merged-adj.ace')
agSize(map)[grepl("E484", agNames(map))] <- 11
# subsetMap to only positioned sera
map_positioned <- remove_na_coords(map)

# reverse drawing order of points
# n_objects <- length(ptDrawingOrder(map_positioned))
# n_sera <- length(srNames(map_positioned))
# ptDrawingOrder(map_positioned)[(n_sera+1):n_objects] <- 1:length(agNames(map_positioned))


map <- map_positioned

lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
lims_no_zoom <- Racmacs:::mapPlotLims(map, sera = TRUE)
#xlim_zoom <- read.csv("./data/metadata/xlim_zoom.csv")$x
#ylim_zoom <- read.csv("./data/metadata/ylim_zoom.csv")$x
xlim_zoom <- round(lims$xlim)
ylim_zoom <- round(lims$ylim)

# ylim_zoom[1] <- ylim_zoom[1]-0.2
#lims <- Racmacs:::mapPlotLims(map_old, sera = TRUE)
#xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
#ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

xlim_no_zoom <- round(lims_no_zoom$xlim)
ylim_no_zoom <- round(lims_no_zoom$ylim)

# Setup plotting function
doplot <- function(map, xlims, ylims, show_labels = TRUE) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  # Plot the regular map
  srOutlineWidth(map) <- 0.8
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.6)
  plot(map, xlim = xlims, 
       ylim =ylims, fill.alpha = 0.9, plot_stress = TRUE)
  
  # Plot labels
  label_adjustments <- matrix(0, numAntigens(map), 2)
  rownames(label_adjustments) <- agNames(map)
  label_adjustments["B.1.351",] <- c(0.9, 0)
  label_adjustments["P.1.1",] <- c(-0.9, 0)
  label_adjustments["B.1.1.7+E484K",] <- c(0.9, -0.5)
  label_adjustments["BA.1",] <- c(-0.4, 0.7)
  label_adjustments["BA.2",] <- c(0, -0.6)
  label_adjustments["B.1.1.7",] <- c(0.3, -0.6)
  label_adjustments["D614G",] <- c(0, -0.5)
  label_adjustments["B.1.617.2",] <- c(0,-0.6)
  label_adjustments["BA.5.3.2",] <- c(0, 0.7)
  
  labels <- agNames(map)
  names(labels) <- agNames(map)
  labels["B.1.351"] <- "beta\n(B.1.351)"
  labels["P.1.1"] <- "gamma\n(P.1.1)"
  labels["BA.1"] <- "BA.1 omicron\n(B.1.1.529+BA.1)"
  labels["BA.2"] <- "BA.2 omicron\n(B.1.1.529+BA.2)"
  labels["B.1.617.2"] <- "delta\n(B.1.617.2)"
  labels["B.1.1.7"] <- "alpha\n(B.1.1.7)"
  labels["B.1.1.7+E484K"] <- "alpha + E484K\n(B.1.1.7+E484K)"
  labels["BA.5.3.2"] <- "BA.5.3.2 omicron\n(B.1.1.529+BA.5)"
  
  label_size <- rep(1, numAntigens(map))
  names(label_size) <- agNames(map)
  if(show_labels){
    text(
     agCoords(map) + label_adjustments,
     cex = label_size,
     label = labels,
     font = 1
    )
  }
  
  
}

png("figures/labelled_map/map_zoom_merged.png", 5, 4, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map, xlim_zoom, ylim_zoom, FALSE)
dev.off()

png("figures/labelled_map/map_no_zoom_merged.png", 7, 7, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map, xlim_no_zoom, ylim_no_zoom, FALSE)
dev.off()

Racmacs::view(map)
