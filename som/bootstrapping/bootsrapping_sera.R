library(Racmacs)
set.seed(100)

map_n <- "map-merged-adj.ace"
neut <- read.acmap(paste0("./data/maps/", map_n))


lims_no_zoom <- Racmacs:::mapPlotLims(neut, sera = TRUE)

xlim_no_zoom <- round(lims_no_zoom$xlim)
ylim_no_zoom <- round(lims_no_zoom$ylim)

srGroups(neut)[as.character(srGroups(neut)) %in% c("alpha conv.", "alpha+E484K conv.")] <- "alpha/alpha+E484K conv."


target_srgs <- c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.',
                 'BA.5 conv.', 'CK.2.1.1 conv.', 'D614G conv. (Hamster)', 'delta conv. (Hamster)', 'BA.1 conv. (Hamster)', 'BA.5 conv. (Hamster)', 'XBB.1.5 conv. (Hamster)')

labels <- data.frame(
  row.names = target_srgs,
  val = LETTERS[1:length(target_srgs)]
)

png("som/bootstrapping/bootstrapping-sera.png", width = 8, height = 8, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:18), ncol = 3, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 0.1, 0))

for(srGroup in target_srgs){
    
    print(srGroup)
  
  newMap <- removeSera(neut, srNames(neut)[as.character(srGroups(neut)) == srGroup])
  newMap <- optimizeMap(
    map                     = newMap,
    number_of_dimensions    = 2,
    number_of_optimizations = 1000,
    minimum_column_basis    = "none",
    options =  list(ignore_disconnected = TRUE,
                    dim_annealing = TRUE)
  )
  
  newMap <- realignMap(newMap, neut)
  
  save_text <- strsplit(srGroup, "/")[[1]][1]
  if(srGroup == "AZ/BNT"){
    save_text <- "AZ_BNT"
  }
  save.acmap(map = newMap, filename = paste0("./som/bootstrapping/wo_",save_text,".ace"))
  #srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  title_text <- srGroup
  if(title_text == "mRNA1273/mRNA1273") {
    title_text <- "mRNA-1273/mRNA-1273"
  } else if (title_text == "WT conv.") {
    title_text <- "Ancestral virus conv."
  } else if (title_text == "BA.1 conv.") {
    title_text <- "BA.1 omicron conv."
  } else if (title_text == "BA.2 conv.") {
    title_text <- "BA.2 omicron conv."
  } else if (title_text == "alpha/alpha+E484K conv.") {
    title_text <- "alpha conv."
  } else if(title_text == "BA.5 conv"){
    title_text <- "BA.5 omicron conv."
  }
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, labels[srGroup, ], cex = 1.2)
}

dev.off()


# when maps already exist
png("som/bootstrapping/bootstrapping-sera2.png", width = 12, height = 8, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:20), ncol = 5, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 0.1, 0))

for(srGroup in target_srgs){
  print(srGroup)
  
  save_text <- strsplit(srGroup, "/")[[1]][1]
  if(srGroup == "AZ/BNT"){
    save_text <- "AZ_BNT"
  }
  newMap <- read.acmap(paste0("./som/bootstrapping/wo_",save_text,".ace"))
  srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  title_text <- srGroup
  if(title_text == "mRNA1273/mRNA1273") {
    title_text <- "mRNA-1273/mRNA-1273"
  } else if (title_text == "WT conv.") {
    title_text <- "Ancestral virus conv."
  } else if (title_text == "BA.1 conv.") {
    title_text <- "BA.1 omicron conv."
  } else if (title_text == "BA.2 conv.") {
    title_text <- "BA.2 omicron conv."
  } else if (title_text == "alpha/alpha+E484K conv.") {
    title_text <- "alpha conv."
  } else if(title_text == "BA.5 conv."){
    title_text <- "BA.5 omicron conv."
  }
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+1, ylim_no_zoom[2]-1, labels[srGroup, ], cex = 1.2)
  text(xlim_no_zoom[1] + 2, ylim_no_zoom[2]-1, srGroup, cex = 0.8, adj = 0)
}

dev.off()

