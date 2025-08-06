library(Racmacs)

map_n <- "map-merged-adj.ace"
neut <- read.acmap(paste0("./data/maps/", map_n))


lims_no_zoom <- Racmacs:::mapPlotLims(neut, sera = TRUE)

xlim_no_zoom <- round(lims_no_zoom$xlim)
ylim_no_zoom <- round(lims_no_zoom$ylim)

do_bootstraps <- function(type,neut, options_list = list(ignore_disconnected = TRUE,
                                                         dim_annealing = TRUE)){
  
  
  neutBoot <- neut
  if(type == "resample"){
    print(paste(type, "starting"))
    neutBoot <- bootstrapMap(
      neut,
      "resample",
      bootstrap_repeats = 500,
      bootstrap_ags = TRUE,
      bootstrap_sr = TRUE,
      reoptimize = TRUE,
      optimizations_per_repeat = 1000,
      ag_noise_sd = 0.7,
      titer_noise_sd = 0.7,
      options = options_list
    )
    print(paste(type, "done"))
  } else if(type == "noisy"){
    print(paste(type, "starting"))
    neutBoot <- bootstrapMap(
      neut,
      "noisy",
      bootstrap_repeats = 500,
      bootstrap_ags = TRUE,
      bootstrap_sr = TRUE,
      reoptimize = TRUE,
      optimizations_per_repeat = 1000,
      ag_noise_sd = 0.7,
      titer_noise_sd = 0.7,
      options = options_list
    )
    print(paste(type, "done"))
  } else if(type == "bayesian"){
    print(paste(type, "starting"))
    neutBoot <- bootstrapMap(
      neut,
      "bayesian",
      bootstrap_repeats = 500,
      bootstrap_ags = TRUE,
      bootstrap_sr = TRUE,
      reoptimize = TRUE,
      optimizations_per_repeat = 1000,
      ag_noise_sd = 0.7,
      titer_noise_sd = 0.7,
      options = options_list
    )
    print(paste(type, "done"))
  }
  
  return(neutBoot)
  
}

neutBoots <- list()
for(type in c("noisy", "resample", "bayesian") ){
  neutBoots[[type]] <- do_bootstraps(type, neut, list(ignore_disconnected = TRUE,
                                                      dim_annealing = TRUE))
}

saveRDS(neutBoots, paste0("som/bootstrapping/", gsub(".ace", ".rds", map_n)))


# change colour
neutBoots <- readRDS(paste0("som/bootstrapping/", gsub(".ace", ".rds", map_n)))

blobs <- lapply(c("noisy", "resample", "bayesian"), function(x){
  agFill(neutBoots[[x]]) <- agFill(neut)
  srOutline(neutBoots[[x]]) <- srOutline(neut)
  bootstrapBlobs(neutBoots[[x]], conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
})



# Plot figure
png("som/bootstrapping/resample-bootstrap_merged.png", width = 9, height = 3, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1, 2, 3), ncol=3))
par(oma=c(0, 0, 0, 0), mar=c(0.2, 0.2, 0.2, 0.2))
lapply(1:length(blobs), function(b){
  plot(blobs[[b]], xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, LETTERS[b], cex = 1.4)
})
dev.off()

