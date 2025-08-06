#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(ablandscapes)
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

set.seed(100)

fit_ags <- "sub" # options are "all", "sub", "sub_old"
sub_ags <- c("D614G", "BA.1", "BA.2", "BA.5.3.2", "XBB.1.5", "BA.2.86", "JG.3", "JN.1")

source("./functions/map_longinfo.R")
source("./functions/sams_landscape_functions.R")

figure_dir <- file.path("figures", "landscapes", "gmt_landscapes")
suppressWarnings(dir.create(figure_dir, recursive = T))

ags_to_exclude <- c("")
# Read the base map
map <- read.acmap("./data/maps/map-merged-adj.ace")
map <- removeAntigens(map, ags_to_exclude)
lims <- Racmacs:::mapPlotLims(map, sera = FALSE)

if(fit_ags == "sub"){
  padding <- 0.25
  ags_to_fit_lndscp <- agNames(map)[(agNames(map) %in% sub_ags)]
} else if(fit_ags == "sub_old") {
  padding <- 0.25
  sub_ags <- c("D614G", "BA.1", "BA.2", "BA.5.3.2")
  ags_to_fit_lndscp <- agNames(map)[(agNames(map) %in% sub_ags)]
} else {
  ags_to_fit_lndscp <- agNames(map)
}

# read the full map
sr_colors <- read.csv(file = "./data/metadata/sr_group_colors_new.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",", row.names = "SerumGroup")

# do groupings by infection event
event_groupings <- read.csv("data/titer_data/2400917_Omicron VIII_data-event_groupings.csv") %>%
  mutate(sr_group = paste(event, timepoint))
rownames(event_groupings) <- event_groupings$sampleID

# get the exposure groups and titrated variants
map_long <- read.csv("data/titer_data/biv_boost_long.csv") %>%
  filter(sr_group %in% c("XBB Boost", "no XBB Boost", "biv. Boost")) %>%
  mutate(ag_name = ag,
         sr_name = sr_info_complete,
         sr_group_small = sr_group) %>%
  filter(ag_name %in% sub_ags)

map_long$sr_group <- event_groupings[map_long$sample_id_nr, "sr_group"]

map_long %>%
  mutate(sr_name = paste(sample.ID, sr_group, sr_info, sep = "_")) -> map_long

# adjust titers based on map reactivities
ag_reacts <- agReactivityAdjustments(map)
names(ag_reacts) <- agNames(map)

map_long %>%
  mutate(logtiter = ifelse(titer == "<16", log2(0.8), log2(as.numeric(titer)/10)),
         logtiter_adj = ifelse(titer == "<16", log2(0.8) + ag_reacts[ag_name], logtiter + ag_reacts[ag_name]),
         titer = ifelse(titer == "<16", paste0("<", round(2^logtiter_adj*20)), as.character(2^logtiter_adj*10)),
         logtiter = logtiter_adj) -> map_long


  
sr_group_gmt <-  read.csv("data/titer_data/lod2_sr_group_gmt_multi_adj.csv") %>%
  select(!X) %>%
  mutate(sr_group = gsub("_adj", "", sr_group))


map_long %>%
  select(titer, ag_name, sr_name, sr_group) -> titerdata

titerdata %>%
  group_by(
    sr_group
  ) -> titerdata


titerdata %>%
  mutate(n_ags = length(ags_to_fit_lndscp)) %>%
  group_by(sr_name) %>%
  mutate(n_undetectable = length(titer[titer %in% c("*", "<16")]),
         all_undetectable = n_undetectable == n_ags) %>%
  ungroup() %>%
  filter(all_undetectable) %>%
  select(sr_name) %>%
  unique()-> undetectable_sera

titerdata %>%
  filter(!(sr_name %in% undetectable_sera$sr_name)) ->titerdata

titerdata %>%
  group_map(
    get_titertable
  ) -> titertables


lndscp_fits <- lapply(
  titertables,
  function(titertable) {
    
    ablandscape.fit(
      titers = titertable,
      bandwidth = 1,
      degree = 1,
      method = "cone",
      error.sd = 1,
      acmap = map,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
    
  }
)

titertables_groups <- group_data(titerdata)

# Add impulses
titerdata$gmt <- sr_group_gmt$logtiter[match(interaction(titerdata$sr_group, titerdata$ag_name), interaction(sr_group_gmt$sr_group, sr_group_gmt$ag_name))]

saveRDS(list("landscape_fit" = lndscp_fits, "gmt_data" = titerdata), paste0("data/landscape_fit/lod2_biv_boost_fits_ags_", fit_ags,"_hamster_merged_sub_adj.rds"))


# Add impulses
titerdata %>%
  select(!sr_name) %>%
  unique() %>%
  # manually set GMT's that are lower than that to LOD2
  mutate(gmt = ifelse(gmt < log2(0.8), log2(0.8), gmt))-> gmt_data


sr_group_gmt %>%
  mutate(titer = 2^logtiter*10) %>%
  select(sr_group, logtiter, ag_name, titer) %>%
  mutate(gmt = logtiter) %>%
  unique() %>%
  # manually set GMT's that are lower than that to LOD2
  mutate(gmt = ifelse(gmt < log2(0.8), log2(0.8), gmt))-> gmt_data

# angle for html page# angle for html`srGroups<-`() page
angle <- list(
  rotation = c(-1.5, 0, 0),# c(-1.4592, 0.0045, -0.0144)
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.5
  # zoom = 1.1646 # higher is more zoomed out
)


# get correct indices
sr_group_order <- titertables_groups$sr_group

big_gmts <- c(1:2)


lndscp_list <- list()
data3js <- base_plot_data3js(map, lndscp_fits, sub_ags, lims, sub_ags, alternative_ba5 = FALSE, lower_z = log2(0.8), upper_z = log2(1.6)+9, z_by = 1, label_only_subset = TRUE)

for(n in c(0, 2, 4, 6)){
  lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups[big_gmts+n,], lndscp_fits[big_gmts+n], map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors,
                                          hide_buttons = FALSE)
  
  lndscp <-r3js(
    lndscp_3js,
    rotation = angle$rotation,
    zoom = angle$zoom
  )
  lndscp_list[[paste(sr_group_order[n+1], "GMT")]] <- lndscp
}

# make idvl landscape fits
individual_lndscp_fits <- list()

for(l in 1:length(lndscp_fits)){
  
  temp_lndscp <- lndscp_fits[[l]]
  
  temp_table <- temp_lndscp$titers
  
  cone_slope <- temp_lndscp$cone$cone_slope
  
  individual_lndscp_fits[[l]] <- lapply(1:nrow(temp_table), function(i){
    
    cone_coords <- temp_lndscp$cone$cone_coords[i,]
    cone_height <- temp_lndscp$cone$cone_heights[i]
    
    ablandscape.fit(
      titers    = temp_table[i, ags_to_fit_lndscp, drop = F],
      acmap     = map,
      bandwidth = 1,
      degree    = 1,
      method = "cone",
      error.sd = 1,
      control = list(
        optimise.cone.slope = FALSE,
        optimise.cone.coords = FALSE,
        start.cone.coords = cone_coords,
        start.cone.slope = cone_slope
      )
    )
    
  }
  )
  
}


individual_sr_colors <- list()
for(table_nr in 1:length(titertables)){
  individual_sr_colors[[table_nr]] <- sapply(rownames(titertables[[table_nr]]), function(sr){
    "grey80"
  })
}


# plot landscapes
for(srg in 1:length(unique(titertables_groups$sr_group))){
  
  target_rows <- srg
  lndscp_fits_t <- lndscp_fits[target_rows]
  titertables_groups_t <- titertables_groups[target_rows,]
  
  # plot idvl landscapes
  lndscp_fits_idvl <- individual_lndscp_fits[[target_rows]]
  lndscp_colors_idvl <- individual_sr_colors[[target_rows]]
  
  opacity_val <- 0.1
  data3js_idvl <- plot_idvl_landscapes_from_list(data3js, lndscp_fits_idvl, lndscp_colors_idvl)
  
  opacity_val <- 1
  lndscp_3js <- plot_landscapes_from_list(data3js_idvl, titertables_groups_t, lndscp_fits_t, map, gmt_data, agNames(map), agNames(map), lndscp_colors = sr_colors,
                                          hide_buttons = FALSE)
  
  lndscp <-r3js(
    lndscp_3js,
    rotation = angle$rotation,
    zoom = angle$zoom
  )
  
   srg_n <- titertables_groups$sr_group[target_rows]
   srg_n <- gsub("/", "_", srg_n)
   srg_n <- gsub(" ", "", srg_n)
  
  lndscp_list[[srg_n]] <- lndscp
}


saveRDS(lndscp_list, paste0("data/landscape_fit/lod2_biv_boosts_gmt_idvl_ags_",fit_ags,"_hamster_merged_sub_adj.rds"))

