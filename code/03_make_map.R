# load required packages
rm(list = ls())
library(Racmacs)
library(stringr)

# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")

# load data
mapColors <- read.csv(file = './data/metadata/map-colors_new.csv', row.names = 'Antigen', header = TRUE)
srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors_new.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
mapColors[srGroup_colors$SerumGroup, "Color"] <- srGroup_colors$Color 

alignment_map <- read.acmap("./data/maps/roessler_netzl_et_al2023.ace")

# load titer table
titer_table <- read.titerTable("./data/titer_data/titer_table_map_samples.csv")
# ------------------------------------------ FULL MAP --------------------------------------
map <- makeMap(titer_table, nOptimisations = 1500, dilution_stepsize = 0, options = list(ignore_disconnected = TRUE))

agFill(map) <- mapColors[agNames(map),]
agGroups(map) <- agNames(map)

sr_groups <- unlist(lapply(srNames(map), function(x) str_split(x, "_")[[1]][2]))
srGroups(map) <- factor(sr_groups, levels = srGroup_colors$SerumGroup)
srOutline(map) <- mapColors[as.character(srGroups(map)),]

map <- apply_style(map)

map <- realignMap(map, alignment_map)

save.acmap(map, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure.ace")


#----------------- Make hamster map
hamster_table <- read.csv("./data/titer_data/250128_hamster_table_ba2_corrected.csv")
hamster_sr_names <- hamster_table %>%
  pull(X)

hamster_table <- t(hamster_table[,-1])
colnames(hamster_table) <- hamster_sr_names

rownames(hamster_table) <- c("D614G","B.1.1.7","B.1.1.7+E484K","B.1.351","B.1.617.2","BA.1", "BA.2", "BA.5.3.2", "XBB.1.5", "BA.2.86", "JN.1")
hamster_map <- makeMap(hamster_table, nOptimisations = 1500, dilution_stepsize = 0, options = list(ignore_disconnected = TRUE))

sr_groups_hamster <- c("d614g" = "D614G conv. (Hamster)", "delta" = "delta conv. (Hamster)", "ba.1" = "BA.1 conv. (Hamster)",
                       "ba.5" = "BA.5 conv. (Hamster)", "xbb.1.5" = "XBB.1.5 conv. (Hamster)")

srGroups(hamster_map) <- sapply(srNames(hamster_map), function(x){
  srg <- strsplit(x, " ")[[1]][1]
  sr_groups_hamster[tolower(srg)]
})
srGroups(hamster_map)
srGroups(hamster_map) <- factor(srGroups(hamster_map), levels = unname(sr_groups_hamster))

agFill(hamster_map) <- mapColors[agNames(hamster_map), "Color"]
srOutline(hamster_map) <- mapColors[as.character(srGroups(hamster_map)), "Color"]

srShape(hamster_map) <- "TRIANGLE"

hamster_map <- apply_style(hamster_map)
hamster_map <- realignMap(hamster_map, map)
save.acmap(hamster_map, "./data/maps/hamster_map.ace")

# make variants that are not in hamster map more transparent
highlight_variants <- agNames(hamster_map)

transparent_variants <- agNames(map)[!(agNames(map) %in% highlight_variants)]


agFill(map)[agNames(map) %in% transparent_variants] <- adjustcolor(agFill(map)[agNames(map) %in% transparent_variants], alpha.f = 0.25)
agOutline(map)[agNames(map) %in% transparent_variants] <- adjustcolor(agOutline(map)[agNames(map) %in% transparent_variants], alpha.f = 0.25)
agSize(map)[agNames(map) %in% transparent_variants] <- 11

save.acmap(map, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure.ace")

# merge Maps now
merged <- mergeMaps(list(map, hamster_map))

dilutionStepsize(merged) <- 0
merged <- optimize_and_realign_map(merged, map, 2500, 2, option_list = list(dim_annealing = TRUE, ignore_disconnected = TRUE))

ptDrawingOrder(merged) <- c(174:157, ptDrawingOrder(map))

save.acmap(merged, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure-merged2.ace")

