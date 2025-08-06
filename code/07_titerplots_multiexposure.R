rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
library(titertools)
library(scales)
library(ggsci)

set.seed(100)

# define titerplot function
source("functions/long_map_info.R")
source("functions/titer_lineplot_functions.R")
source("./functions/sr_group_color_functions.R")

calc_titertools_gmt <- FALSE

sr_group_colors <- read.csv(file = "./data/metadata/sr_group_colors_new.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",",
                            row.names = "SerumGroup")

map <- read.acmap("./data/maps/map-merged-adj_colour.ace")

sub_ags <- c("D614G", "BA.1", "BA.2", "BA.5.3.2", "XBB.1.5", "JG.3", "BA.2.86", "JN.1")

ag_order <- sub_ags
titer_types <- c("*" = 0, "<16" = 2)

ymax <- 11

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
  mutate(sr_name = paste(sample.ID, sr_group, sr_info, sep = "_"),
         titertype = ifelse(titer %in% c("*", "<16"), titer_types[titer], 1),
         logtiter = ifelse(titer == "<16", log2(0.8), log2(as.numeric(titer)/10)))-> map_long


# adjust titers based on map reactivities
ag_reacts <- agReactivityAdjustments(map)
names(ag_reacts) <- agNames(map)

map_long %>%
  mutate(logtiter_adj = logtiter + ag_reacts[ag_name],
         titer = ifelse(titer == "<16", paste0("<", round(2^logtiter_adj*20)), as.character(2^logtiter_adj*10)),
         logtiter = logtiter_adj,
         sr_group = paste0(sr_group, "_adj")) %>%
  select(!logtiter_adj) -> map_long_adj

map_long <- rbind(map_long, map_long_adj)

## Here we are going to calculate the gmts and sr reactivity effects
# to pass to the function, we will create a map titer table
if(calc_titertools_gmt){
  map_long[!(grepl("adj", map_long$sr_group)),] %>%
    select(sr_name, ag_name, titer) %>%
    pivot_wider(names_from = "sr_name", values_from = "titer") %>%
    column_to_rownames("ag_name") -> titer_table_unadj
  
  gmt_map <- make.acmap(titer_table = titer_table_unadj)
  srGroups(gmt_map) <- sapply(srNames(gmt_map), function(x){
    strsplit(x, "_")[[1]][2]
  })
  
  sr_group_gmt_unadj <- calc_gmt_sr_effects_all_srgs(gmt_map, map_long, 16)
  write.csv(sr_group_gmt_unadj$ag_means, "data/titer_data/titertools_sr_group_gmt_multi_unadj.csv")
  write.csv(sr_group_gmt_unadj$sr_effects, "data/titer_data/titertools_sr_effects_multi_unadj.csv")
  sr_group_gmt_unadj <- read.csv("data/titer_data/titertools_sr_group_gmt_multi_unadj.csv") %>%
    select(!X)
  
  # do the same for the reactivity adjusted
  map_long[(grepl("adj", map_long$sr_group)),] %>%
    select(sr_name, ag_name, titer) %>%
    pivot_wider(names_from = "sr_name", values_from = "titer") %>%
    column_to_rownames("ag_name") -> titer_table_unadj
  
  gmt_map <- make.acmap(titer_table = titer_table_unadj)
  srGroups(gmt_map) <- sapply(srNames(gmt_map), function(x){
    paste0(strsplit(x, "_")[[1]][2], "_adj")
  })
  
  sr_group_gmt_adj <- calc_gmt_sr_effects_all_srgs(gmt_map, map_long, 16)
  write.csv(sr_group_gmt_adj$ag_means, "data/titer_data/titertools_sr_group_gmt_multi_adj.csv")
  write.csv(sr_group_gmt_adj$sr_effects, "data/titer_data/titertools_sr_effects_multi_adj.csv")
  sr_group_gmt_adj <- read.csv("data/titer_data/titertools_sr_group_gmt_multi_adj.csv") %>%
    select(!X)
} else {
  
  sr_group_gmt <- map_long %>%
    group_by(sr_group, ag_name) %>%
    reframe(all_below_thresh = length(titer[titer == "<16"]) == length(titer),
            lower = Rmisc::CI(logtiter)["lower"],
              upper = Rmisc::CI(logtiter)["upper"],
              logtiter = Rmisc::CI(logtiter)["mean"])
  
  sr_group_gmt_adj <- sr_group_gmt[grepl("adj", sr_group_gmt$sr_group),] %>%
    mutate(sr_group = gsub("_adj", "", sr_group))
  
  write.csv(sr_group_gmt_adj, "data/titer_data/lod2_sr_group_gmt_multi_adj.csv")
  
  sr_group_gmt_unadj <- sr_group_gmt[!grepl("adj", sr_group_gmt$sr_group),]
  
  write.csv(sr_group_gmt_unadj, "data/titer_data/lod2_sr_group_gmt_multi_unadj.csv")
  
}

# mutate(sr_group = gsub("WT", "Anc. virus", sr_group))
# gmts_sr_effects <- calc_gmt_sr_effects_all_srgs(map, data_long, 16)

ag_order <- sub_ags
target_sr_groups <- unique(map_long$sr_group)
fac_levels <- c("no event first", "no event second",
                "infection first", "infection second",
                "vaccine first", "vaccine second",
                "vaccine + inf first", "vaccine + inf second")

facet_labels <- rep(c("first TP", "second TP"), 4)
names(facet_labels) <- fac_levels

titerplot20 <- do_titer_plot_fc_label(map, 2, thresh = 16, fc_label = FALSE, adj_titers = FALSE,
                                      sr_group_gmt_plotdata = sr_group_gmt_unadj %>%
                                        filter(sr_group %in% fac_levels),
                                      fc_df = NULL,
                                      target_sr_groups = fac_levels,
                                      ag_order = ag_order,
                                      sr_group_colors = sr_group_colors,
                                      color_var = "sr_group",
                                      facet_var = "sr_group",
                                      data_long = map_long %>%
                                        filter(sr_group %in% fac_levels),
                                      titertools_gmt = FALSE,
                                      titertools_fc = FALSE,
                                      show_idvls = TRUE,
                                      show_gmt_conf = TRUE,
                                      show_gmt = TRUE,
                                      facet_levels = fac_levels,
                                      ymax = ymax) + 
  facet_wrap(~sr_group,
             labeller = labeller(sr_group = facet_labels),
             ncol = 2)


# do it here for adjusted
sr_group_gmt_adj %>%
  mutate(sr_group = gsub("_adj", "", sr_group)) ->sr_group_gmt_adj

map_long_adj %>%
  mutate(sr_group = gsub("_adj", "", sr_group)) ->map_long_adj

titerplot20_adj <- do_titer_plot_fc_label(map, 2, thresh = 16, fc_label = FALSE, adj_titers = FALSE,
                                      sr_group_gmt_plotdata = sr_group_gmt_adj %>%
                                        filter(sr_group %in% fac_levels),
                                      fc_df = NULL,
                                      target_sr_groups = fac_levels,
                                      ag_order = ag_order,
                                      sr_group_colors = sr_group_colors,
                                      color_var = "sr_group",
                                      facet_var = "sr_group",
                                      data_long = map_long_adj %>%
                                        filter(sr_group %in% fac_levels),
                                      titertools_gmt = FALSE,
                                      titertools_fc = FALSE,
                                      show_idvls = TRUE,
                                      show_gmt_conf = TRUE,
                                      show_gmt = TRUE,
                                      facet_levels = fac_levels,
                                      ymax = ymax) + 
  facet_wrap(~sr_group,
             labeller = labeller(sr_group = facet_labels),
             ncol = 2)


# here fold change plot
# Need to change the sr group order
map_long$first_tp <- grepl("first", map_long$sr_group)
map_long$sample_comb <- sapply(map_long$sr_name, function(x){
  strsplit(x, "_")[[1]][1]
})

map_long %>%
  select(sample_comb, sr_group, ag_name, logtiter, first_tp) %>%
  mutate(sr_group = gsub(" first| second", "", sr_group)) %>%
  pivot_wider(names_from = "first_tp", values_from = "logtiter") %>%
  mutate(fc = `FALSE` - `TRUE`) -> fc_df

fc_df <- fc_df[!grepl("adj", fc_df$sr_group),] %>%
  mutate(sr_group = factor(sr_group, levels = c("no event", "infection",
                                                "vaccine", "vaccine + inf")))

fc_df %>%
  group_by(sr_group, ag_name) %>%
  summarize(lower = Rmisc::CI(fc)["lower"],
            upper = Rmisc::CI(fc)["upper"],
            fc = Rmisc::CI(fc)["mean"]) %>%
  mutate(sample_comb = "gmt")-> fc_gmt

plot_colors <- sr_group_colors$Color
names(plot_colors) <- rownames(sr_group_colors)

fc_labels <- c("no event" = "no event",
               "infection" = "infection",
               "vaccine" = "XBB.1.5 boost",
               "vaccine + inf" = "XBB.1.5 boost + inf")


fc_gmt %>%
  ggplot(aes(x = ag_name, y = fc, group = sample_comb, color = sr_group)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_line(data = fc_df, alpha = 0.2) + 
  geom_point(data = fc_df, alpha = 0.2) +
  geom_line() + 
  geom_pointrange(aes(ymin = lower, ymax = upper)) +
  geom_point(color = "white", size = 1) +
  scale_color_manual(values = plot_colors,
                     name = "Serum group",
                     labels = fc_labels) +
  facet_wrap(~sr_group, ncol = 1,
             labeller = labeller(sr_group = fc_labels)) + 
  scale_x_discrete(limits = sub_ags,
                   name = "Antigen variant",
                   expand = expansion(add = 0)) + 
  scale_y_continuous(labels = function(x) paste0(round(2^x, 2), "x"),
                     name = "FC from first to second TP",
                     breaks = seq(-4,8,2)) +
  theme_bw() + 
  coord_cartesian(ylim = c(-5, 8)) +
  theme(strip.background.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") -> p

for (n in sub_ags) {
  p <- p +
    annotate(
      "tile",
      x = n,
      y = -5.5,
      height = 1,
      fill = agFill(map)[grepl(n, agNames(map), fixed = TRUE)][1],
      color = NA
    ) 
}  

titerplot20 + p + plot_layout(widths = c(2, 1)) -> sr_group_plots_unadj
titerplot20_adj + p + plot_layout(widths = c(2, 1)) -> sr_group_plots_adj


# set colors for sr groups
gmt_cols <- plot_colors
fc_gmt %>%
  ggplot(aes(x = ag_name, y = fc, group = sr_group, color = sr_group )) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) + 
  geom_point(color = "white", size = 0.5,position = position_dodge(width = 0.2)) +
  scale_x_discrete(limits = sub_ags,
                   name = "Antigen variant",
                   expand = expansion(add = 0)) + 
  scale_y_continuous(labels = function(x) paste0(round(2^x, 2), "x"),
                     name = "FC from first to second TP",
                     breaks = seq(-4,8,2)) +
  scale_color_manual(values = gmt_cols,
                     name = "Serum group",
                     labels = fc_labels)+
  theme_bw() + 
  coord_cartesian(ylim = c(-5, 8)) +
  theme(strip.background.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") -> p_gmt

for (n in sub_ags) {
  p_gmt <- p_gmt +
    annotate(
      "tile",
      x = n,
      y = -5.5,
      height = 1,
      fill = agFill(map)[grepl(n, agNames(map), fixed = TRUE)][1],
      color = NA
    ) 
}  


# make here adjusted plots with sr groups
sr_group_gmt_adj %>%
  mutate(sr_name = gsub(" first| second","", sr_group),
         sr_group = paste0(gsub("vaccine |vaccine \\+ inf |no event |infection ", "", sr_group), " TP"),
         sr_name = factor(sr_name, levels = c("no event", "infection",
                                              "vaccine", "vaccine + inf"))) %>%
  mutate(titertype = 1) -> comb_gmt_adj


do_titerplot_gmt(comb_gmt_adj, map, sub_ags, 16, gmt_cols) + 
  scale_color_manual(values = plot_colors,
                     labels = fc_labels,
                     name = "Serum group") +
  theme(legend.position = "top") -> p_titer


(p_titer +  p_gmt + plot_layout(widths = c(2, 1))) / guide_area() + plot_layout(guides = "collect", heights = c(1, 0.1)) -> gmts_adj_plot
  
sr_group_gmt_unadj %>%
  mutate(sr_name = gsub(" first| second","", sr_group),
         sr_group = paste0(gsub("vaccine |vaccine \\+ inf |no event |infection ", "", sr_group), " TP"),
         sr_name = factor(sr_name, levels = c("no event", "infection",
                                                "vaccine", "vaccine + inf"))) %>%
  mutate(titertype = 1) -> comb_gmt_unadj


do_titerplot_gmt(comb_gmt_unadj, map, sub_ags, 16, gmt_cols) + 
  scale_color_manual(values = plot_colors,
                     labels = fc_labels,
                     name = "Serum group") +
  theme(legend.position = "top") -> p_titer_unadj

(p_titer_unadj +  p_gmt + plot_layout(widths = c(2, 1))) / guide_area() + plot_layout(guides = "collect", heights = c(1, 0.1)) -> gmts_unadj_plot


sr_group_plots_unadj / gmts_unadj_plot + plot_layout(heights = c(4, 1.4)) + plot_annotation(tag_levels = c("A")) -> comb_plot_unadj

sr_group_plots_adj / gmts_adj_plot + plot_layout(heights = c(4, 1.4)) + plot_annotation(tag_levels = c("A")) -> comb_plot_adj

if(calc_titertools_gmt){
  ggsave("figures/titerplots/multi_exposures_unadj.png", comb_plot_unadj, width = 7, height = 11, dpi = 300)
  
  ggsave("figures/titerplots/multi_exposures_adj.png", comb_plot_adj, width = 7, height = 11, dpi = 300)
} else {
  
  ggsave("figures/titerplots/multi_exposures_unadj_lod2.png", comb_plot_unadj, width = 7, height = 11, dpi = 300)
  
  ggsave("figures/titerplots/multi_exposures_adj_lod2.png", comb_plot_adj, width = 7, height = 11, dpi = 300)
  
  
}

