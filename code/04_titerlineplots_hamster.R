rm(list = ls())
library(Racmacs)
library(tidyverse)
library(patchwork)
library(titertools)

# define titerplot function
source("functions/long_map_info.R")
source("functions/titer_lineplot_functions.R")
source("./functions/sr_group_color_functions.R")

calc_titertools_gmt <- TRUE

sr_group_colors <- read.csv(file = "./data/metadata/sr_group_colors_new.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",",
                            row.names = "SerumGroup")

map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-merged2.ace")


titer_types <- c("*" = 0, "<16" = 2)

ymax <- 11

ag_pretty <- data.frame(
  row.names = c("D614G", "B.1.1.7", "B.1.1.7+E484K", "P.1.1", "B.1.351", "B.1.617.2", "BA.1", "BA.2", 
                'CB.1', 'BR.3', 'CH.1.1', 'BA.5.3.2', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF',
                "JG.3", "BA.2.86", "JN.1"),
  val = c('D614G', 'alpha', 'alpha+E484K', 'gamma', 'beta', 'delta', 'BA.1', 'BA.2',
          'CB.1', 'BR.3', 'CH.1.1', 'BA.5', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF', "JG.3", "BA.2.86", "JN.1")
)

agNames(map) <- ag_pretty[agNames(map),1]

data_long <- long_map_info(map) %>%
  mutate(sr_group = gsub("WT", "Anc. virus", sr_group))

if(calc_titertools_gmt){
  
  if(file.exists("data/titer_data/titertools_sr_group_gmt_single_w_sr_effects.csv")){
    sr_group_gmt <- read.csv("data/titer_data/titertools_sr_group_gmt_single_w_sr_effects.csv")
  } else {
    sr_group_gmt_comb <- calc_gmt_sr_effects_all_srgs(map, data_long, 16)
    write.csv(sr_group_gmt_comb$ag_means, "data/titer_data/titertools_sr_group_gmt_single_w_sr_effects.csv", row.names = FALSE)
    write.csv(sr_group_gmt_comb$sr_effects, "data/titer_data/titertools_sr_effects_single.csv", row.names = FALSE)
    sr_group_gmt <- sr_group_gmt_comb$ag_means
  }
  
  
}

map_adj <- read.acmap("data/maps/map-merged-adj.ace")
agNames(map_adj) <- ag_pretty[agNames(map_adj),1]

data_long_adj <- long_map_info(map_adj) %>%
  mutate(sr_group = gsub("WT", "Anc. virus", sr_group))

if(calc_titertools_gmt){
  
  if(file.exists("data/titer_data/titertools_sr_group_gmt_single_w_sr_effects_adj.csv")){
    sr_group_gmt_adj <- read.csv("data/titer_data/titertools_sr_group_gmt_single_w_sr_effects_adj.csv")
  } else {
    sr_group_gmt_comb_adj <- calc_gmt_sr_effects_all_srgs(map_adj, data_long_adj, 16)
    write.csv(sr_group_gmt_comb_adj$ag_means, "data/titer_data/titertools_sr_group_gmt_single_w_sr_effects_adj.csv", row.names = FALSE)
    write.csv(sr_group_gmt_comb_adj$sr_effects, "data/titer_data/titertools_sr_effects_single_adj.csv", row.names = FALSE)
    sr_group_gmt_adj <- sr_group_gmt_comb_adj$ag_means
  }
  
  
}


fac_levels <- c('mRNA1273/mRNA1273','BNT/BNT', 'AZ/BNT','AZ/AZ', 'Anc. virus conv.',
                'alpha/alpha+E484K conv.', 'alpha conv.', 'alpha+E484K conv.',
                'beta conv.', 'delta conv.',
                'BA.1 conv.', 'BA.2 conv.', 'BA.5 conv.', 
                'CK.2.1.1 conv.', 
                "D614G conv. (Hamster)" ,"delta conv. (Hamster)" ,  "BA.1 conv. (Hamster)" ,   "BA.5 conv. (Hamster)" , "XBB.1.5 conv. (Hamster)")


# show ba.2.86 and JN.1 titers in ALL serum groups, make serum group x axis and facet by antigen
sub_data <- data_long %>%
  filter(ag_name %in% c("BA.2.86", "JN.1", "BA.1")) %>% 
  rbind(., data_long_adj %>%
          filter(ag_name %in% c("JN.1", "BA.1")) %>%
          mutate(ag_name = paste(ag_name, "adj."))) %>%
  mutate(sr_group = factor(sr_group, levels = fac_levels)) %>%
  mutate(logtiter = logtiter_adjusted,
         logtiter = ifelse(logtiter < log2(1.6), log2(0.8), logtiter),
         logtiter = ifelse(ag_name == "BA.1 adj." & logtiter < log2(3.64895), log2(3.64895/2), logtiter))


# calculate gmts
gmt_sub <- sub_data %>%
  filter(!is.na(logtiter)) %>%
  group_by(ag_name, sr_group) %>%
  summarize(lower = Rmisc::CI(logtiter)["lower"],
         upper = Rmisc::CI(logtiter)["upper"],
         logtiter = Rmisc::CI(logtiter)["mean"]) %>%
  ungroup() %>%
  mutate(logtiter = ifelse(logtiter < log2(1.6), log2(0.8), logtiter),
         sr_name = paste("GMT", sr_group))

# gmt_sub <- sr_group_gmt %>%
#   filter(ag_name %in% c("BA.2.86", "JN.1", "BA.1")) %>%
#   rbind(., sr_group_gmt_adj %>%
#           filter(ag_name %in% c("JN.1", "BA.1")) %>%
#           mutate(ag_name = paste(ag_name, "adj."))) %>%
#   mutate(sr_group = factor(sr_group, levels = fac_levels)) %>%
#   mutate(sr_name = paste("GMT", sr_group)) %>%
#   filter(!is.na(logtiter)) %>%
#   mutate(logtiter_estimate = logtiter,
#          logtiter = ifelse(logtiter < log2(1.6), log2(0.8), logtiter))



plot_colors <- sr_group_colors$Color
names(plot_colors) <- rownames(sr_group_colors)

sub_data %>%
  filter(!ag_name %in% c("BA.1", "BA.1 adj.")) %>%
  filter(!is.na(logtiter)) %>%
  do_titer_barplot(., gmt_sub %>%
                     filter(!ag_name %in% c("BA.1", "BA.1 adj.")) %>%
                     filter(!is.na(logtiter)),
                   ymax = 7) -> jn1_p


sub_data %>%
  filter(ag_name %in% c("BA.1", "BA.1 adj.")) %>%
  filter(!is.na(logtiter)) %>%
  do_titer_barplot(., gmt_sub %>%
                     filter(ag_name %in% c("BA.1", "BA.1 adj.")) %>%
                     filter(!is.na(logtiter)),
                   ymax = 7) + 
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) -> ba1_p


jn1_p + ba1_p + plot_layout(width = c(3, 2)) -> barp_comb

ggsave("figures/titerplots/ag_sub_barplots_lod2.png", barp_comb, dpi = 300, width = 10, height = 5)

# ----------------- Titer line plots here
ag_order <- c("BA.1 adj.", "BA.1", "BA.2.86", "JN.1", "JN.1 adj.")
target_sr_groups <- unique(data_long$sr_group)

# combine alpha and alpha/e484K serum groups
sub_data$sr_group[grepl("alpha", sub_data$sr_group)] <- "alpha/alpha+E484K conv."
# combine ba.5 and ck.2.1.1
sub_data$sr_group[grepl("CK", sub_data$sr_group)] <- "BA.5 conv."

# only do sr groups with jn.1 titartions and show xbb.1.5 hamsters extra
sub_data %>%
  filter(!ag_name %in% c("BA.1", "BA.1 adj.")) %>%
  filter(!is.na(logtiter)) %>%
  pull(sr_group) %>%
  unique() -> target_srgs

target_srgs <- target_srgs[!grepl("XBB", target_srgs)]

titerplot_sub <- do_titer_plot_fc_label(map, 6, thresh = 16, fc_label = FALSE, adj_titers = FALSE,
                                            sr_group_gmt_plotdata = NULL,
                                            fc_df = NULL,
                                            target_sr_groups = target_sr_groups[target_sr_groups %in% target_srgs],
                                            ag_order = ag_order,
                                            sr_group_colors = sr_group_colors,
                                            color_var = "sr_group",
                                            facet_var = "sr_group",
                                            data_long = sub_data,
                                            titertools_gmt = FALSE,
                                            titertools_fc = FALSE,
                                            show_idvls = TRUE,
                                            show_gmt_conf = TRUE,
                                            show_gmt = TRUE,
                                            facet_levels = fac_levels,
                                            ymax = 9,
                                        do_ag_tiles = FALSE) + 
  expand_limits(x = c(0.5, 5.5)) + 
  facet_wrap(~sr_group, labeller = label_wrap_gen(width=15), ncol = 6)


ag_order_hamster <- c("D614G", "BA.1 adj.", "BA.1", "BA.2", "BA.5", "XBB.1.5", "BA.2.86", "JN.1", "JN.1 adj.")

hamster_xbb15 <- data_long %>%
  rbind(., data_long_adj %>%
          filter(ag_name %in% c("JN.1", "BA.1")) %>%
          mutate(ag_name = paste(ag_name, "adj."))) %>%
  filter(sr_group == "XBB.1.5 conv. (Hamster)") %>%
  mutate(logtiter = logtiter_adjusted,
         logtiter = ifelse(logtiter < log2(1.6), log2(0.8), logtiter),
         logtiter = ifelse(ag_name == "BA.1 adj." & logtiter < log2(3.64895), log2(0.8), logtiter))

titerplot_xbb15 <- do_titer_plot_fc_label(map, 1, thresh = 16, fc_label = FALSE, adj_titers = FALSE,
                                        sr_group_gmt_plotdata = NULL,
                                        fc_df = NULL,
                                        target_sr_groups =  "XBB.1.5 conv. (Hamster)",
                                        ag_order = ag_order_hamster,
                                        sr_group_colors = sr_group_colors,
                                        color_var = "sr_group",
                                        facet_var = "sr_group",
                                        data_long = hamster_xbb15,
                                        titertools_gmt = FALSE,
                                        titertools_fc = FALSE,
                                        show_idvls = TRUE,
                                        show_gmt_conf = TRUE,
                                        show_gmt = TRUE,
                                        facet_levels = fac_levels,
                                        ymax = 9,
                                        do_ag_tiles = FALSE) + 
  expand_limits(x = c(0.5, 5.5))

titerplot_sub + titerplot_xbb15 + plot_layout(widths = c(3, 1)) -> comb_sub

ggsave("figures/titerplots/single_exposures_sub_lines.png", comb_sub, dpi = 300, width = 12, height = 5)

# -------------- SOM Figures for full plots
# GMT is set to LOD/2 here

ag_order <- ag_pretty$val[ag_pretty$val %in% agNames(map)]
target_sr_groups <- unique(data_long$sr_group)


# no hamster plot
titerplot20_human <- do_titer_plot_fc_label(map, 4, thresh = 16, fc_label = FALSE, adj_titers = FALSE,
                                      sr_group_gmt_plotdata = NULL,
                                      fc_df = NULL,
                                      target_sr_groups = target_sr_groups[!grepl("amster", target_sr_groups)],
                                      ag_order = ag_order,
                                      sr_group_colors = sr_group_colors,
                                      color_var = "sr_group",
                                      facet_var = "sr_group",
                                      data_long = data_long,
                                      titertools_gmt = FALSE,
                                      titertools_fc = FALSE,
                                      show_idvls = TRUE,
                                      show_gmt_conf = TRUE,
                                      show_gmt = TRUE,
                                      facet_levels = fac_levels,
                                      ymax = ymax)


hamster_ags <- data_long %>%
  filter(sr_group == "D614G conv. (Hamster)") %>%
  filter(titer != "*") %>%
  pull(ag_name) %>%
  unique()

titerplot20_hamster <- do_titer_plot_fc_label(map, 5, thresh = 16, fc_label = FALSE, adj_titers = FALSE,
                                            sr_group_gmt_plotdata = NULL,
                                            fc_df = NULL,
                                            target_sr_groups = target_sr_groups[grepl("amster", target_sr_groups)],
                                            ag_order = hamster_ags,
                                            sr_group_colors = sr_group_colors,
                                            color_var = "sr_group",
                                            facet_var = "sr_group",
                                            data_long = data_long,
                                            titertools_gmt = FALSE,
                                            titertools_fc = FALSE,
                                            show_idvls = TRUE,
                                            show_gmt_conf = TRUE,
                                            show_gmt = TRUE,
                                            facet_levels = fac_levels,
                                            ymax = ymax)



titerplot20_human / titerplot20_hamster + plot_layout(heights = c(4, 1)) + plot_annotation(tag_levels = "A") -> comb_unadj

ggsave("figures/titerplots/single_exposures_unadj.png", comb_unadj, dpi = 300, width = 14, height = 12)


# same for adjusted titers

titerplot20_adj <- do_titer_plot_fc_label(map_adj, 4, thresh = 16, fc_label = FALSE, adj_titers = TRUE,
                                      sr_group_gmt_plotdata = NULL,
                                      fc_df = NULL,
                                      target_sr_groups = target_sr_groups[!grepl("amster", target_sr_groups)],
                                      ag_order = ag_order,
                                      sr_group_colors = sr_group_colors,
                                      color_var = "sr_group",
                                      facet_var = "sr_group",
                                      data_long = data_long_adj,
                                      titertools_gmt = FALSE,
                                      titertools_fc = FALSE,
                                      show_idvls = TRUE,
                                      show_gmt_conf = TRUE,
                                      show_gmt = TRUE,
                                      facet_levels = fac_levels,
                                      ymax = ymax)

titerplot20_adj_hamster <- do_titer_plot_fc_label(map_adj, 5, thresh = 16, fc_label = FALSE, adj_titers = TRUE,
                                          sr_group_gmt_plotdata = NULL,
                                          fc_df = NULL,
                                          target_sr_groups = target_sr_groups[grepl("amster", target_sr_groups)],
                                          ag_order = hamster_ags,
                                          sr_group_colors = sr_group_colors,
                                          color_var = "sr_group",
                                          facet_var = "sr_group",
                                          data_long = data_long_adj,
                                          titertools_gmt = FALSE,
                                          titertools_fc = FALSE,
                                          show_idvls = TRUE,
                                          show_gmt_conf = TRUE,
                                          show_gmt = TRUE,
                                          facet_levels = fac_levels,
                                          ymax = ymax)

titerplot20_adj / titerplot20_adj_hamster + plot_layout(heights = c(4, 1)) + plot_annotation(tag_levels = "A") -> comb_adj

ggsave("figures/titerplots/single_exposures_adj.png", comb_adj, dpi = 300, width = 14, height = 12)


