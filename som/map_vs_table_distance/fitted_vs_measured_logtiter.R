
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(ggplot2)
library(patchwork)

map <- read.acmap("./data/maps/map-merged-adj.ace")

get_dist_df <- function(map) {
  df <- data.frame("map_d" = as.numeric(as.vector(mapDistances(map))), 
                   "table_d" = as.vector(tableDistances(map)),
                   "fitted_titer" = as.vector(sapply(1:ncol(mapDistances(map)), function(x) {
                     colBases(map)[x] - as.numeric(mapDistances(map)[,x])
                   })),
                   "fitted_titer_w_below" = as.vector(sapply(1:ncol(mapDistances(map)), function(x) {
                     mapd <- as.numeric(gsub(">", "", mapDistances(map)[,x]))
                     colBases(map)[x] - mapd
                   })),
                   "measured_titer" = as.vector(adjustedLogTiterTable(map))) %>%
    mutate(below_detectable = grepl(">", table_d), 
           table_d = as.numeric(gsub(">", "", table_d)),
           difference = table_d - map_d,
           titer_diff = measured_titer - fitted_titer)
  
  return(df)
}

plot_map_vs_table_dist <- function(map) {
  
  df <- get_dist_df(map)
  
  ggplot() +
    geom_abline(intercept =0 , slope = , linetype = "dashed", color = "black") +
    geom_point(data = df %>% filter(!(below_detectable)), aes(x = map_d, y = table_d), color = unique(srOutline(map))) +
    labs(x = "Log2 map distance", y = "Log2 table distance", title = paste0(unique(srGroups(map)), " sera")) +
    scale_y_continuous(breaks = c(0:6), limits = c(0,6)) +
    scale_x_continuous(breaks = c(0:6), limits = c(0,6)) +
    
    theme_bw()-> p
  
 return(p)
}

plot_map_vs_table_titers <- function(map) {
  
  df <- get_dist_df(map)
  
  if(unique(as.character(srGroups(map))) == "mRNA1273/mRNA1273") {
    sr <- "mRNA-1273/mRNA-1273"
  }else if(unique(as.character(srGroups(map))) == "WT conv."){
    sr <- "Ancestral virus conv."
  }else if(unique(as.character(srGroups(map))) == "BA.1 conv."){
    sr <- "BA.1 omicron conv."
  }else if(unique(as.character(srGroups(map))) == "BA.2 conv."){
    sr <- "BA.2 omicron conv."
  } else if(unique(as.character(srGroups(map))) == "alpha/alpha+E484K conv."){
    sr <- "alpha conv."
  } else if(unique(as.character(srGroups(map))) == "BA.5 conv."){
    sr <- "BA.5 omicron conv."
  }  else {
    sr <- unique(as.character(srGroups(map)))
  }
  ggplot() +
    geom_abline(intercept =0 , slope = , linetype = "dashed", color = "black") +
    geom_point(data = df %>% filter(!(below_detectable)), aes(x = fitted_titer, y = measured_titer), color = unique(srOutline(map)), alpha = 0.7) +
    labs(title = paste0(sr, " sera")) +
    scale_y_continuous(labels = function(x) round(2^x*10),
                       breaks = c(0:9),
                       limits = c(0.5, 9),
                       name = "Measured titers") +
    scale_x_continuous(labels = function(x) round(2^x*10),
                       breaks = c(0:9),
                       limits = c(0.5, 9),
                       name = "Fitted titers") +
    theme_bw() + 
    theme(plot.title = element_text(size = 8),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 5))-> p
  
  return(p)
}

plot_residual_distance <- function(map) {
  
  colors_fill = c("TRUE" = "grey", "FALSE" =  unique(srOutline(map)))
  df <- get_dist_df(map)
  
  m <- round(mean(df$difference[!df$below_detectable], na.rm = T),2)
  sd <- round(sd(scale(df$difference[!df$below_detectable], scale = F), na.rm = T), 2)
  
  ggplot() +
    geom_histogram(data = df, aes(x = difference, fill= below_detectable), color = "grey90") +
    scale_fill_manual(values = colors_fill) +
    labs(x = "Measured - fitted log2 distance", y = "Count") +
    xlim(c(-4, 4)) + 
    ylim(c(0,16)) + 
    theme_bw() +
    geom_text(mapping = aes(x = 3, y = 15, label = paste('mu', "==", m)), parse = TRUE, color = "black") +
    geom_text(mapping = aes(x = 3, y = 14, label = paste('sigma', "==", sd)), parse = TRUE, color = "black") +
    
    theme(legend.position = "none") -> p
  
  return(p)
}

plot_residual_titers <- function(map) {
  
  colors_fill = c("TRUE" = "grey", "FALSE" =  unique(srOutline(map)))
  df <- get_dist_df(map)
  
  m <- round(mean(df$titer_diff[!df$below_detectable], na.rm = T),2)
  sd <- round(sd(scale(df$titer_diff[!df$below_detectable], scale = F), na.rm = T), 2)
  
  ggplot() +
    geom_histogram(data = df, aes(x = titer_diff, fill= below_detectable), color = "grey90") +
    scale_fill_manual(values = colors_fill) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
    labs(x = "Measured - fitted log2 titers", y = "Count") +
    xlim(c(-4, 4)) + 
    ylim(c(0,45)) + 
    theme_bw() +
    geom_text(mapping = aes(x = 4, y = 40, label = "Detectable titers:", hjust =1), color = "black", size = 3) +
    geom_text(mapping = aes(x = 4, y = 35, hjust =1, label = paste('mu', "==", m)), parse = TRUE, color = "black", size = 3) +
    geom_text(mapping = aes(x = 4, y = 30, hjust =1, label = paste('sigma', "==", sd)), parse = TRUE, color = "black", size = 3) +
    theme(legend.position = "none",
         axis.title = element_text(size = 8),
         
         ) -> p
  
  return(p)
} 

map_vs_table_plots <- list()
difference_hist <- list()
titer_difference_hist <- list()
titers_plots <- list()

srGroups(map)[as.character(srGroups(map)) %in% c("alpha conv.", "alpha+E484K conv.")] <- "alpha/alpha+E484K conv."


target_srgs <- c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.',
                 'BA.5 conv.', 'CK.2.1.1 conv.', 'D614G conv. (Hamster)', 'delta conv. (Hamster)', 'BA.1 conv. (Hamster)', 'BA.5 conv. (Hamster)', 'XBB.1.5 conv. (Hamster)')

labels <- data.frame(
  row.names = target_srgs,
  val = LETTERS[1:length(target_srgs)]
)

for(sr in target_srgs){
   
  map_vs_table_plots[[sr]] <- plot_map_vs_table_dist(subsetMap(map, sera = srNames(map)[as.character(srGroups(map)) == sr]))
  titers_plots[[sr]] <- plot_map_vs_table_titers(subsetMap(map, sera = srNames(map)[as.character(srGroups(map)) == sr]))
  difference_hist[[sr]] <- plot_residual_distance(subsetMap(map, sera = srNames(map)[as.character(srGroups(map)) == sr]))
  titer_difference_hist[[sr]] <- plot_residual_titers(subsetMap(map, sera = srNames(map)[as.character(srGroups(map)) == sr]))
}


ggpubr::ggarrange(plotlist = titers_plots[1:5], labels = labels$val[1:5], nrow = 1) -> titers_top
ggpubr::ggarrange(plotlist = titer_difference_hist[1:5],labels = rep(" ",5), nrow = 1) -> titers_res_top 

ggpubr::ggarrange(plotlist = titers_plots[6:10], labels = labels$val[6:10], nrow = 1) -> titers_mid
ggpubr::ggarrange(plotlist = titer_difference_hist[6:10],labels = rep(" ",5), nrow = 1) -> titers_res_mid


ggpubr::ggarrange(plotlist = titers_plots[11:15], labels = labels$val[11:15], nrow = 1) -> titers_bottom
ggpubr::ggarrange(plotlist = titer_difference_hist[11:15],labels = rep(" ", 5), nrow = 1) -> titers_res_bottom


ggpubr::ggarrange(plotlist = titers_plots[16:17], labels = labels$val[16:17], nrow = 1, ncol = 5) -> titers_final
ggpubr::ggarrange(plotlist = titer_difference_hist[16:17],labels = rep(" ", 5), nrow = 1, ncol = 5) -> titers_res_final


# titers_top / titers_res_top / titers_bottom / titers_res_bottom -> combo

ggpubr::ggarrange(titers_top, titers_res_top, titers_mid, titers_res_mid, titers_bottom, titers_res_bottom, titers_final, titers_res_final, nrow = 8, ncol = 1) -> combo
ggsave(filename = "som/map_vs_table_distance/map_vs_table_titer_plot.png", plot = combo, width=10, height=14, dpi = 300)

 
