# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
library(dplyr)
library(ggplot2)



# Read the map
map_n <- "map-merged-adj.ace"
map <- read.acmap(paste0("./data/maps/", map_n))
results <- readRDS(paste0("./som/cross_validation/cross_validation_titer_", gsub(".ace", ".rds", map_n)))

agFillScale <- function(map) {
  fills <- agFill(map)
  names(fills) <- agNames(map)
  fills
}


# Set detectable results subset
detectable_results <- filter(results, measured_titer_type == 1)


detectable_rmse <- sqrt(mean(detectable_results$residual^2))

hist(detectable_results$residual, 40, freq = F, xlab = "Measured log titer - Predicted log titer", 
       main = sprintf("Detectable titer rmse = %s", round(detectable_rmse, 2)))

# do histogram of pred - measured
mean <- round(mean(detectable_results$residual, na.rm = T),2)
sd <- round(sd(scale(detectable_results$residual, scale = F), na.rm = T), 2)

hist_diff <- ggplot(detectable_results) +
  geom_histogram(aes(x = residual), fill = "grey50", alpha = 0.8, bins = 100) +
  xlim(c(-15, 15)) +
  geom_vline(xintercept = mean, linetype = "dashed") +
  labs(x= "Measured - predicted log2 titers", y = "Count", title = paste0("Mean = ", mean, "; SD = ", sd)) +
  theme_bw()

do_hist_single_sr_group <- function(data, sr_group_spec, color) {
  
  data %>%
    filter(sr_group == sr_group_spec) -> data_sub
  
  mean <- round(mean(data_sub$residual, na.rm = T),2)
  sd <- round(sd(scale(data_sub$residual, scale = F), na.rm = T), 2)
  
  data_sub %>%
    ggplot() +
    geom_histogram(aes(x = residual), fill = color, alpha = 0.8, bins = 100) +
    geom_vline(xintercept = mean, linetype = "dashed") +
    xlim(c(0, 300)) +
    ylim(c(0, 3000)) +
    labs(x= "Measured - predicted log2 titers", y = "Count", title = paste0("Mean = ", mean, "; SD = ", sd)) +
    theme_bw() -> plot
  
  return(plot)
}


ggsave(plot = hist_diff, filename = "./som/cross_validation/histogram_residuals.png", width = 5, height = 4, dpi = 300)


# have to change it to include new antigens and sera
ag_pretty <- data.frame(
  row.names = c("D614G", "B.1.1.7", "B.1.1.7+E484K", "P.1.1", "B.1.351", "B.1.617.2", "BA.1", "BA.2", 
                'CB.1', 'BR.3', 'CH.1.1', 'BA.5.3.2', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF', "JG.3", "BA.2.86", "JN.1"),
  val = c('D614G', 'alpha', 'alpha+E484K', 'gamma', 'beta', 'delta', 'BA.1 omicron', 'BA.2 omicron',
          'CB.1', 'BR.3', 'CH.1.1', 'BA.5 omicron', 'BA.5.2.1', 'BE.1.1', 'BF.7', 'BQ.1.3', 'BQ.1.1', 'BQ.1.18', 'XBB.1', 'XBB.1.5', 'XBF', "JG.3", "BA.2.86", "JN.1")
)


sr_pretty <- data.frame(
  row.names = c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.', 'BA.5 conv.', "CK.2.1.1 conv.",
                'D614G conv. (Hamster)', 'delta conv. (Hamster)', 'BA.1 conv. (Hamster)', 'BA.5 conv. (Hamster)', 'XBB.1.5 conv. (Hamster)'),
  val = c('2xmRNA-1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"Anc. virus conv.", 'alpha conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.', 'BA.5 conv.', "CK.2.1.1 conv.",
          'D614G conv. (Hamster)', 'delta conv. (Hamster)', 'BA.1 conv. (Hamster)', 'BA.5 conv. (Hamster)', 'XBB.1.5 conv. (Hamster)')
)

detectable_results$sr_group[detectable_results$sr_group %in% c("alpha conv.", "alpha+E484K conv.")] <- 'alpha/alpha+E484K conv.'


detectable_results$ag_pretty <- factor(ag_pretty[as.character(detectable_results$ag_name),], levels = ag_pretty$val)
detectable_results$sr_pretty <- factor(sr_pretty[as.character(detectable_results$sr_group),], levels = sr_pretty$val)

# Antigen and serum group tab
        detectable_results %>%
          ggplot(
            aes(
              x = predicted_logtiter,
              y = measured_logtiter,
              color = ag_name
            )
          ) +
          labs(x = "Predicted log2 titer",
               y = "Measured log2 titer") +
          # geom_smooth() +
          geom_point(
            shape = 21,
            alpha = 1,
            size = 0.5
          ) +
          geom_abline(
            slope = 1,
            intercept = 0,
            linetype = "dashed"
          ) +
          scale_color_manual(
            # values = Racmacs:::srGroupOutline(map)
            values = agFillScale(map)
          ) +
          xlim(c(-10,10))+
          ylim(c(-10, 10))+
          facet_grid(
            rows = vars(sr_pretty),
            cols = vars(ag_pretty),
            labeller = label_wrap_gen(multi_line = TRUE,
                                      width = 10)
          ) +
          theme_bw() +
          theme(legend.position = "none",
               strip.text.x = element_text(size = 6),
               strip.text.y = element_text(size = 6))-> gp

 
ggsave(plot = gp, filename = "./som/cross_validation/scatter_pred_vs_measured.png", width = 14, height = 10, dpi = 300)
        
