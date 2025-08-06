# investigating hamster repeat data
rm(list = ls())
library(Racmacs)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)

source("./functions/map_functions.R")
# set seed
set.seed(100)

# first BA.2 titrations were much higher than homologous in XBB.1.5 conv hamsters
# Repeat titrations were lower, now adjustment of repeat magniture using repeat BA.5 and XBB.1.5 titrations

mapColors <- read.csv(file = './data/metadata/map-colors.csv', row.names = 'Antigen', header = TRUE)

hamster_data <- read_excel("data/titer_data/250128_Hamster update_w_repeat_titrations.xlsx", range = "C4:P23")
colnames(hamster_data)[1:3] <- hamster_data[1, 1:3]
hamster_data <- hamster_data[-1,]
hamster_data$sr_name <- paste(hamster_data$Sample, hamster_data$`nAB ID`, sep = "_")
hamster_data$titration <- 1

hamster_repeats <- read_excel("data/titer_data/250128_Hamster update_w_repeat_titrations.xlsx", range = "C24:P27")
hamster_repeats$sr_name <- sapply(1:3, function(x) paste(hamster_repeats[x,2], hamster_repeats[x,1], sep = "_"))
hamster_repeats$titration <- 2
colnames(hamster_repeats) <- colnames(hamster_data)

hamster_data <- rbind(hamster_data, hamster_repeats)

colnames(hamster_data)[grepl("XBB.1.5.1", colnames(hamster_data))] <- "XBB.1.5"

hamster_data %>%
  filter(`nAB ID` %in% c("H207", "H208", "H209")) %>%
  select(sr_name, titration, XBB.1.5, BA.5, BA.2, D614G, BA.1, BA.2.86, JN.1) %>%
  pivot_longer(cols = c("XBB.1.5", "BA.5", "BA.2", "D614G", "BA.1", "BA.2.86", "JN.1"), names_to = "ag_name", values_to = "titer", values_transform = as.numeric) %>%
  mutate(logtiter = log2(titer/10)) -> hamster_long

data_repeat <-  hamster_long %>% 
  filter(titration == 2)

hamster_long %>%
  filter(titration == 1) %>%
  ggplot(aes(x = ag_name, y = logtiter, group = sr_name), color = "grey20") + 
  geom_line() + 
  geom_point() + 
  geom_line(data = data_repeat, color = "skyblue") + 
  geom_point(data = data_repeat, color = "skyblue") + 
  scale_x_discrete(limits = c("D614G", "BA.1", "BA.2", "BA.5", "XBB.1.5", "BA.2.86", "JN.1")) + 
  theme_bw() +
  xlab("Variant") + 
  ylab("Log2(Titer/10)") + 
  ylim(c(-1,8)) + 
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = log2(16/10),
    fill = "grey50",
    color = NA,
    alpha = 0.3
  ) -> p_raw


hamster_data %>%
  select(sr_name, titration, XBB.1.5, BA.5, BA.2) %>%
  filter(!is.na(XBB.1.5)) %>%
  pivot_longer(cols = c("XBB.1.5", "BA.5", "BA.2"), names_to = "ag_name", values_to = "titer", values_transform = as.numeric) %>%
  mutate(logtiter = log2(titer/10)) -> repeat_long

repeat_wide <- hamster_long %>% 
  select(!titer) %>%
  pivot_wider(names_from = titration, values_from = logtiter) %>%
  mutate(ba2 = ag_name == "BA.2")

plot_fills <- mapColors[c("BA.2", "BA.5.3.2", "XBB.1.5"), "Color"]
names(plot_fills) <- c("BA.2", "BA.5", "XBB.1.5")

repeat_wide %>%
  ggplot(aes(x = `2`, y = `1`, color = ba2)) + 
  geom_point(aes(fill = ag_name), shape = 21, size = 3) + 
  geom_smooth(method = "lm") + 
  ggpubr::stat_regline_equation() +
  ylim(c(1, 8)) + 
  ylab("Log2(Titer/10) Repeat 1") + 
  xlab("Log2(Titer/10) Repeat 2") + 
  xlim(c(1,8)) + 
  scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "tomato"),
                     name = "BA.2") +
  scale_fill_manual(values = plot_fills, name = "Variant") +
  theme_bw() -> point_raw


repeat_wide %>%
  mutate(`1` = ifelse(ba2, `1` - 1.8, `1`)) %>%
  ggplot(aes(x = `2`, y = `1`, color = ba2)) + 
  geom_point(aes(fill = ag_name), shape = 21, size = 3) + 
  geom_smooth(method = "lm") + 
  ggpubr::stat_regline_equation() +
  ylim(c(1, 8)) + 
  ylab("Log2(Titer/10) Repeat 1") + 
  xlab("Log2(Titer/10) Repeat 2") + 
  xlim(c(1,8)) + 
  scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "tomato"),
                     name = "BA.2") +
  scale_fill_manual(values = plot_fills, name = "Variant") +
  theme_bw() -> point_corrected
# transform BA.2 titer accordingly

repeat_wide %>%
  mutate(`1` = ifelse(ba2, `1` - 1.8, `1`)) %>%
  filter(ba2) %>%
  mutate(titer = 2^`1`*10)-> transformed_ba2

hamster_corrected <- hamster_data %>%
  filter(titration == 1)

# change titers here and correct them
hamster_corrected[grep("XBB.1.5", hamster_corrected$Sample), "BA.2"] <- transformed_ba2 %>% pull(titer) %>%
  as.character


hamster_corrected %>%
  filter(`nAB ID` %in% c("H207", "H208", "H209")) %>%
  select(sr_name, titration, XBB.1.5, BA.5, BA.2, D614G, BA.1, BA.2.86, JN.1) %>%
  pivot_longer(cols = c("XBB.1.5", "BA.5", "BA.2", "D614G", "BA.1", "BA.2.86", "JN.1"), names_to = "ag_name", values_to = "titer", values_transform = as.numeric) %>%
  mutate(logtiter = log2(titer/10)) -> hamster_long_corrected


hamster_long %>%
  filter(titration == 1) %>%
  ggplot(aes(x = ag_name, y = logtiter, group = sr_name), color = "grey20") + 
  geom_line(data = hamster_long_corrected, color = "skyblue") + 
  geom_point(data = hamster_long_corrected, color = "skyblue") + 
  geom_line() + 
  geom_point() + 
  scale_x_discrete(limits = c("D614G", "BA.1", "BA.2", "BA.5", "XBB.1.5", "BA.2.86", "JN.1")) + 
  theme_bw() + 
  xlab("Variant") + 
  ylab("Log2(Titer/10)") + 
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = log2(16/10),
    fill = "grey50",
    color = NA,
    alpha = 0.3
  ) +
  ylim(c(-1,8))-> p_corrected


(p_raw + point_raw) / (p_corrected + point_corrected) + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") -> comb_plot


ggsave("som/hamster_repeat_inv/repeat_correction.png", comb_plot, dpi = 300, width = 10, height = 6)

## create hamster table here
hamster_table <- hamster_corrected %>%
  select(colnames(hamster_corrected)[4:(ncol(hamster_corrected)-1)]) %>%
  column_to_rownames("sr_name")

hamster_table <- as.matrix(hamster_table)
# set <16 and "*"
hamster_table[as.numeric(hamster_table) <16] <- "<16"
hamster_table[is.na(hamster_table)] <- "*"

write.csv(hamster_table, "./data/titer_data/250128_hamster_table_ba2_corrected.csv")
