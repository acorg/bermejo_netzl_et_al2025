# format titer table from excel sheet
rm(list = ls())
library(tidyverse)
library(stringr)

# vaccine status for sheet I
vacc_stat <- c("n" = "",
               "1" = "ChAdOx1-S/ChAdOx1-S", 
               "2" = "ChAdOx1-S/BNT162b2",
               "3" = "BNT162b2/BNT162b2",
               "4" = "mRNA-1273/mRNA-1273",
               "d" = "divers",
               "-" = "")

# function for sr_info

sr_info_function <- function(x, kimpel1) {
  stat <- kimpel1[[x, "Immune status"]]
  vacc <- kimpel1[[x,4]]
  var <- kimpel1[[x, "Variant"]]

  if(var == "firstwave"){
    var <- "WT"
  }

  if(stat == "convalescent" | stat == "unvaccinated") {
    res <- paste0(var, " conv.")
    
  } else if(stat == "3x vaccinated"){
    res <- vacc
  } else if(str_detect(stat, "con/vacc")) {
    
    if(str_detect(vacc, "d")) {
      res <- paste0(var, "+", substring(vacc, 4, nchar(vacc)-1))
    } else {
      res <- paste0(var, "+", vacc_stat[vacc])
    }
    
  } else if(str_detect(stat, "vacc/con")) {
    
    if(str_detect(vacc, "d")) {
      res <- paste0(substring(vacc, 4, nchar(vacc)-1), "+", var)
    } else {
      res <- paste0(vacc_stat[vacc],"+", var)
    }
    
  } else {
    res <- vacc_stat[vacc]
  }
  
  return(res)
}

shorten_vacc_names <- function(sr_group) {
  sr_group <- gsub("ChAdOx1-S", "AZ", sr_group)
  sr_group <- gsub("ChAdOx1", "AZ", sr_group)
  sr_group <- gsub("BNT162b2", "BNT", sr_group)
  sr_group <- gsub("BNT162b", "BNT", sr_group)
  sr_group <- gsub("Ad26.COV2.S", "J\\&J", sr_group)
  sr_group <- gsub("mRNA\\-1273", "mRNA1273", sr_group)
  
  sr_group <- gsub("B.1.617.2", "delta", sr_group)
  sr_group <- gsub("B.1.351", "beta", sr_group)
  sr_group <- gsub("B.1.1.7\\/B.1.1.7\\+E484K", "alpha\\/alpha\\+E484K", sr_group)
  sr_group <- gsub("alpha-E484K", "alpha+E484K", sr_group)
  
  return(sr_group)
}

set_threshold <- function(tab, thresh) {
  tab[as.numeric(tab) < as.numeric(thresh)] <- paste0("<", thresh)
  tab[is.na(tab)] <- "*"
  
  return(tab)
}


# =============================== Omicron I sheet, published in NEJM
kimpel1 <- readxl::read_excel("./data/titer_data/241008_map_update-2.xlsx", sheet = 1)
colnames(kimpel1) <- kimpel1[2,]
kimpel1 <- kimpel1[c(3:nrow(kimpel1)),2:ncol(kimpel1)]

# remove all NA rows
kimpel1 <- kimpel1[rowSums(is.na(kimpel1)) != ncol(kimpel1),]

kimpel1$Variant <- gsub(" ", "", kimpel1$Variant)
kimpel1$Variant <- gsub("\\(BA.5\\)", "", kimpel1$Variant)


all_ags <- c("D614G", "B.1.1.7", "B.1.1.7+E484K", "B.1.351", "P.1.1", "B.1.617.2", "BA.1",
  "BA.2", "CB.1", "BR.3","CH.1.1","BA.5.3.2", "BA.5.2.1", "BE.1.1", "BF.7", "BQ.1.3", "BQ.1.1", "BQ.1.18", "XBB.1", "XBB.1.5", "XBF",
  "BA.2.86", "JG.3", "JN.1")

# get XBF data
xbf_data <- readxl::read_excel("./data/titer_data/Archive/2304/230411_Omicron V_update.xlsx", sheet = 1, skip = 8) %>%
  select(`sample ID`, XBF)

kimpel1$XBF <- xbf_data$XBF[match(kimpel1$`sample ID`, xbf_data$`sample ID`)]
# add sample info
# save sample ID, then sr group, then sr info including time since dose
# for now keep old samples from duplicates as there are some differing titers and the new ones were not titrated against
# all variants, so this could lead to distotion
kimpel1$`sample ID` <- sapply(kimpel1$`sample ID`, function(x) strsplit(x, " ")[[1]][1])

kimpel1 <- kimpel1 %>%
  group_by(`sample ID`) %>%
  mutate(sample_nr = 1:length(Variant),
         sample_id_nr = paste0(`sample ID`,':', sample_nr)) %>%
  ungroup() %>%
  mutate(
    sr_info = unlist(lapply(1:nrow(kimpel1), function(x) sr_info_function(x, kimpel1)
  )),
  sr_group = shorten_vacc_names(sr_info),
  sr_info_complete = paste(`sample ID`,sr_group, sr_info, sep = "_")
  )

# set serum_info_column that contains all info for map
# pivot wider to titer table
kimpel1 %>%
  filter(sample_nr == 1) %>%
  column_to_rownames("sr_info_complete") %>%
  select(all_of(all_ags[all_ags %in% colnames(kimpel1)])) -> kimpel1_table

kimpel1_table[is.na(kimpel1_table)] <- "*"
# add lower threshold at the end
kimpel1_table <- set_threshold(as.matrix(kimpel1_table), 16)

kimpel1_table <- t(kimpel1_table)
# remove samples from titer table 
kimpel1_table <- kimpel1_table[, !(grepl("C711|C860|G776|G780", colnames(kimpel1_table)))]

# remove double samples 
kimpel1_table <- kimpel1_table[, !(grepl("G935|G803|G898", colnames(kimpel1_table)))]


write.csv(kimpel1_table, "./data/titer_data/titer_table_map_samples.csv")

#======================================  Hybrid Immunity

kimpel <- readxl::read_excel("./data/titer_data/241009_Omicron VIII_data.xlsx", sheet = 1, range = "A1:P47")

colnames(kimpel)[1:2] <- c("Vaccination", "sample ID")

kimpel %>%
  mutate(Vaccination = gsub(0, NA, Vaccination)) %>%
  fill(Vaccination) %>%
  mutate(sr_group = Vaccination,
         sr_group = ifelse(grepl("kein", `biv boost`), "no XBB Boost", Vaccination),
         prior_inf = ifelse(`N ELISA` >= 0.9, "inf", "nonInf"),
         latest_boost = ifelse(grepl("XBB", sr_group), "XBB", ifelse(grepl("BA.1", `biv boost`), "BA.1", "BA.5")),
         sr_info = paste0(sr_group, "_", prior_inf, "_",latest_boost, "_", `Days since last vacc`, '_', Sex, '_', Age),
         sample_id_nr = `sample ID`,
         sample_nr = 1) -> kimpel

kimpel$`sample ID` <- sapply(kimpel$`sample ID`, function(x) strsplit(x, " ")[[1]][1])

n_samples_biv <- nrow(kimpel %>% 
                        filter(sr_group == "biv. Boost"))

kimpel$`sample ID`[1:n_samples_biv] <- sapply(c(1:n_samples_biv), function(x){
  
  paste(kimpel$`sample ID`[x], kimpel$`sample ID`[x+n_samples_biv], sep = "/")
  
})

kimpel$`sample ID`[(n_samples_biv+1) : nrow(kimpel)] <- kimpel$`sample ID`[1: n_samples_biv]
kimpel$sample_nr[(n_samples_biv+1) : nrow(kimpel)] <-2

colnames(kimpel) <- gsub("XBB.1.5.1", "XBB.1.5", colnames(kimpel))

kimpel %>%
  mutate(sr_info_complete = paste(`sample ID`,sr_group, sr_info, sep = "_")) %>%
  select(!Vaccination) %>%
  select(!`N ELISA`) %>%
  select(!prior_inf) %>%
  plyr::rbind.fill(kimpel1 %>% select(!colnames(kimpel1)[2:4])) -> comb

comb %>%
  pivot_longer(cols = all_of(all_ags), names_to = "ag", values_to = "titer", values_transform = as.character) %>%
  mutate(titer = ifelse(is.na(titer), "*", ifelse(as.numeric(titer) < 16, "<16", titer))) %>%
  filter(!`sample ID` %in% c("C711 - leer", "C860", "G776", "G780", "G935", "G803", "G898","C711"))-> comb_long

write.csv(comb_long, "data/titer_data/biv_boost_long.csv", row.names = FALSE)

# here also duplicate samples in it. For now just for testing!
comb_long %>%
  select(sr_info_complete, ag, titer) %>%
  pivot_wider(names_from = ag, values_from = titer) %>%
  column_to_rownames("sr_info_complete") -> comb_table

comb_table <- t(as.matrix(comb_table))
idx <- which(is.na(comb_table), arr.ind=TRUE)

write.csv(comb_table, "./data/titer_data/titer_table_all_samples.csv")

