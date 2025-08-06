library(chromote)
library(htmlwidgets)

lndscp_list <-  readRDS(file.path("data", "landscape_fit", "biv_boosts_gmt_idvl_ags_sub_hamster_merged_sub_adj.rds"))

full_path <- file.path(getwd(), "figures", "landscapes", "202506_manuscript")
dir.create(full_path)

path_to_figures <- full_path

screenshot_landscapes <- function(lndscp_list, path_to_figures){
  
  
  target_scps <- names(lndscp_list)[grepl("GMT", names(lndscp_list))]
  for(name_scp in target_scps){
 
    lndscp <- lndscp_list[[name_scp]]
    file_path <- file.path(path_to_figures, paste0(name_scp, ".html"))
    saveWidget(lndscp, file_path, selfcontained = FALSE)
    
    b <- chromote::ChromoteSession$new()
    b$Page$navigate(paste0("file://", normalizePath(file_path, winslash = "/")))
    Sys.sleep(0.3)
    
    b$set_viewport_size(width = 1400,
                        height =1000,
                        zoom = 1)
    
    file_save <- gsub(".html", ".png", file_path)
    
    b$screenshot(
      filename = file_save,
      cliprect = c(50, 220, 1100, 600), #x, y, width, height
      scale = 1,
      show = FALSE,
      wait_ = TRUE
    )
    
    b$close()
    
    file.remove(file_path)
    
    unlink(gsub(".html", "_files", file_path))
  }
  
  
}


screenshot_landscapes(lndscp_list, full_path)


  
