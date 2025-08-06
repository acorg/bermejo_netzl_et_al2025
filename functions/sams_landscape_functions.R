padding <- 0.1
opacity_val <- 0.9


# Functions to remove buttons
addObject3js <- function(
    data3js,
    object,
    number_of_ids = 1
){
  
  # Generate an object ID
  if(is.null(data3js$lastID)){ data3js$lastID <- 0 }
  object$ID <- max(data3js$lastID) + seq_len(number_of_ids)
  
  # If object is interactive and highlighted add a reference to itself to
  # it's highlight group by default
  if(!is.null(object$properties$interactive)){
    object$group <- object$ID
  }
  
  # Add the object to the plot data
  data3js$plot[[length(data3js$plot)+1]] <- object
  
  # Update the ID of the last object added
  data3js$lastID <- object$ID
  
  # Return the new data
  data3js
  
}

remove_buttons <- function(data3js){
  
  new_data3js = data3js
  
  new_data3js = data3js
  
  new_data3js[['lastID']] = 0
  new_data3js[['plot']] = list()
  
  N = data3js[['lastID']] 
  
  
  
  
  for (i in 1:N)
  {
    obj = data3js[['plot']][[i]]
    
    
    
    if ('toggle' %in% names(obj[['properties']])){
      obj[['properties']][['toggle']] <- NULL
    }
    
    new_data3js = addObject3js(new_data3js,obj)
    
    
  }
  
  
  
  return (new_data3js)
  
}


base_plot_data3js <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1,
                              add_border = TRUE, add_axis = TRUE, upper_z = 10, z_by = 2, add_ag_label = TRUE, lower_z = 0,
                              label_only_subset = FALSE){
  
  if(label_only_subset){
    x_coords <- c(agCoords(map)[, 1])
    y_coords <- c(agCoords(map)[, 2])
    z_coords <- rep(lower_z + 0.02, length(agNames(map))) # was 0.02
    ag_point_size <- c(rep(8, length(agNames(map)))) / 5
    ag_point_size[match(highlighted_ags, agNames(map))] <- 14/5
    ag_col <- c(agOutline(map))
    ag_fill <- c(agFill(map))
    labels <- rep("", length(agNames(map)))
    labels[match(highlighted_ags, agNames(map))] <- ag_plot_names
    
  } else {
    x_coords <- c(agCoords(map)[match(highlighted_ags, agNames(map)), 1])
    y_coords <- c(agCoords(map)[match(highlighted_ags, agNames(map)), 2])
    z_coords <- rep(lower_z + 0.02, length(highlighted_ags)) # was 0.02
    ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
    ag_col <- c(agOutline(map)[match(highlighted_ags, agNames(map))])
    ag_fill <- c(agFill(map)[match(highlighted_ags, agNames(map))])
    labels <- ag_plot_names
  }
  
  border_col <- "grey50" 
  z_lims <- c(lower_z,upper_z)
  axis_at <- seq(z_lims[1]+1, z_lims[2],z_by)
  # Setup plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid =  0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "log"
    ),
    show.axis = FALSE
  )
  
  if(add_axis){
    
    axis_labels <- 2^axis_at*10
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
      # labeloffset = 0.11,
      cornerside = "f",
      size = 20,
      alignment = "right"
    )
  }
  
  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  data3js <- r3js::points3js(
    data3js,
    x          = x_coords,
    y          = y_coords,
    z          = z_coords,
    size       = ag_point_size,
    col        = ag_col,
    fill       = ag_fill,
    lwd        = 0.5,
    opacity    = 1,
    highlight  = list(col = "red"),
    label      = labels,
    toggle     = "Basepoints",
    depthWrite = FALSE,
    shape      = "circle filled"
  )
  
  if(add_ag_label){
    
    x_coords[grep("XBB.1.5", labels)] <- x_coords[grep("XBB.1.5", labels)] - 0.7
    x_coords[grep("BQ.1.1", labels)] <- x_coords[grep("BQ.1.1", labels)] + 0.7
    
    data3js <- r3js::text3js(
      data3js,
      x          = x_coords,
      y          = y_coords - 0.5,
      z          = z_coords,
      text       = labels,
      #  toggle     = "Labels",
      size       = c(rep(14*0.02, length(x_coords))), #agSize(map)[agNames(map) %in% highlighted_ags]*0.02,
      alignment  = "center"
    )
  }
  
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[1]), y = c(lims$ylim[1], lims$ylim[2]), z = rep(lower_z, 2),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[2]), z = rep(lower_z, 2),
                        lwd = 1.2, col = border_col)
    
    # y border
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[1]), z = rep(lower_z, 2),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[2], lims$ylim[2]), z = rep(lower_z, 2),
                        lwd = 1.2, col = border_col)
    
    data3js <- r3js::box3js(
      data3js,
      col   = border_col
    )
    
  }
  
  return(data3js)
}

plot_idvl_landscapes_from_list <- function(data3js, idvl_landscapes, sr_colors){
 
  for (x in seq_along(idvl_landscapes)) {
    
    surface_options <- list()
    surface_options$col.surface = sr_colors[x]
    surface_options$col.surface.grid = adjustcolor(
      "grey",
      red.f = 0.25,
      green.f = 0.25,
      blue.f = 0.25
    )
    surface_options$opacity.surface = 0.2
    
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = idvl_landscapes[[x]],
      toggle = x,
      options = surface_options,
      crop2chull = FALSE,
      grid_spacing = 0.5,
      padding = padding
    )
    
  }
  
  return(data3js)
}


plot_landscapes_from_list <- function(data3js, titertables_groups, lndscp_fits,map, gmt_data, highlighted_ags,
                                      ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1, hide_buttons = TRUE, add_ag_label = FALSE, lndscp_colors,
                                      show_gmts = TRUE){
  
  if(alternative_ba5){
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags)), 12) / 5
    text_x <- c(x_coords[1:6] - ag_point_size[1:6]*0.2)
    text_y <- c(y_coords[1:6] - ag_point_size[1:6]*0.2)
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]], "BA.4/BA.5(2)")
    
  } else {
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
    text_x <- c(x_coords[1:5] - ag_point_size[1:5]*0.2)
    text_y <- c(y_coords[1:5] - ag_point_size[1:5]*0.2)
  }

  if(add_ag_label){
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]])
  } else {
    text_plot <- rep("", length(highlighted_ags))
  }
  
  if(length(lndscp_fits) > 1){
    min_offset <- -0.1
  max_offset <- 0.1
  } else {
    min_offset <- 0
  max_offset <- 0
  }
  offset <- seq(from = min_offset, to = max_offset, by = (max_offset - min_offset)/length(lndscp_fits))
    

  for (i in seq_along(lndscp_fits)) {
    
   # message(i)
    srg <- as.character(titertables_groups$sr_group[i])
    lndscp_fit <- lndscp_fits[[i]]
    
    coords <- cbind(x_coords, y_coords)
    
    coords <- coords[!is.na(x_coords),]
    # Add titers
    gmts <- filter(gmt_data, sr_group == srg)
    
    gmts <- gmts[match(rownames(coords), gmts$ag_name),]
   
    for (j in seq_len(nrow(coords))) {
      
      if(show_gmts){

      data3js <- r3js::lines3js(
        data3js,
        x = rep(coords[j, 1], 2),
        y = rep(coords[j, 2], 2),
        z = c(0, gmts$gmt[j]),
        col = "grey50",
        highlight = list(col = "red"),
        interactive = FALSE,
        toggle = sprintf("GMT, %s", srg),
        geometry = TRUE,
        opacity = 0.7,
        lwd = 0.2 #was 0.4
      )

      
      data3js <- r3js::points3js(
        data3js,
        x         = coords[j, 1],# + offset[i],
        y         = coords[j, 2],
        z         = gmts$gmt[j],
      #  size      = 0.7, #was 2
        size      = 0.9, #was 0.9
      #  col       = "grey50",
        col  = lndscp_colors[srg, 'Color'],
        highlight = list(col = "red"),
   #     label     = gmts$variant[j],
        toggle = sprintf("GMT, %s", srg),
        opacity   = 1 # was 1
     
      )
      }
     # text_x <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 1] + agSize(map)[agNames(map) %in% highlighted_ags[1:4]]*0.025), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 1]- agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
    #  text_y <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 2]), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 2] - agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
      
      # set points and coordinates of highlighted ags

      
      data3js <- r3js::text3js(
        data3js,
        x          = text_x,
        y          = text_y,
        z          = z_coords,
        text       = text_plot,
        size       = c(rep(14*0.02, length(text_x))), #agSize(map)[agNames(map) %in% highlighted_ags]*0.02,
        alignment  = "right"
      )
      
    }
    
    # Add landscapes
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = lndscp_fit,
      # zlim = c(0, 10),
      crop2chull = FALSE,
      # crop2base = TRUE,
      toggle = sprintf("Landscape, %s", srg),
      grid_spacing = 0.5,
      padding = padding,
      options = list(
        col.surface = lndscp_colors[srg, 'Color'],
       # opacity.surface = 0.5
        opacity.surface = opacity_val
      )
    )
    
  }
  
  if(hide_buttons){
    data3js <- remove_buttons(data3js)
  }
 
  
  
  return(data3js)
}

# set orientation
set_r3js_orentation <- function(data3js, angle =  list(
  rotation = c(-1.4681,0.004,-0.0162),
  translation = c(0, 0.05,0.1), 
  zoom = 1.45
)){
  r3js(
    data3js,
    rotation = angle$rotation,
    zoom = angle$zoom
  )
}

# sams landscape functions to add landscape from lndscp fits list
get_titertable <- function(data, group) {
  
  data %>% 
    select(
      ag_name,
      sr_name,
      titer
    ) %>%
    mutate(
      titer = replace(titer, is.na(titer), "*")
    ) %>%
    pivot_wider(
      id_cols = sr_name,
      names_from = ag_name,
      values_from = titer
    ) %>% 
    as.matrix() -> titermatrix
  
  attr(titermatrix, "sr_group") <- group$sr_group
  rownames(titermatrix) <- titermatrix[,"sr_name"]
  titermatrix <- titermatrix[,-1]
  
#  print(titermatrix)
#  print(titermatrix[titermatrix[,"BA.4/BA.5"] != "*",,drop=F])
#  titermatrix[titermatrix[,"BA.4/BA.5"] != "*",,drop=F]
#  titermatrix[titermatrix[,"BA.4/BA.5(2)"] != "*",,drop=F]
  
  return(titermatrix)
  
}


plot_single_landscape_panel_webshot <- function(landscape, label, label_size = 10, label_x_pos = 2, label_y_pos = 9,
                                        sr_group_label = "", sr_group_y_pos = 0, sr_group_size = 3, show_border = FALSE,
                                        delete_html = TRUE, save_name = "temp"){
  
  
  to_save <- file.path(paste0(save_name, ".html"))
  png_save <- gsub(".html", ".png", to_save)
  saveWidget(landscape, to_save, selfcontained = FALSE)
  webshot(url=to_save,file = png_save)
  temp_plot <- readPNG(png_save)
  
  qplot(c(1:10),c(1:10), geom="blank") +
    annotation_custom(rasterGrob(temp_plot, height = unit(0.7, "npc")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    annotate(geom="text", x=label_x_pos, y=label_y_pos, label=label,size= label_size, hjust = 0) + 
    annotate(geom="text", x=label_x_pos, y=sr_group_y_pos, label=sr_group_label,size= sr_group_size, hjust = 0) +
    theme_void() -> p
  
  if(show_border) {
    p + theme(panel.border = element_rect(color = "grey50",
                                          fill = NA,
                                          size = 0.5))-> p
  }
  
  if(delete_html){
    if (file.exists(to_save)) {
      #Delete file if it exists
      file.remove(to_save)
    }
    if (file.exists(png_save)) {
      #Delete file if it exists
      file.remove(png_save)
    }
  }
  
  return(p) 
}


plot_single_landscape_panel <- function(landscape, label, label_size = 10, label_x_pos = 2, label_y_pos = 9,
                                        sr_group_label = "", sr_group_y_pos = 0, sr_group_size = 3, show_border = FALSE,
                                        delete_html = TRUE, save_name = "temp"){
  
  
  to_save <- file.path(paste0(save_name, ".html"))
  png_save <- gsub(".html", ".png", to_save)
  saveWidget(landscape, to_save, selfcontained = FALSE)
  
}


# Residual investigation
residuals_to_long <- function(residuals, values_name = "residuals"){
  
  ag_names <- colnames(residuals)
  as.data.frame(residuals) %>%
    rownames_to_column(var = "sr_name") %>%
    pivot_longer(cols = all_of(ag_names), names_to = "ag_name", values_to = values_name) -> residuals_long
  
  
  return(residuals_long)
  
  
}

combine_residuals <- function(fit, sr_group_fields = 1){
  
  
  residuals <- fit$residuals
  less_thans <- fit$residuals.lessthan
  more_thans <- fit$residuals.morethan
  
  long_res <- residuals_to_long(residuals, "residuals")
  long_less <- residuals_to_long(less_thans, "less_than")
  long_more <- residuals_to_long(more_thans, "more_than")
  
  serum_group <- paste0(str_split(long_res$sr_name[1], "_")[[1]][sr_group_fields], collapse = "_")
  
  comb <- long_res %>%
    left_join(long_less, by = c("ag_name", "sr_name")) %>%
    left_join(long_more, by = c("ag_name", "sr_name")) %>%
    mutate(measured = ifelse(is.na(residuals), ifelse(is.na(less_than), "more_than", "less_than"), "detectable"),
           residuals = ifelse(is.na(residuals), ifelse(is.na(less_than), more_than, less_than), residuals)) %>%
    select(!less_than:more_than) %>%
    mutate(sr_group = serum_group)
  
  return(comb)
}


combine_landscape_and_calculated_gmt <- function(lndscp_fits, gmt_data, sr_group_fields = 1){
  
  lndscp_gmts <- lapply(lndscp_fits, function(x){
    
    gmts <- data.frame(logtiter = x$fitted.value,
                       ag_name = names(x$fitted.value),
                       sr_group = paste0(str_split(rownames(x$logtiters)[1], "_")[[1]][sr_group_fields], collapse = "_"))
    return(gmts)
  })
  
  lndscp_gmts <- do.call(rbind, lndscp_gmts)
  
  
  ## copmare lndscp gmts and claculated gmts
  comb_gmt <- rbind(lndscp_gmts %>%
                      mutate(Data = "Fitted Landscape GMT"),
                    gmt_data %>%
                      mutate(logtiter = gmt) %>%
                      select(sr_group, ag_name, logtiter) %>%
                      unique() %>%
                      mutate(Data = "Calculated GMT"))
  
  
  return(comb_gmt)
}


plot_lndscp_calculated_gmt_lineplot <- function(comb_data, ag_order, lower_lim = -1, plot_colors = NULL){
  
  comb_data %>%
    ggplot(aes(x = ag_name, y = logtiter, color = Data, fill = Data)) + 
    geom_line(aes(group = Data), position = position_dodge(width = 0.3)) + 
    geom_point(shape = 21, color = "grey20", position = position_dodge(width = 0.3)) + 
    # scale_fill_manual(values = plot_colors,
    #                   name = "Serum group") +
    # scale_color_manual(values = plot_colors,
    #                    name = "Serum group") +
    scale_x_discrete(name = "Variant",
                     limits = ag_order) + 
    scale_y_continuous(limits = c(lower_lim, NA),
                       labels = function(x) round(2^x*10,2), 
                       breaks = c(-3:10),
                       name = "GMT") + 
    annotate(
      "rect",
      xmin= -Inf,
      xmax = Inf,
      ymin = -Inf, 
      ymax = log2(16/10),
      #  x = sub_ags,
      #  ymin = 2-1, #min(c(-2, log2(20/10))),
      #  height = 2*(1 + log2(20/10)),
      fill = "grey50",
      color = NA,
      alpha = 0.2
    ) +
    facet_wrap(~sr_group,
               labeller = label_wrap_gen(multi_line = TRUE),
               ncol = 2) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7)) -> p
  
  if(!is.null(plot_colors)){
    p <- p +
      scale_color_manual(values = plot_colors) + 
      scale_fill_manual(values = plot_colors)
  }
  
  return(p)
  
}

rmse_per_variant <- function(lndscp_fits, sr_group_fields = 1){
  all_residuals <- lapply(lndscp_fits, function(x) combine_residuals(x, sr_group_fields))
  
  all_residuals <- do.call(rbind, all_residuals)
  
  all_residuals %>%
    filter(!is.na(residuals)) %>%
    group_by(sr_group) %>%
    summarize(ag_name = "Total", 
              rmse = sqrt(sum(residuals^2, na.rm = TRUE)/(length(residuals))), # Residual standard error
              residual_type = "All variants") -> total_ag
  
  all_residuals %>%
    filter(!is.na(residuals)) %>%
    group_by(ag_name, sr_group) %>%
    mutate(rmse = sqrt(sum(residuals^2, na.rm = TRUE)/(length(residuals))),
           residual_type = "By variant") %>%
    plyr::rbind.fill(., total_ag) -> ssr
  
  return(ssr)
}


plot_total_residuals_by_base_map <- function(sub_residuals, fit_order, ag_colors, target_group,
                                             ymax = 3){
  
  sub_residuals <- sub_residuals %>%
    filter(sr_group == target_group)
  
  fit_order <- fit_order %>%
    filter(sr_group == target_group) %>%
    pull(Data)
  
  sub_residuals %>%
    ggplot(aes(x = Data, y = rmse, color = ag_name, fill = ag_name, group = ag_name)) + 
    # geom_line(position = position_dodge(width = 0.2)) +
    geom_point(color = "grey20", shape = 21, position = position_dodge(width = 0.2)) + 
    scale_x_discrete(limits = fit_order,
                     name = "Base Map") +
    ylab("RMSE") + 
    ylim(c(0,ymax)) +
    facet_wrap(~sr_group) + 
    scale_fill_manual(values = ag_colors, name = "Variant") +
    scale_color_manual(values = ag_colors, name = "Variant") +
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) -> residuals_plot_base_map
  
  return(residuals_plot_base_map)
}

plot_gmt_diff_per_base_map <- function(gmt_diff, fit_order, ag_colors, target_group,
                                       ymin = -2, ymax = 2){
  
  gmt_diff <- gmt_diff %>%
    filter(Data != "Calculated GMT") %>%
    filter(sr_group == target_group)
  
  fit_order <- fit_order %>%
    filter(sr_group == target_group) %>%
    pull(Data)
  
  gmt_diff %>%
    ggplot(aes(x = Data, y = gmt_diff, color = ag_name, fill = ag_name, group = ag_name)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(color = "grey20", shape = 21, position = position_dodge(width = 0.2)) + 
    scale_x_discrete(limits = fit_order,
                     name = "Base Map") +
    scale_y_continuous(name = "GMT difference",
                       limits = c(ymin, ymax),
                       breaks = seq(ymin, ymax, 0.5),
                       labels = function(x) paste0(round(2^x, 1), "x")) +
    facet_wrap(~sr_group) + 
    scale_fill_manual(values = ag_colors, name = "Variant") +
    scale_color_manual(values = ag_colors, name = "Variant") +
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1)) -> p_gmt_diff
  
  return(p_gmt_diff)
  
}

landscape_lapply_fit <- function(titertables, map, ags_to_fit_lndscp, 
                                 method = "cone",
                                 error_sd = 1,
                                 bandwidth = 1,
                                 degree = 1,
                                 control_list = list(optimise.cone.slope = TRUE,
                                                     optimise.cone.coords = TRUE)){
  
  lndscp_fits <- lapply(
    titertables,
    function(titertable) {
      if(is.null(dim(titertable))){
        t_table <- matrix(titertable, nrow = 1, ncol = length(titertable))
        colnames(t_table) <- names(titertable)
        
      } else {
        t_table <- titertable
      }
      
      ablandscape.fit(
        titers = t_table[,ags_to_fit_lndscp],
        bandwidth = bandwidth,
        degree = degree,
        method = method,
        error.sd = error_sd,
        acmap = map,
        control = control_list
      )
      
      
    }
  )
  
  return(lndscp_fits)
  
}


fit_cone_landscapes_stan_single_sr_group <- function(model, titerdata_b, target_sr_group, map, n_cones = 2, 
                                                     init_cone_centers = matrix(rep(0, n_cones*2), n_cones, 2),
                                                     init_cone_heights = rep(6, n_cones),
                                                     init_cone_slopes = rep(0.8, n_cones),
                                                     sum_cone_eval = 0, 
                                                     lower_lod = 16,
                                                     upper_lod = 16*2^10,
                                                     std_dev_titer = 0.25,
                                                     std_dev_slope = 0.5,
                                                     std_dev_height = 3,
                                                     fit_adjusted_titer = TRUE,
                                                     iterations_sampling = 1000,
                                                     set_adapt_delta = 0.95,
                                                     upper_cone_slope = 3,
                                                     max_cone_height = 12,
                                                     std_dev_cone_centers = 1){
  
  titerdata_b %>%
    filter(sr_group == target_sr_group) -> titerdata
  
  if(fit_adjusted_titer){
    # let's get it in titertable form
    titerdata %>%
      select(ag_name, sr_name, logtiter_adjusted) %>%
      pivot_wider(names_from = "sr_name", values_from = "logtiter_adjusted") %>%
      column_to_rownames("ag_name")-> titer_table_b
  } else {
    
    titerdata %>%
      select(ag_name, sr_name, logtiter_adjusted) %>%
      pivot_wider(names_from = "sr_name", values_from = "logtiter") %>%
      column_to_rownames("ag_name")-> titer_table_b
    
  }
  
  
  n_ag <- nrow(titer_table_b)
  n_sera <- ncol(titer_table_b)
  titer_table <- titer_table_b[1:n_ag, 1:n_sera]
  n_data <- n_ag*n_sera
  # get initial position of cone center, option to provide coordinates but here I'll set it to 0 for now
  # assumption here is all sera have the same apex position. might need to change that. 
  # same height is OK because we fit the bias
  ag_coords <- agCoords(map)[1:n_ag,]
  lower_lod <- log2(lower_lod/10)
  upper_lod <- log2(upper_lod/10)
  apex_min <- min(ag_coords)-0.5
  apex_max <- max(ag_coords)+0.5
  observed_data <- as.matrix(titer_table)
  

  if(n_cones ==1){
    

    stan_init <- list(
      apex_x_coords = array(init_cone_centers[,1], dim = 1),
      apex_y_coords = array(init_cone_centers[,2], dim = 1),
      cone_slopes = array(rep(init_cone_slopes[1], n_cones), dim = 1), # if it is a single cone, it needs to be passed as array(value, dim = 1)like this otherwise interpreted by stan as scalar
      cone_heights =array(rep(init_cone_heights[1], n_cones), dim = 1),
      bias = rep(0, n_sera)
    )
    
  } else {
    # initial values for parameters
    stan_init <- list(
      apex_x_coords = init_cone_centers[,1],
      apex_y_coords = init_cone_centers[,2],
      cone_slopes = init_cone_slopes, 
      cone_heights = init_cone_heights,
      bias = rep(0, n_sera)
    )
  }
  
  x_range <- range(ag_coords[,1])
  x_range[2] <- x_range[2]-1
  
  y_range <- range(ag_coords[,2])
  
  # Set input data
  stan_data <- list(
    Ndata =nrow(titerdata),
    ncones = n_cones,
    nag = nrow(titer_table),
    nsera = ncol(titer_table),
    init_centers = init_cone_centers,
    init_heights = init_cone_heights,
    init_slopes = init_cone_slopes,
    ag_coordinates = ag_coords,
    apex_x_range = x_range,
    apex_y_range = y_range,
    lower_lod = lower_lod,
    upper_lod = upper_lod,
    apex_min = apex_min,
    apex_max = apex_max,
    cone_slope_max = upper_cone_slope,
    cone_height_max = max_cone_height,
    std_dev = std_dev_titer,
    std_dev_slope = std_dev_slope,
    std_dev_height = std_dev_height,
    std_dev_apex_centers = std_dev_cone_centers,
    observed = observed_data,
    sum_cone_eval = sum_cone_eval # FALSE FOR max eval
  )
  
  # let's try to fit it
  fit <- model$sample(
    data = stan_data,
    seed = 123,
    init = list(
      stan_init,
      stan_init,
      stan_init,
      stan_init
    ),
    save_warmup = FALSE,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = iterations_sampling, # was 10000
    adapt_delta = set_adapt_delta,
    max_treedepth = 20,
    refresh = 500 # print update every 500 iters,
  )
  
  
  # this is the fit above, need to return it for diagnostics
  
  # but also need to assign sr name, etc., to bias, residuals, etc
  cone_vals <- list()
  cone_vals[[target_sr_group]] <- list("cone_coords" = matrix(fit$summary("apex_centres")$mean, ncol = 2, byrow = FALSE),
                                       "cone_heights" = fit$summary("cone_heights")$mean,
                                       "cone_slopes" = fit$summary("cone_slopes")$mean,
                                       "combined_cone_fit" = fit$summary(variables = "combined_cone_y")$mean,
                                       "sum_cone_eval" = sum_cone_eval
  )
  
  # organize results test
  fit$summary("bias") %>%
    mutate(sr_name = colnames(titer_table_b),
           variable = "sr_bias",
           sr_bias = mean) %>%
    select(sr_name, variable, sr_bias) -> sr_bias_df
  
  fit$summary("residuals") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    separate(
      indices,
      sep = ",",
      into = c("ag_ind", "sr_ind")
    ) %>%
    mutate(ag_name = rownames(titer_table_b)[as.numeric(ag_ind)],
           sr_name = colnames(titer_table_b)[as.numeric(sr_ind)],
           residuals = mean) %>%
    select(ag_name, sr_name, residuals) -> residuals_df
  
  fit$summary("combined_cone_y") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    mutate(ag_name = rownames(titer_table_b)[as.numeric(indices)],
           mean_logtiter_fitted = mean) %>%
    select(ag_name, mean_logtiter_fitted) -> fitted_logtiter_df
  
  titerdata %>%
    left_join(sr_bias_df, by = "sr_name") %>%
    left_join(residuals_df, by = c("sr_name", "ag_name")) %>%
    left_join(fitted_logtiter_df, by = "ag_name") %>%
    mutate(cone_nr = n_cones,
           log_likelihood = fit$summary("lp__")$mean,
           sum_cone_eval = sum_cone_eval) -> titerdata
  
  return_list <- list("fit" = fit,
                      "cone_vals" = cone_vals,
                      "titerdata" = titerdata)
  return(return_list)
}


fit_cone_landscapes_stan_single_sr_group_gmt <- function(model, titerdata_b, target_sr_group, map, n_cones = 2, 
                                                     init_cone_centers = matrix(rep(0, n_cones*2), n_cones, 2),
                                                     init_cone_heights = rep(6, n_cones),
                                                     init_cone_slopes = rep(0.8, n_cones),
                                                     sum_cone_eval = 0, 
                                                     lower_lod = 16,
                                                     upper_lod = 16*2^10,
                                                     std_dev_titer = 0.25,
                                                     std_dev_slope = 0.5,
                                                     std_dev_height = 3,
                                                     fit_adjusted_titer = TRUE,
                                                     iterations_sampling = 1000,
                                                     set_adapt_delta = 0.95,
                                                     upper_cone_slope = 3,
                                                     max_cone_height = 12,
                                                     std_dev_cone_centers = 1,
                                                     gmt_data = NULL){
  
  titerdata_b %>%
    filter(sr_group == target_sr_group) -> titerdata
  
  if(fit_adjusted_titer){
    # let's get it in titertable form
    titerdata %>%
      select(ag_name, sr_name, logtiter_adjusted) %>%
      pivot_wider(names_from = "sr_name", values_from = "logtiter_adjusted") %>%
      column_to_rownames("ag_name")-> titer_table_b
    
    if(is.null(gmt_data)){
      titerdata %>%
        select(ag_name, logtiter_adjusted) %>%
        group_by(ag_name) %>%
        summarize(logtiter_adjusted = mean(logtiter_adjusted)) %>%
        column_to_rownames("ag_name")  -> gmt_data
    }
    
    gmt_data <- gmt_data[rownames(titer_table_b),1]
    
  } else {
    titerdata %>%
      select(ag_name, sr_name, logtiter) %>%
      pivot_wider(names_from = "sr_name", values_from = "logtiter") %>%
      column_to_rownames("ag_name")-> titer_table_b
    
    if(is.null(gmt_data)){
      titerdata %>%
        select(ag_name, logtiter) %>%
        group_by(ag_name) %>%
        summarize(logtiter_adjusted = mean(logtiter)) %>%
        column_to_rownames("ag_name")  -> gmt_data
    }
    
    gmt_data <- gmt_data[rownames(titer_table_b),1]
    
    
  }
  
  
  n_ag <- nrow(titer_table_b)
  n_sera <- ncol(titer_table_b)
  titer_table <- titer_table_b[1:n_ag, 1:n_sera]
  n_data <- n_ag*n_sera
  # get initial position of cone center, option to provide coordinates but here I'll set it to 0 for now
  # assumption here is all sera have the same apex position. might need to change that. 
  # same height is OK because we fit the bias
  ag_coords <- agCoords(map)[1:n_ag,]
  lower_lod <- log2(lower_lod/10)
  upper_lod <- log2(upper_lod/10)
  apex_min <- min(ag_coords)-0.5
  apex_max <- max(ag_coords)+0.5
  observed_data <- as.matrix(titer_table)
  
  
  if(n_cones ==1){
    
    
    stan_init <- list(
      apex_x_coords = array(init_cone_centers[,1], dim = 1),
      apex_y_coords = array(init_cone_centers[,2], dim = 1),
      cone_slopes = array(rep(init_cone_slopes[1], n_cones), dim = 1), # if it is a single cone, it needs to be passed as array(value, dim = 1)like this otherwise interpreted by stan as scalar
      cone_heights =array(rep(init_cone_heights[1], n_cones), dim = 1),
      bias = rep(0, n_sera)
    )
    
  } else {
    # initial values for parameters
    stan_init <- list(
      apex_x_coords = init_cone_centers[,1],
      apex_y_coords = init_cone_centers[,2],
      cone_slopes = init_cone_slopes, 
      cone_heights = init_cone_heights,
      bias = rep(0, n_sera)
    )
  }
  
  x_range <- range(ag_coords[,1])
 # x_range[2] <- x_range[2]-1
  
  y_range <- range(ag_coords[,2])
  
  # Set input data
  stan_data <- list(
    Ndata =nrow(titerdata),
    ncones = n_cones,
    nag = nrow(titer_table),
    nsera = ncol(titer_table),
    init_centers = init_cone_centers,
    init_heights = init_cone_heights,
    init_slopes = init_cone_slopes,
    ag_coordinates = ag_coords,
    apex_x_range = x_range,
    apex_y_range = y_range,
    lower_lod = lower_lod,
    upper_lod = upper_lod,
    apex_min = apex_min,
    apex_max = apex_max,
    cone_slope_max = upper_cone_slope,
    cone_height_max = max_cone_height,
    std_dev = std_dev_titer,
    std_dev_slope = std_dev_slope,
    std_dev_height = std_dev_height,
    std_dev_apex_centers = std_dev_cone_centers,
    observed = observed_data,
    observed_gmt = gmt_data,
    sum_cone_eval = sum_cone_eval # FALSE FOR max eval
  )
  
  # let's try to fit it
  fit <- model$sample(
    data = stan_data,
    seed = 123,
    init = list(
      stan_init,
      stan_init,
      stan_init,
      stan_init
    ),
    save_warmup = FALSE,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = iterations_sampling, # was 10000
    adapt_delta = set_adapt_delta,
    max_treedepth = 20,
    refresh = 500 # print update every 500 iters,
  )
  
  
  # this is the fit above, need to return it for diagnostics
  
  # but also need to assign sr name, etc., to bias, residuals, etc
  cone_vals <- list()
  cone_vals[[target_sr_group]] <- list("cone_coords" = matrix(fit$summary("apex_centres")$mean, ncol = 2, byrow = FALSE),
                                       "cone_heights" = fit$summary("cone_heights")$mean,
                                       "cone_slopes" = fit$summary("cone_slopes")$mean,
                                       "combined_cone_fit" = fit$summary(variables = "combined_cone_y")$mean,
                                       "sum_cone_eval" = sum_cone_eval
  )
  
  # organize results test
  fit$summary("bias") %>%
    mutate(sr_name = colnames(titer_table_b),
           variable = "sr_bias",
           sr_bias = mean) %>%
    select(sr_name, variable, sr_bias) -> sr_bias_df
  
  fit$summary("residuals") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    separate(
      indices,
      sep = ",",
      into = c("ag_ind", "sr_ind")
    ) %>%
    mutate(ag_name = rownames(titer_table_b)[as.numeric(ag_ind)],
           sr_name = colnames(titer_table_b)[as.numeric(sr_ind)],
           residuals = mean) %>%
    select(ag_name, sr_name, residuals) -> residuals_df
  
  fit$summary("combined_cone_y") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    mutate(ag_name = rownames(titer_table_b)[as.numeric(indices)],
           mean_logtiter_fitted = mean) %>%
    select(ag_name, mean_logtiter_fitted) -> fitted_logtiter_df
  
  titerdata %>%
    left_join(sr_bias_df, by = "sr_name") %>%
    left_join(residuals_df, by = c("sr_name", "ag_name")) %>%
    left_join(fitted_logtiter_df, by = "ag_name") %>%
    mutate(cone_nr = n_cones,
           log_likelihood = fit$summary("lp__")$mean,
           sum_cone_eval = sum_cone_eval) -> titerdata
  
  return_list <- list("fit" = fit,
                      "cone_vals" = cone_vals,
                      "titerdata" = titerdata)
  return(return_list)
}


norm_data_to_gmt <- function(titerdata, gmt_data){
  
  gmt_to_join <- gmt_data %>%
    mutate(gmt = logtiter) %>%
    select(gmt, ag_name, sr_group)

  titerdata %>%
    left_join(.,gmt_to_join, by = c("ag_name", "sr_group")) -> titerdata_gmt
  
  titerdata_gmt %>%
    group_by(sr_name) %>%
    mutate(mean_sr_titer = mean(logtiter)) %>%
    ungroup() %>%
    group_by(sr_group) %>%
    mutate(mean_gmt_group = mean(gmt)) %>%
    ungroup() %>%
    mutate(sr_norm = mean_gmt_group - mean_sr_titer,
           raw_logtiter = logtiter,
           logtiter = logtiter + sr_norm) -> titerdata_gmt
  
  return(titerdata_gmt)
  
}


fit_cone_landscapes_stan_single_sr_group_gmt_studentT <- function(model, titerdata_b, target_sr_group, map, n_cones = 2, 
                                                         init_cone_centers = matrix(rep(0, n_cones*2), n_cones, 2),
                                                         init_cone_heights = rep(6, n_cones),
                                                         init_cone_slopes = rep(0.8, n_cones),
                                                         sum_cone_eval = 0, 
                                                         lower_lod = 16,
                                                         upper_lod = 16*2^10,
                                                         std_dev_cone_centers = 0.5,
                                                         init_val_noise = 3,
                                                         std_dev_titer = 0.25,
                                                         std_dev_slope = 0.25,
                                                         std_dev_height = 0.5,
                                                         init_val_nu = 10,
                                                         init_sd_nu = 3,
                                                         sr_bias_sd = 3,
                                                         fit_adjusted_titer = TRUE,
                                                         iterations_sampling = 1000,
                                                         set_adapt_delta = 0.95,
                                                         upper_cone_slope = 3,
                                                         max_cone_height = 12,
                                                         gmt_data = NULL,
                                                         normalize_data = FALSE,
                                                         geometric_sum_eval = 0){
  
  if(normalize_data){
    sr_bias_sd <- 1
  }
  
  titerdata_b %>%
    filter(sr_group == target_sr_group) -> titerdata

  
  if(fit_adjusted_titer){
    # let's get it in titertable form
    titerdata %>%
      mutate(logtiter = logtiter_adjusted) -> titerdata
    
  }

  if(is.null(gmt_data)){
    titerdata %>%
      select(ag_name, logtiter, sr_group) %>%
      group_by(ag_name) %>%
      summarize(logtiter = mean(logtiter)) -> gmt_data
  } else {
    gmt_data %>%
      filter(sr_group == target_sr_group) -> gmt_data
  }
  
  
  if(normalize_data){
    titerdata <- norm_data_to_gmt(titerdata, gmt_data)
  }
  
  titerdata %>%
    select(ag_name, sr_name, logtiter) %>%
    pivot_wider(names_from = "sr_name", values_from = "logtiter") %>%
    column_to_rownames("ag_name")-> titer_table_b
  
  gmt_data <- gmt_data %>%
    select(ag_name, logtiter) %>%
    unique() %>%
    column_to_rownames("ag_name")
  gmt_data <- gmt_data[rownames(titer_table_b),1]
  
  
  n_ag <- nrow(titer_table_b)
  n_sera <- ncol(titer_table_b)
  titer_table <- titer_table_b[1:n_ag, 1:n_sera]
  n_data <- n_ag*n_sera
  # get initial position of cone center, option to provide coordinates but here I'll set it to 0 for now
  # assumption here is all sera have the same apex position. might need to change that. 
  # same height is OK because we fit the bias
  ag_coords <- agCoords(map)[1:n_ag,]
  lower_lod <- log2(lower_lod/10)
  upper_lod <- log2(upper_lod/10)
  apex_min <- min(ag_coords)-0.5
  apex_max <- max(ag_coords)+0.5
  observed_data <- as.matrix(titer_table)
  
  
  if(n_cones ==1){
    
    
    stan_init <- list(
      apex_x_coords = array(init_cone_centers[,1], dim = 1),
      apex_y_coords = array(init_cone_centers[,2], dim = 1),
      cone_slopes = array(rep(init_cone_slopes[1], n_cones), dim = 1), # if it is a single cone, it needs to be passed as array(value, dim = 1)like this otherwise interpreted by stan as scalar
      cone_heights =array(rep(init_cone_heights[1], n_cones), dim = 1),
      bias = rep(0, n_sera),
      nu = init_val_nu,
      std_dev = std_dev_titer
    )
    
  } else {
    # initial values for parameters
    stan_init <- list(
      apex_x_coords = init_cone_centers[,1],
      apex_y_coords = init_cone_centers[,2],
      cone_slopes = init_cone_slopes, 
      cone_heights = init_cone_heights,
      bias = rep(0, n_sera),
      nu = init_val_nu,
      std_dev = std_dev_titer
    )
  }
  
  x_range <- range(ag_coords[,1])
  # x_range[2] <- x_range[2]-1
  
  y_range <- range(ag_coords[,2])

  
  # Set input data
  stan_data <- list(
    Ndata =nrow(titerdata),
    ncones = n_cones,
    nag = nrow(titer_table),
    nsera = ncol(titer_table),
    init_centers = init_cone_centers,
    init_heights = init_cone_heights,
    init_slopes = init_cone_slopes,
    ag_coordinates = ag_coords,
    apex_x_range = x_range,
    apex_y_range = y_range,
    lower_lod = lower_lod,
    upper_lod = upper_lod,
    apex_min = apex_min,
    apex_max = apex_max,
    cone_slope_max = upper_cone_slope,
    cone_height_max = max_cone_height,
    init_titer_mu = init_val_noise,
    init_titer_std_dev = std_dev_titer,
    init_nu = init_val_nu,
    init_nu_std_dev = init_sd_nu,
    sr_bias_std_dev = sr_bias_sd,
    std_dev_slope = std_dev_slope,
    std_dev_height = std_dev_height,
    std_dev_apex_centers = std_dev_cone_centers,
    observed = observed_data,
    observed_gmt = gmt_data,
    sum_cone_eval = sum_cone_eval, # FALSE FOR max eval
    geometric_sum_eval = geometric_sum_eval
  )
  
  # let's try to fit it
  fit <- model$sample(
    data = stan_data,
    seed = 123,
    init = list(
      stan_init,
      stan_init,
      stan_init,
      stan_init
    ),
    save_warmup = FALSE,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = iterations_sampling, # was 10000
    adapt_delta = set_adapt_delta,
    max_treedepth = 20,
    refresh = 500 # print update every 500 iters,
  )
  
  
  # this is the fit above, need to return it for diagnostics
  
  # but also need to assign sr name, etc., to bias, residuals, etc
  cone_vals <- list()
  cone_vals[[target_sr_group]] <- list("cone_coords" = matrix(fit$summary("apex_centres")$mean, ncol = 2, byrow = FALSE),
                                       "cone_heights" = fit$summary("cone_heights")$mean,
                                       "cone_slopes" = fit$summary("cone_slopes")$mean,
                                       "combined_cone_fit" = fit$summary(variables = "combined_cone_y")$mean,
                                       "sum_cone_eval" = sum_cone_eval,
                                       "geometric_sum_eval" = geometric_sum_eval
  )
  
  # organize results test
  fit$summary("bias") %>%
    mutate(sr_name = colnames(titer_table_b),
           variable = "sr_bias",
           sr_bias = mean) %>%
    select(sr_name, sr_bias) -> sr_bias_df
  
  fit$summary("residuals") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    separate(
      indices,
      sep = ",",
      into = c("ag_ind", "sr_ind")
    ) %>%
    mutate(ag_name = rownames(titer_table_b)[as.numeric(ag_ind)],
           sr_name = colnames(titer_table_b)[as.numeric(sr_ind)],
           residuals = mean) %>%
    select(ag_name, sr_name, residuals) -> residuals_df
  
  fit$summary("combined_cone_y") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    mutate(ag_name = rownames(titer_table_b)[as.numeric(indices)],
           mean_logtiter_fitted = mean) %>%
    select(ag_name, mean_logtiter_fitted) -> fitted_logtiter_df
  
  fit$summary(c("nu", "std_dev")) %>%
    select(variable, mean) %>%
    pivot_wider(names_from = "variable", values_from = "mean")-> nu_sd_df
  
  titerdata %>%
    left_join(sr_bias_df, by = "sr_name") %>%
    left_join(residuals_df, by = c("sr_name", "ag_name")) %>%
    left_join(fitted_logtiter_df, by = "ag_name") %>%
    mutate(cone_nr = n_cones,
           log_likelihood = fit$summary("lp__")$mean,
           sum_cone_eval = sum_cone_eval,
           geometric_sum_eval = geometric_sum_eval,
           nu = nu_sd_df$nu,
           titer_noise = nu_sd_df$std_dev) -> titerdata
  
  if(normalize_data){
    
    titerdata %>%
      mutate(sr_bias_normalized = sr_bias,
             sr_bias = sr_bias - sr_norm) -> titerdata
    
  }
  
  return_list <- list("fit" = fit,
                      "cone_vals" = cone_vals,
                      "titerdata" = titerdata)
  return(return_list)
}



fit_cone_landscapes_stan_single_sr_group_serum_apex_coords <- function(model, titerdata_b, target_sr_group, map, n_cones = 2, 
                                                     init_cone_slopes = rep(0.6, n_cones),
                                                     sum_cone_eval = 0, 
                                                     lower_lod = 16,
                                                     upper_lod = 16*2^10,
                                                     std_dev_cone = 0.25,
                                                     fit_adjusted_titer = TRUE){
  
 
  titerdata_b %>%
    filter(sr_group == target_sr_group) -> titerdata
  
  if(fit_adjusted_titer){
    # let's get it in titertable form
    titerdata %>%
      select(ag_name, sr_name, logtiter_adjusted) %>%
      pivot_wider(names_from = "sr_name", values_from = "logtiter_adjusted") %>%
      column_to_rownames("ag_name")-> titer_table_b
  } else {
    
    titerdata %>%
      select(ag_name, sr_name, logtiter_adjusted) %>%
      pivot_wider(names_from = "sr_name", values_from = "logtiter") %>%
      column_to_rownames("ag_name")-> titer_table_b
    
  }
  
  n_ag <- nrow(titer_table_b)
  n_sera <- ncol(titer_table_b)
  titer_table <- titer_table_b[1:n_ag, 1:n_sera]
  n_data <- n_ag*n_sera
  ag_coords <- agCoords(map)[1:n_ag,]
  lower_lod <- log2(lower_lod/10)
  upper_lod <- log2(upper_lod/10)
  st_d <- std_dev_cone
  observed_data <- as.matrix(titer_table)
  
  # init data
  init_cone_x <- matrix(0, nrow = n_cones, ncol = n_sera)
  init_cone_y <- matrix(0, nrow = n_cones, ncol = n_sera)
  init_cone_heights <- matrix(6, nrow = n_cones, ncol = n_sera)
  
  stan_init <- list(
    apex_x_coords = init_cone_x,
    apex_y_coords = init_cone_y,
    cone_slopes = init_cone_slopes,
    cone_heights = init_cone_heights
  )
  
  if(n_cones ==1){
    
    init_cone_x <- matrix(ag_coords["D614G",1], nrow = n_cones, ncol = n_sera)
    init_cone_y <- matrix(ag_coords["D614G",2], nrow = n_cones, ncol = n_sera)
    stan_init <- list(
      apex_x_coords = init_cone_x,
      apex_y_coords = init_cone_y,
      cone_slopes = array(rep(init_cone_slopes[1], n_cones), dim = 1), # if it is a single cone, it needs to be passed as array(value, dim = 1)like this otherwise interpreted by stan as scalar
      cone_heights = init_cone_heights
    )
    
  } 
  
  # Set input data
  stan_data <- list(
    Ndata =nrow(titerdata),
    ncones = n_cones,
    nag = nrow(titer_table),
    nsera = ncol(titer_table),
    init_x_centers = init_cone_x,
    init_y_centers = init_cone_y,
    init_heights = init_cone_heights,
    ag_coordinates = ag_coords,
    lower_lod = lower_lod,
    upper_lod = upper_lod,
    std_dev = st_d,
    observed = observed_data,
    sum_cone_eval = sum_cone_eval # FALSE FOR max eval
  )
  
  # let's try to fit it
  fit <- model$sample(
    data = stan_data,
    seed = 123,
    init = list(
      stan_init,
      stan_init,
      stan_init,
      stan_init
    ),
    save_warmup = FALSE,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 1000, # was 10000
    max_treedepth = 20,
    refresh = 500 # print update every 500 iters,
  )
  
  # this is the fit above, need to return it for diagnostics
  
  # but also need to assign sr name, etc., to bias, residuals, etc
  cone_vals_curr <- list("cone_x_coords" = matrix(fit$summary("apex_x_coords")$mean, nrow = n_cones, byrow = FALSE),
                    "cone_y_coords" = matrix(fit$summary("apex_y_coords")$mean, nrow = n_cones, byrow = FALSE),
                    "cone_heights_per_serum" = matrix(fit$summary("cone_heights")$mean, nrow = n_cones, byrow = FALSE),
                    "cone_slopes" = fit$summary("cone_slopes")$mean,
                    "combined_cone_fit_per_serum" = matrix(fit$summary("combined_cone_y_per_serum")$mean, nrow = n_ag, byrow = FALSE),
                    "sum_cone_eval" = sum_cone_eval
                    )
  
  # add mean
  cone_vals_curr$combined_cone_fit <- rowMeans(cone_vals_curr$combined_cone_fit_per_serum)
  cone_vals_curr$cone_coords <- matrix(c(rowMeans(cone_vals_curr$cone_x_coords), rowMeans(cone_vals_curr$cone_y_coords)), nrow = n_cones, ncol = 2)
  cone_vals_curr$cone_heights <- rowMeans(cone_vals_curr$cone_heights_per_serum)
  
  
  cone_vals <- list()
  cone_vals[[target_sr_group]] <- cone_vals_curr
  
  # organize results test
  fit$summary("residuals") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    separate(
      indices,
      sep = ",",
      into = c("ag_ind", "sr_ind")
    ) %>%
    mutate(ag_name = rownames(titer_table_b)[as.numeric(ag_ind)],
           sr_name = colnames(titer_table_b)[as.numeric(sr_ind)],
           residuals = mean) %>%
    select(ag_name, sr_name, residuals) -> residuals_df
  
  # add titer per serum to df
  fit$summary("combined_cone_y_per_serum") %>%
    mutate(
      indices = gsub("^.*\\[", "[", variable),
      indices = gsub("\\[", "", indices),
      indices = gsub("\\]", "", indices)
    ) %>% 
    separate(
      indices,
      sep = ",",
      into = c("ag_ind", "sr_ind")
    ) %>%
    mutate(ag_name = rownames(titer_table_b)[as.numeric(ag_ind)],
           sr_name = colnames(titer_table_b)[as.numeric(sr_ind)],
           logtiter_fitted = mean) %>%
    select(ag_name, sr_name, logtiter_fitted) -> fitted_logtiter_df
  
  
  titerdata %>%
    left_join(residuals_df, by = c("sr_name", "ag_name")) %>%
    left_join(fitted_logtiter_df, by =  c("sr_name", "ag_name")) %>%
    mutate(cone_nr = n_cones,
           log_likelihood = fit$summary("lp__")$mean,
           sum_cone_eval = sum_cone_eval) -> titerdata
  
  return_list <- list("fit" = fit,
                      "cone_vals" = cone_vals,
                      "titerdata" = titerdata)
  
  
  return(return_list)
  
}
  




fit_lndscp_val <- function(x, y, cone_coords, cone_heights, cone_slopes, sr_bias = 0, sum_cone_eval = 0, geometric_sum_eval = 0) {
  
  vals <- lapply(1:nrow(cone_coords), function(cone) {
    norm_dists <- as.matrix(dist(rbind(c(x, y), cone_coords[cone,])))[1,-1]
    
    val <- (cone_heights[cone] - norm_dists*cone_slopes[cone]) + sr_bias
  })
  
  
  
  if(sum_cone_eval){
    if(geometric_sum_eval){
      return(sum(unlist(vals)))
    } else {
      return(log2(sum(2^unlist(vals))))
    }
  } else {
    return(max(unlist(vals)))
  }
  
  
}

calculate_gmt_landscape_from_idvl_landscapes <- function(grid_z_matrices){
  
  grid_z_matrices <- lapply(grid_z_matrices, function(mat){
    
    as.data.frame(mat) %>%
      rownames_to_column(var = "x") %>%
      pivot_longer(cols = colnames(.)[colnames(.) != "x"], names_to = "y", values_to = "value")
    
  })
  
  grid_z_comb <- do.call(rbind, grid_z_matrices)
  grid_z_comb %>%
    group_by(x,  y) %>%
    summarize(value = mean(value)) %>%
    mutate(x = as.numeric(x),
           y = as.numeric(gsub("V", "", y))) -> grid_z_comb
  
  grid_z_comb %>%
    arrange(x, y) %>%
    pivot_wider(names_from = "y", values_from = "value") %>%
    ungroup() %>%
    select(!x) -> mean_grid_z
  
  return(as.matrix(mean_grid_z))
  
}


plot_cone_landscapes <- function(data3js, map, cone_vals, 
                                 sr_groups, # serum group for landscape color
                                 lndscp_colors, # colors for lndscp
                                 sum_cone_eval = 0, # maximum or sum cones
                                 sr_bias = 0, # bias for individual samples
                                 show_cone_coordinates = TRUE,
                                 show_cone_heights = TRUE,
                                 show_fitted_values = FALSE,
                                 show_measured_values = TRUE,
                                 gmt_data = NULL,
                                 show_idvl_landscapes = FALSE,
                                 calculate_gmt_from_idvl_landscapes = FALSE
                                 ){
  

  lims <- Racmacs:::mapPlotLims(map, sera = FALSE)
  lndscp_xlim <- lims$xlim
  lndscp_ylim <- lims$ylim
  
  grid_x_coords <- seq(from = lndscp_xlim[1], to = lndscp_xlim[2], by = 0.25)
  grid_y_coords <- seq(from = lndscp_ylim[1], to = lndscp_ylim[2], by = 0.25)
  grid_x_matrix <- matrix(grid_x_coords, length(grid_y_coords), length(grid_x_coords), byrow = T)
  grid_y_matrix <- matrix(grid_y_coords, length(grid_y_coords), length(grid_x_coords), byrow = F)
  grid_z_matrix <- matrix(NA, length(grid_y_coords), length(grid_x_coords))
  
  ag_coords <- agCoords(map)
  
  for(srg in sr_groups){
    
    temp_cone <- cone_vals[[srg]]
    cone_slopes <- temp_cone$cone_slopes
    combined_cone_fit <- temp_cone$combined_cone_fit
    
    if(show_idvl_landscapes | calculate_gmt_from_idvl_landscapes){
      
      if(is.null(temp_cone$cone_x_coords)){
        warning("You want to show individual landscapes but did not provide the data.")
      }
      
      
      grid_z_matrices <- list()
      
      for(sr in 1:ncol(temp_cone$cone_x_coords)) {
        
      
        cone_coords <- matrix(c(temp_cone$cone_x_coords[,sr], temp_cone$cone_y_coords[,sr]), nrow = length(temp_cone$cone_slopes), byrow = FALSE)
        cone_heights <- temp_cone$cone_heights_per_serum[,sr]
        
  
        # calculate landscape values
        grid_z_matrix[] <- vapply(
          seq_len(length(grid_z_matrix)),
          \(n) {
            fit_lndscp_val(
              grid_x_matrix[n], grid_y_matrix[n],
              cone_coords = cone_coords, 
              cone_heights = cone_heights,
              cone_slopes = cone_slopes, 
              sr_bias = sr_bias, 
              sum_cone_eval = temp_cone$sum_cone_eval,
              geometric_sum_eval = temp_cone$geometric_sum_eval
            )
          }, numeric(1)
        )
        
        grid_z_matrices[[sr]] <- grid_z_matrix
        
        if(show_idvl_landscapes){
          
          # plot the fitted landscape
          data3js <- r3js::surface3js(
            data3js,
            x = grid_x_matrix,
            y = grid_y_matrix,
            z = grid_z_matrix,
            col = "grey70",
            opacity = 0.3,
            toggle = sprintf("Idvl, %s", srg),
            wireframe = FALSE,
            doubleSide = FALSE
          )
        }
        
      
      }
      
    }
    
    if(calculate_gmt_from_idvl_landscapes){
      
      grid_z_matrix <- calculate_gmt_landscape_from_idvl_landscapes(grid_z_matrices)
      
    } else {
      
      cone_coords <- cone_vals[[srg]]$cone_coords
      cone_heights <- cone_vals[[srg]]$cone_heights
      
      # calculate landscape values
      grid_z_matrix[] <- vapply(
        seq_len(length(grid_z_matrix)),
        \(n) {
          fit_lndscp_val(
            grid_x_matrix[n], grid_y_matrix[n],
            cone_coords = cone_coords, 
            cone_heights = cone_heights,
            cone_slopes = cone_slopes, 
            sr_bias = sr_bias, 
            sum_cone_eval = cone_vals[[srg]]$sum_cone_eval,
            geometric_sum_eval = cone_vals[[srg]]$geometric_sum_eval
          )
        }, numeric(1)
      )
      
    }
    
    
    # plot the fitted landscape
    data3js <- r3js::surface3js(
      data3js,
      x = grid_x_matrix,
      y = grid_y_matrix,
      z = grid_z_matrix,
      col = lndscp_colors[srg, "Color"],
      opacity = ifelse(show_idvl_landscapes,1, 0.8),
      toggle = sprintf("GMT, %s", srg),
      wireframe = FALSE,
      doubleSide = FALSE
    )
    
    # options
    if(show_cone_coordinates){
      
      data3js <- r3js::points3js(
        data3js,
        x         = cone_coords[,1],
        y         = cone_coords[,2],
        z         = rep(0, nrow(cone_coords)),
        size      = 2,
        col  = lndscp_colors[srg, "Color"],
        toggle = sprintf("%s cone coordinates", srg),
        opacity   = 1,
        shape      = "triangle open"
      )
    }
    
    if(show_cone_heights){
      
      data3js <- r3js::points3js(
        data3js,
        x         = cone_coords[,1],
        y         = cone_coords[,2],
        z         = cone_heights,
        size      = 0.5,
        col  = lndscp_colors[srg, "Color"],
        toggle = sprintf("%s cone heights", srg),
        opacity   = 1,
        shape      = "sphere"
      )
      
    }
    
    if(show_fitted_values){
      data3js <- r3js::points3js(
        data3js,
        x         = ag_coords[,1],
        y         = ag_coords[,2],
        z         = combined_cone_fit,
        size      = 0.5,
        col  = lndscp_colors[srg, "Color"],
        toggle = sprintf("%s fitted titers", srg),
        opacity   = 1,
        shape      = "sphere"
      )
    }
    
    if(show_measured_values){
      
      if(is.null(gmt_data)){
        warning("You want to plot measured GMTs but have not provided them")
      }
      
      temp_gmt <- gmt_data %>%
        filter(sr_group == srg)
      
      temp_gmt <- temp_gmt[match(agNames(map), temp_gmt$ag_name),]
      
      for(ag in 1:nrow(ag_coords)){
        data3js <- r3js::lines3js(
          data3js,
          x = rep(ag_coords[ag,1], 2),
          y = rep(ag_coords[ag,2], 2),
          z = c(0, temp_gmt$logtiter[ag]),
          col = "grey50",
          toggle = sprintf("GMT, %s", srg),
          geometry = TRUE,
          opacity = 1,
          lwd = 0.2 
        )
      }
      
      
      data3js <- r3js::points3js(
        data3js,
        x         = ag_coords[,1],
        y         = ag_coords[,2],
        z         = temp_gmt$logtiter,
        size      = 1,
        col  = lndscp_colors[srg, "Color"],
        toggle = sprintf("GMT, %s", srg),
        opacity   = 1,
        shape      = "sphere"
      )
      
      
    }
      
      
  }
  
  data3js <- set_r3js_orentation(data3js)
  
  return(data3js)
    
  
}



#---------------- Bayesian simulation evaluation
# ------------------- Titer comparison
compare_simulated_fitted_titer <- function(titerdata_test, titerdata_true, fitted_titers, ag_order){
  
  titerdata_comb <- rbind(titerdata_test %>%
                            mutate(Data = "Simulation") %>%
                            mutate(sr_bias = as.numeric(sr_name)),
                          titerdata_true %>%
                            mutate(Data = "Simulation") %>%
                            mutate(sr_name = "GMT") %>%
                            mutate(sr_bias = 0),
                          fitted_titers %>%
                            mutate(Data = "Fit") %>%
                            mutate(logtiter_adjusted = mean_logtiter_fitted + sr_bias) %>%
                            select(sr_name, sr_bias, logtiter_adjusted, sr_group, ag_name, Data),
                          fitted_titers %>%
                            mutate(Data = "Fit") %>%
                            mutate(logtiter_adjusted = mean_logtiter_fitted) %>%
                            mutate(sr_name = "GMT",
                                   sr_bias = 0) %>%
                            select(sr_name, sr_bias, logtiter_adjusted, sr_group, ag_name, Data) %>%
                            unique()
  )
  
  titerdata_comb %>%
    mutate(GMT = sr_name == "GMT",
           sr_name = paste(sr_name, Data)) %>% 
    ggplot(aes(x = ag_name, y = logtiter_adjusted, group = sr_name, color = Data, alpha = GMT)) + 
    geom_point() + 
    geom_line() + 
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
    scale_x_discrete(limits = ag_order,
                     name = "Antigen variant") +
    ylab("Log2(Titer/10)") + 
    theme_bw() -> titers_overview
  
  titerdata_comb %>%
    group_by(sr_name, ag_name) %>%
    summarize(titer_diff = 2^logtiter_adjusted[Data == "Simulation"] - 2^logtiter_adjusted[Data == "Fit"],
              sr_bias_diff = 2^sr_bias[Data == "Simulation"] - 2^sr_bias[Data == "Fit"]) -> titer_diff
  
  titer_diff %>%
    filter(sr_name == "GMT") %>%
    ggplot(aes(x = ag_name, y = titer_diff)) + 
    geom_point() + 
    scale_x_discrete(limits = ag_order,
                     name = "Antigen variant") +
    ylab("True - fitted titer") +
    theme_bw() -> titer_diff_plot
  
  
  titer_diff %>%
    filter(sr_name != "GMT") %>%
    ggplot(aes(x = sr_name, y = sr_bias_diff)) + 
    geom_point() + 
    ylab("True - fitted serum bias") +
    theme_bw() -> sr_diff_plot
  
  titerdata_comb %>%
    ggplot(aes(x = ag_name, y = logtiter_adjusted, color = Data, group = Data)) + 
    geom_line() + 
    geom_point() + 
    theme_bw() + 
    scale_x_discrete(limits = ag_order,
                     name = "Antigen variant") +
    ylab("Log2(Titer/10)") +
    facet_wrap(~sr_name) -> sr_titers_plot
  
  titerdata_comb %>%
    ggplot(aes(x = sr_name, y = sr_bias, color = Data, group = Data)) + 
    geom_line() + 
    geom_point() + 
    ylab("Log2(Serum bias)") +
    theme_bw() -> sr_bias_plot
  
  
  return(list("titer_overview" = titers_overview,
              "titers_by_serum" = sr_titers_plot,
              "titers_differences" = titer_diff_plot,
              "sr_bias_difference" = sr_diff_plot,
              "sr_bias_per_serum" = sr_bias_plot,
              "combined_titerdata" = titerdata_comb
  ))
  
}

compare_measured_fitted_titer <- function(titerdata_test, titerdata_true, fitted_titers, ag_order){
  
  if(length(fitted_titers) == 1){
    
    fitted_titers <- fitted_titers[[1]]
    
    
  } else {
    
    fitted_titers <- lapply(names(fitted_titers), function(x){
      fitted_titers[[x]] %>%
        mutate(Data = x)
    })
    
    fitted_titers <- do.call(rbind, fitted_titers)
  }
  
  titerdata_comb <- rbind(titerdata_test %>%
                            mutate(Data = "Measured") %>%
                            select(sr_name, logtiter_adjusted, sr_group, ag_name, Data),
                          titerdata_true %>%
                            mutate(Data = "Measured") %>%
                            mutate(sr_name = "GMT") %>%
                            select(sr_name, logtiter_adjusted, sr_group, ag_name, Data) %>%
                            unique(),
                          fitted_titers %>%
                            mutate(logtiter_adjusted = mean_logtiter_fitted + sr_bias) %>%
                            select(sr_name, logtiter_adjusted, sr_group, ag_name, Data),
                          fitted_titers %>%
                            mutate(logtiter_adjusted = mean_logtiter_fitted) %>%
                            mutate(sr_name = "GMT") %>%
                            select(sr_name, logtiter_adjusted, sr_group, ag_name, Data) %>%
                            unique()
  )
  
  different_data <- unique(fitted_titers$Data)
  overview_by_data <- lapply(different_data, function(x){
    
    titerdata_comb %>%
      filter(Data %in% c("Measured", x)) %>%
      mutate(GMT = sr_name == "GMT",
             sr_name = paste(sr_name, Data)) %>% 
      ggplot(aes(x = ag_name, y = logtiter_adjusted, group = sr_name, color = Data, alpha = GMT)) + 
      geom_point() + 
      geom_line() + 
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
      scale_x_discrete(limits = ag_order,
                       name = "Antigen variant") +
      scale_y_continuous(name = "Titer",
                         limits = c(-2, 10),
                         breaks = c(0:10),
                         labels = function(x) round(2^x*10)) + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } )
  
  
  titerdata_comb %>%
    group_by(sr_name, ag_name) %>%
    summarize(titer_diff = 2^logtiter_adjusted[Data == "Measured"] - 2^logtiter_adjusted[Data == "Fit"]) -> titer_diff
  
  titer_diff %>%
    filter(sr_name == "GMT") %>%
    ggplot(aes(x = ag_name, y = titer_diff)) + 
    geom_point() + 
    scale_x_discrete(limits = ag_order,
                     name = "Antigen variant") +
    ylab("True - fitted titer") +
    theme_bw() -> titer_diff_plot
  
  sr_levels <- unique(titerdata_comb$sr_name)
  sr_levels <- c(sr_levels[grep("GMT", sr_levels)], sr_levels[!grepl("GMT", sr_levels)])
  
  titerdata_comb$sr_name <- factor(titerdata_comb$sr_name, levels = sr_levels)
  
  titerdata_comb %>%
    ggplot(aes(x = ag_name, y = logtiter_adjusted, color = Data, group = Data)) + 
    geom_line() + 
    geom_point() + 
    theme_bw() + 
    scale_x_discrete(limits = ag_order,
                     name = "Antigen variant") +
    facet_wrap(~sr_name) + 
    scale_y_continuous(name = "Titer",
                       limits = c(-2, 10),
                       breaks = c(0:10),
                       labels = function(x) round(2^x*10)) + 
    theme(strip.background.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))  -> sr_titers_plot
  
  
  titerdata_comb %>%
    filter(sr_name == "GMT") %>%
    ggplot(aes(x = ag_name, y = logtiter_adjusted, color = Data, group = Data)) + 
    geom_line() + 
    geom_point() + 
    theme_bw() + 
    scale_x_discrete(limits = ag_order,
                     name = "Antigen variant") +
    facet_wrap(~sr_group) + 
    scale_y_continuous(name = "Titer",
                       limits = c(-2, 10),
                       breaks = c(0:10),
                       labels = function(x) round(2^x*10)) + 
    theme(strip.background.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))-> only_gmts
  
  return(list("titer_overview" = overview_by_data,
              "titers_by_serum" = sr_titers_plot,
              "gmts_overview" = only_gmts,
              "titers_differences" = titer_diff_plot
  ))
  
}


plot_cone_comparisons <- function(map, fitted_cone, cone_centers, cone_heights, cone_slopes){
  
  n_cones <- nrow(fitted_cone$cone_coords)
  
  ag_coords_df <- as.data.frame(agCoords(map)) 
  for(n in 1:n_cones){
    ag_coords_df[paste0("True center", n),] <- cone_centers[n,]
    ag_coords_df[paste0("Fitted center", n),] <- fitted_cone$cone_coords[n,]
    
  }
  ag_coords_df %>%
    rownames_to_column(var = "Antigen") -> ag_coords_df
  
  ag_coords_df %>%
    mutate(map_coords = !grepl("center", Antigen),
           Antigen = gsub(paste(n_cones, collapse = "|"), "",Antigen)) -> ag_coords_df
  
  ag_colors <- c(agFill(map), "red", "blue")
  names(ag_colors) <- c(agNames(map), "True center", "Fitted center")
  
  ag_coords_df %>%
    ggplot(aes(x = V1, y = V2, fill = Antigen, size = map_coords, alpha = map_coords)) + 
    geom_point(shape = 21) +
    scale_alpha_manual(values = c("TRUE" = 0.6, "FALSE" = 1)) +
    scale_fill_manual(values = ag_colors) +
    xlab("x") + 
    ylab("y") +
    theme_bw() -> map_2d_plot
  
  cone_values_df <- data.frame("x" = c(cone_centers[,1], fitted_cone$cone_coords[,1]),
                               "y" = c(cone_centers[,2],  fitted_cone$cone_coords[,2]),
                               "slope" = c(cone_slopes, fitted_cone$cone_slopes),
                               "height" =  c(cone_heights, fitted_cone$cone_heights),
                               "cone_nr" = rep(c(1:n_cones), 2),
                               "Data" = c(rep("Simulation", n_cones), rep("Fit", n_cones)))
  
  cone_values_df %>%
    pivot_longer(cols = c("x", "y", "slope", "height"), names_to = "variable", values_to = "value") -> cone_values_df
  
  cone_values_df %>%
    group_by(variable, cone_nr) %>%
    mutate(Difference = value[Data == "Simulation"] - value[Data == "Fit"],
           norm_difference = abs(Difference)/abs(value[Data == "Simulation"])) -> cone_values_df
  
  diff_threshold <- 0.7
  cone_values_df %>%
    ggplot(aes(x = variable, y = value, color = Data, group = factor(cone_nr))) + 
    geom_point(position = position_dodge(width = 0.4)) + 
    geom_point(data = cone_values_df %>% filter(Data == "Fit"), aes(y = Difference), position = position_dodge(width = 0.4), size = 2, shape = 24) +
    theme_bw() -> cone_values_plot
  
  return(list("cones_in_map" = map_2d_plot,
              "cone_values_plot" = cone_values_plot))
  
}
