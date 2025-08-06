# label map

ggplot(map) +
  ggrepel::geom_label_repel(
    data = filter(gp$data, type == "AG"),
    mapping = aes(
      x = x,
      y = y,
      label = text,
      point.size = size
    ),
    min.segment.length = 0,
    label.padding = 0.2,
    label.size = NA,
    fill = NA
  )