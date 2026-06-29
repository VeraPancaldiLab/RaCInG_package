library(ggplot2)
library(hexSticker)
library(showtext)
library(sysfonts)

font_add_google("Roboto Condensed", "roboto")
showtext_auto()

set.seed(42)

# Cell colors (different cell types communicating)
col_blue   <- "#2e86c1"
col_green  <- "#27ae60"
col_red    <- "#cb4335"
col_purple <- "#8e44ad"
col_orange <- "#e67e22"

# Cell positions in a network layout
cells <- data.frame(
  x    = c(0.5, 1.5, 1.8, 1.2, 0.2),
  y    = c(0.85, 0.85, 0.45, 0.15, 0.45),
  size = c(5, 5, 4.5, 4.5, 4.5),
  type = c("A", "B", "C", "D", "E")
)

cell_colors <- c(
  "A" = col_blue,
  "B" = col_green,
  "C" = col_red,
  "D" = col_purple,
  "E" = col_orange
)

# Communication edges (directed arrows between cells)
edges <- data.frame(
  x    = c(0.5, 1.5, 1.8, 1.2, 0.2, 1.5, 0.5),
  xend = c(1.5, 1.8, 1.2, 0.2, 0.5, 0.2, 1.2),
  y    = c(0.85, 0.85, 0.45, 0.15, 0.45, 0.85, 0.85),
  yend = c(0.85, 0.45, 0.15, 0.45, 0.85, 0.45, 0.15)
)

p <- ggplot() +
  geom_segment(data = edges, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
               linewidth = 0.6, color = "#f39c12", alpha = 0.6) +
  geom_point(data = cells, aes(x = x, y = y, color = type),
             size = cells$size, alpha = 0.9, show.legend = FALSE) +
  geom_point(data = cells, aes(x = x, y = y),
             size = cells$size * 0.3, color = "white", alpha = 0.4, show.legend = FALSE) +
  scale_color_manual(values = cell_colors) +
  xlim(-0.1, 2.1) +
  ylim(-0.05, 1.05) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

sticker(
  subplot = p,
  package = "RaCInG",
  p_size = 17,
  p_y = 1.52,
  p_color = "#FFFFFF",
  p_family = "roboto",
  p_fontface = "bold",
  s_x = 1.0,
  s_y = 0.82,
  s_width = 1.7,
  s_height = 1.05,
  h_fill = "#1a1a2e",
  h_color = "#2e86c1",
  h_size = 1.8,
  url = "",
  u_size = 0,
  filename = "man/figures/logo.png",
  dpi = 300
)

cat("Logo saved to man/figures/logo.png\n")
