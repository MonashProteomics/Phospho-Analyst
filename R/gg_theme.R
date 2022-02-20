# ggplot theme 1
theme_DEP1 <- function() {
  # Use theme_bw() as default
  basesize <- 12
  theme <- ggplot2::theme_bw(base_size = basesize)
  
  # Change plot title appearance
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 1
  theme$plot.title$hjust <- 0.5
  
  # Change axis title appearance
  theme$axis.title.x$size <- basesize + 1
  
  theme$axis.title.y$size <- basesize + 1
  
  # Change axis text appearance
  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"
  
  # Change legend title appearance
  theme$legend.title$size <- basesize + 1
  
  # Change legend text appearance
  theme$legend.text$size <- basesize
  
  # Change strip text (facet headers) appearance
  theme$strip.text$face <- "plain"
  # theme$strip.text$size <- basesize + 1
  theme$strip.text$colour <- "black"
  
  return(theme)
}

# ggplot theme 2
theme_DEP2 <- function() {
  # Get vertical x-axis labels
  theme <- theme_DEP1()
  theme$axis.text.x$angle <- 90
  theme$axis.text.x$hjust <- 1
  theme$axis.text.x$vjust <- 0.5
  return(theme)
}