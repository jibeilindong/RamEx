###
Spec_all <- function(hyperspec)
{
  hyperspec_frame <- cbind(dplyr::select(as.data.frame(hyperspec@data),-spc),as.data.frame(hyperspec$spc))
  hyperspec_frame$Number <- 1:nrow(hyperspec_frame)
  hyperspec_melt <- spc_melt(hyperspec_frame,c("Group","Number"),floor(min(hyperspec@wavelength)),ceiling(max(hyperspec@wavelength)))
  
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot_hyperspec <- ggplot(data = hyperspec_melt,aes(x = wavenumber, y = value, group = as.factor(Number))) +
    geom_line(aes(color = as.factor(Number)),size = 0.8) + theme_bw() + facet_wrap(. ~ Group,nrow = 1) +#,
    labs(y = "Intensity (cnt)") + xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15,family = "myFont"),
      legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color = "black"),
      axis.title.y = element_text(size = 20, angle = 90,family = "myFont"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_markdown(size = 20,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title = element_text(family = "myFont" ,size = 20))
  return(plot_hyperspec)
}
