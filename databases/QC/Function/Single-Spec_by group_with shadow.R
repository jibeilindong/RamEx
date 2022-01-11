
Spec_shadow <- function(hyperspec,Meta_group)
{
  hyperspec_frame <- cbind(as.data.frame(hyperspec)[,1:6],as.data.frame(hyperspec$spc))
  hyperspec_melt <- spc_melt(hyperspec_frame,Meta_group,450,3000)
  hyperspec_melt_summary <- summarySE(hyperspec_melt, measurevar = "value", groupvars = c(Meta_group,"wavenumber"))
  #hyperspec_melt_summary$group_C <- factor(hyperspec_melt_summary$group_C,levels = c("0h","2h","4h","8h","12h"))
  #hyperspec_melt_summary <- filter(hyperspec_melt_summary,group_B != "bg")
  
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot_hyperspec <- ggplot(data = hyperspec_melt_summary,aes(x = wavenumber, y = value, group = group_B)) +
    geom_ribbon(aes(ymin = value - sd, ymax = value + sd, fill = group_B),alpha = 0.5) +
    geom_line(aes(color = group_B),size = 0.8) + theme_bw() + facet_grid(group_C ~ group_B) +
    labs(y = "Normalized Intensity (a.u.)") + xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15,family = "myFont"),
      #legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color = "black"),
      axis.title.y = element_text(size = 20, angle = 90,family = "myFont"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      strip.text = element_markdown(size = 20,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title = element_text(family = "myFont" ,size = 20))
  plot_hyperspec
  return(plot_hyperspec)
}
