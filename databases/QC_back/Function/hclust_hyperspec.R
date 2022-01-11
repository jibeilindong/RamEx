hclust_hyperspec_drawing <- function(data,result_hc_ID,name,outpath)
{
  data$group <- result_hc_ID
  hyperspec <- data
  Meta_group <- c("group_A","group_B","group")
  hyperspec_frame <- cbind(as.data.frame(hyperspec)[,c(1:6,9)],as.data.frame(hyperspec$spc))
  hyperspec_melt <- spc_melt(hyperspec_frame,Meta_group,450,3000)
  hyperspec_melt_summary <- summarySE(hyperspec_melt, measurevar = "value", groupvars = c(Meta_group,"wavenumber"))
  return_data <- select(hyperspec_frame,filename,group_A,group_C,group)
  
  windowsFonts(myFont = windowsFont("Times New Roman"))
  plot_hyperspec <- ggplot(data = hyperspec_melt_summary,aes(x = wavenumber, y = value, group = group_A,color = group_A)) +
    geom_ribbon(aes(ymin = value - sd, ymax = value + sd, fill = group_A),alpha = 0.5) +
    geom_line(size = 0.8) + theme_bw() + facet_grid(group ~ .) +#,nrow = 3,ncol = 5
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
      strip.text = element_text(size = 20,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title = element_text(family = "myFont" ,size = 20))
  plot_hyperspec
  ggsave(filename=paste(outpath,name,"_hclust-hyperspec.png", sep=""),plot=plot_hyperspec, limitsize=T,width = 8,height = 15)
  return(return_data)
}
