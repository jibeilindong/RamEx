Real_time_box_plot <- function(alldata_n)
{
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot<-ggplot(alldata_n,aes(x=Group,y=CD_ratio,color = Group))+
    geom_boxplot(aes(group=alldata_n$Group),alpha=1/5)+#
    geom_point(aes(group=Group,colour=Group))+
    facet_grid( . ~ N,scales = "free_x" )+
    stat_smooth(aes(group=N),method = "loess",size=1,level = 0.9,col = "black")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15,family = "myFont"),
      legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color = "black"),
      axis.title.y = element_text(size = 20, angle = 90,family = "myFont"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 15,angle = 45,family = "myFont",hjust = 1,vjust = 1),
      axis.text.y = element_text(size = 15,family = "myFont"),
      strip.background = element_rect(fill = "white"),
      strip.text = element_markdown(size = 20,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title = element_text(family = "myFont" ,size = 20))+
    ylab('CD-ratio Value')
  return(plot)
}
