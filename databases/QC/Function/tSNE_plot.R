tSNE_drawing <- function(data)
{
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot <- ggplot(data,aes(tSNE1,tSNE2,color = as.factor(label)))+ 
    geom_point() + theme_bw() + stat_ellipse(level = 0.8)+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15,family = "myFont"),
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "myFont" ,size=15))
}
