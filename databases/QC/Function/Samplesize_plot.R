Samplesize_plot <- function(Pvalue_data)
{
  plot<-ggplot(Pvalue_data,aes(x=Number_sample,y=Pvalue,color = N))+
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
    geom_point(aes(group=N,color = N))+
    facet_grid( Group ~ N,scales = "free_x" )+
    geom_hline(aes(yintercept=0.95), colour="red", linetype="dashed")+
    geom_vline(aes(xintercept=Number_sample_max_Pvalue), colour="red", linetype="dashed")+
    geom_text(aes(x=Number_sample_max_Pvalue,y=0.74,label=paste(Number_sample_max_Pvalue+20,"cells",sep=" ")),colour="#990000",size=5)+
    xlab("Sample size")+
    ylab(expression(paste('Probability (RE'[mean],'<0.95)',sep="")) )
  return(plot)
}
