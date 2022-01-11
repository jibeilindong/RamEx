
###预载函数
Sampling_depth <- function(CD_ratio_group)
{
  Time_sample <- 1000
  Pvalue <- Mean_sample <- c()
  if (length(CD_ratio_group) < 300)
  { Number_sample <- c(1:length(CD_ratio_group)) }
  else {Number_sample <- c(1:280) }
  Sample_1 <- sample(CD_ratio_groups$CD_ratio,50)
  for (n in Number_sample )#length(CD_ratio_groups$CD_ratio)
  {
    Sample_x <- c(Sample_1,sample(CD_ratio_groups$CD_ratio,n))
    New_Mean_sample_x <- mean(Sample_x)
    for (j in 1:Time_sample)
    { 
      Sample_y <- sample(Sample_x,n + 19)
      factor_a <- mean(Sample_y)
      Mean_sample <- c(Mean_sample,factor_a)
    }
    Pvalue <- c(Pvalue,sum(abs(1 - Mean_sample/New_Mean_sample_x) < 0.05)/Time_sample)
    Mean_sample <- c()
  }
  return(Pvalue)
}

###画平均光谱
Figure_1_A_mean_spc_draw <- function(cluster.meansd.dataframe_melt,name,outpath,width,height)
{
  #barplot( 1: 9, col = Color_3h_persister) 
  Color_1_Candida <- brewer.pal( 9, 'YlGn')[c(9,8,7,6,5,4)]
  Color_2_Candida <- brewer.pal( 9, 'YlOrRd')[c(3,5,9)]
  Color_4h_Candida <- c(Color_1_Candida,Color_2_Candida)
  
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot_hyperspec <- ggplot(data = cluster.meansd.dataframe_melt,aes(x = wavenumber, y = value, group= factor(V2),color= factor(V2)))+
    geom_line(size = 0.8)+ theme_bw()+ 
    labs(y = "Normalization Intensity (a.u.)")+   xlab(expression(paste('Wavenumber (cm'^{-1},')',sep="")))+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15,family = "myFont"),
      #legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.title.y = element_text(size = 20, angle = 90,family = "myFont"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      strip.text = element_text(size = 20,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "myFont" ,size=20))
  plot_hyperspec
  ggsave(filename=paste(outpath,"//",name,".png", sep=""),plot=plot_hyperspec,
         limitsize=T,width = width,height = height,dpi = 1200)
}

###画平均光谱
Figure_1_zoom_mean_spc_draw <- function(cluster.meansd.dataframe_melt,name,outpath,width,height)
{
  for(i in levels(factor(cluster.meansd.dataframe_melt$V1)))
  {
    if (i == levels(factor(cluster.meansd.dataframe_melt$V1))[1])
    {
      Color <- Color_4h_Candida[1:6]
    }else
    {
      Color <- Color_4h_Candida[c(1,5,6,7,8,9)]
    }
    cluster.meansd.dataframe_melt_zoom <- filter(cluster.meansd.dataframe_melt,V1 == i,wavenumber < 2500)
    #windowsFonts(myFont = windowsFont("Times New Roman"))
    myFont="Times New Roman"
    plot_hyperspec <- ggplot(data = cluster.meansd.dataframe_melt_zoom,aes(x = wavenumber, y = value, group= factor(V2),color= factor(V2)))+
      geom_line(size = 0.8)+ theme_bw()+ 
      scale_color_manual("Micafungin (μg/mL)",values = Color,aesthetics = "color")+
      scale_y_continuous(breaks = seq(0,0.3,0.1),limits = c(-0.01,0.3))+
      theme(
        panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_markdown(size = 15,family = "myFont"),
        legend.text = element_markdown(size = 15,family = "myFont"),
        legend.position = "none",
        legend.background = element_blank(),
        text = element_text(color="black"),
        axis.line = element_line(size = 0.5),
        axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
        axis.text.y = element_text(size = 15,family = "myFont"),
        strip.text=element_blank(),
        axis.ticks = element_line(size = 1),
        axis.ticks.length = unit(0.4,"lines"),
        axis.title= element_blank())
    plot_hyperspec
    ggsave(filename=paste(outpath,"//",i,"_",name,".tiff", sep=""),plot=plot_hyperspec,
           limitsize=T,width = width,height = height,dpi = 1200)
  }
  
}

growth_curve <- function(all_data_summary_Drug_i,Color_4h_Candida)
{
  plot_growth_curve <- 
    ggplot(all_data_summary_Drug_i,aes(x=Time,y=ODvalue))+ 
    geom_ribbon(aes(ymin = ODvalue-se, ymax = ODvalue+se,fill = Conc),alpha =0.2)+
    geom_line(aes(colour=Conc),size=1.5)+
    scale_fill_manual("Concentration of Drug (μg/mL)",values = Color_4h_Candida,aesthetics = c("fill","color"))+ #
    #scale_x_continuous(breaks = seq(10,50,5),limits = c(10,45))+
    theme(axis.ticks = element_blank())+  theme_bw()+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank()) +  
    theme(axis.text = element_text(size=15))+ theme(legend.text = element_text(size=15))+
    theme(legend.title = element_text(size=20))+
    theme(axis.title= element_text( family = "myFont" ,size=20))+
    ylab(expression('OD'[600]))+
    xlab(expression(paste('Time(h)',sep="")))+
    theme(legend.title = element_text(size=15))#+
  #theme(legend.justification=c(0.95,0.05),legend.position=c(0.95,0.05))
  return(plot_growth_curve)
}

growth_curve_D2O <- function(all_data_summary_Drug_i,Color_4h_Candida)
{
  plot_growth_curve <- 
    ggplot(all_data_summary_Drug_i,aes(x=Time,y=ODvalue))+ 
    geom_ribbon(aes(ymin = ODvalue-se, ymax = ODvalue+se,fill = Conc),alpha =0.2)+
    geom_line(aes(colour=Conc),size=1.5)+ facet_grid(Cellconc ~ D2O,scales = "free")+
    scale_fill_manual("Concentration of Drug (μg/mL)",values = Color_4h_Candida,aesthetics = c("fill","color"))+ #
    #scale_x_continuous(breaks = seq(10,50,5),limits = c(10,45))+
    theme(axis.ticks = element_blank())+  theme_bw()+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank()) +  
    theme(axis.text = element_text(size=15))+ theme(legend.text = element_text(size=15))+
    theme(legend.title = element_text(size=20))+
    theme(axis.title= element_text( family = "myFont" ,size=20))+
    ylab(expression('OD'[600]))+
    xlab(expression(paste('Time(h)',sep="")))+
    theme(legend.title = element_text(size=15))#+
  #theme(legend.justification=c(0.95,0.05),legend.position=c(0.95,0.05))
  return(plot_growth_curve)
}
###tSNE_plot
tSNE_plot_pdf <- function(tsne_plot,name,width = 10,height = 5,outpath)
{
  plot <- ggplot(tsne_plot,aes(tSNE1,tSNE2,color=label)) + 
    geom_point() + theme_bw() + stat_ellipse(level = 0.9)+
    #scale_color_manual(values = c("#7CAE00","#00BFC4","#F8766D","#C77CFF"))+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.position = "none",
      legend.text = element_markdown(size = 15,family = "Times"),
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.text.x = element_text(size = 15,angle = 0,family = "Times"),
      axis.text.y = element_text(size = 15,family = "Times"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "Times" ,size=15))
  ggsave(filename=paste(outpath,"//",name,"_tSNE.pdf", sep=""),plot=plot,
         limitsize=T,width = width,height = height,dpi = 1200)
}

###tSNE_plot
tSNE_plot <- function(tsne_plot,name,width = 10,height = 5,outpath)
{
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot <- ggplot(tsne_plot,aes(tSNE1,tSNE2,color=label)) + 
    geom_point() + theme_bw() + 
    stat_ellipse(level = 0.9)+
    scale_color_manual(values = c("#7CAE00","#00BFC4","#F8766D","#C77CFF"))+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15,family = "myFont"),
      #legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.title.y = element_text(size = 20, angle = 90,family = "myFont"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      strip.text = element_text(size = 20,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "myFont" ,size=20))
  ggsave(filename=paste(outpath,"//",name,"_tSNE.tiff", sep=""),plot=plot, limitsize=T,width = width,height = height,dpi = 1200)
}


