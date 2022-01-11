default_theme <- function()
{
  windowsFonts(myFont = windowsFont("Times New Roman"))
  theme_bw()+
  theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 17,family = "myFont"),
      legend.text = element_text(size = 15,family = "myFont"),
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
}

library("rlang")
draw_spc <- function(data,name,color_text = "Drug(μg/mL)",group,color,facet_grid_x,facet_grid_y,weight = 12,height = 4,outpath)
{
  if (facet_grid_x == ".")
  {
    if (facet_grid_y == ".")
    {
      plot_hyperspec1 <- ggplot(data = data,aes(x = wavenumber, y = value,group= {{group}},color={{color}}))+
        #facet_grid( B ~ .,scales = "free")+
        geom_line(size = 1)+ 
        labs(
          x = expression(paste('Wave Number (cm'^{-1},')',sep="")),
          y = "Intensity(a.u.)",
          color = color_text)+
        default_theme()
    }else
    {
      plot_hyperspec1 <- ggplot(data = data,aes(x = wavenumber, y = value,group= {{group}},color={{color}}))+
        facet_grid( . ~ {{facet_grid_y}},scales = "free")+
        geom_line(size = 1)+ 
        labs(
          x = expression(paste('Wave Number (cm'^{-1},')',sep="")),
          y = "Intensity(a.u.)",
          color = color_text)+
        default_theme()
    }
  }else if (facet_grid_y == ".")
  {
    plot_hyperspec1 <- ggplot(data = data,aes(x = wavenumber, y = value,group= {{group}},color={{color}}))+
      facet_grid( {{facet_grid_x}} ~ .,scales = "free")+
      geom_line(size = 1)+ 
      labs(
        x = expression(paste('Wave Number (cm'^{-1},')',sep="")),
        y = "Intensity(a.u.)",
        color = color_text)+
      default_theme()
  }else
  {
    plot_hyperspec1 <- ggplot(data = data,aes(x = wavenumber, y = value,group= {{group}},color={{color}}))+
      facet_grid( {{facet_grid_x}} ~ {{facet_grid_y}},scales = "free")+
      geom_line(size = 1)+ 
      labs(
        x = expression(paste('Wave Number (cm'^{-1},')',sep="")),
        y = "Intensity(a.u.)",
        color = color_text)+
      default_theme()    
  }
    ggsave(filename=paste(outpath,"//",name,"_spc.wmf", sep=""),plot=plot_hyperspec1,
           limitsize=T,width = weight,height = height)
}

###画箱式图-CDR
draw_boxplot <- function(data,name,xlab_title,ylab_title = "CD/(CD+CH)",color = "Drug(μg/mL)",width = 8,height = 4,outpath)
{
  windowsFonts(myFont = windowsFont("Times New Roman"))
  plot_hyperspec1 <- ggplot(data = data,aes(x = C, y = CD_ratio))+
    facet_grid( . ~ A)+
    geom_boxplot(aes(color=B)) + 
    stat_smooth(aes(group=B,color=B),method = "loess",size=1.5,level = 0.9,position=position_dodge(0.2))+
    labs(
        x = xlab_title,
        y = ylab_title,
        color = color)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 17,family = "myFont"),
      legend.text = element_text(size = 15,family = "myFont"),
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.title.y = element_text(size = 20, angle = 90,family = "myFont"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      strip.text = element_text(size = 20,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "myFont" ,size=20))
  plot_hyperspec1 
  ggsave(filename=paste(outpath,"//",name,"_CDR_boxplot.png", sep=""),plot=plot_hyperspec1, limitsize=T,width = width,height = height)
}

###画密度图-CDR
draw_density <- function(data,name,xlab_title = "CD/(CD+CH)",ylab_title,color = "Drug(μg/mL)",width = 10,height = 4,outpath)
{
    windowsFonts(myFont = windowsFont("Times New Roman"))
    plot_hyperspec1 <- ggplot(data = data,aes(x = CD_ratio,group = C,fill =C))+
        facet_grid( C ~ A)+
        geom_density() + 
        labs(
            x = xlab_title,
            y = ylab_title,
            color = color)+
        #scale_x_continuous(breaks = seq(2050,2300,80),limits = c(2050,2300))+
        theme_bw()+
        theme(
            panel.grid = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_text(size = 17,family = "myFont"),
            legend.text = element_text(size = 15,family = "myFont"),
            legend.background = element_blank(),
            text = element_text(color="black"),
            axis.title.y = element_text(size = 20, angle = 90,family = "myFont"),
            axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
            axis.text.y = element_text(size = 15,family = "myFont"),
            strip.text = element_text(size = 20,family = "myFont"),
            axis.ticks = element_line(size = 1),
            axis.ticks.length = unit(0.4,"lines"),
            axis.title= element_text(family = "myFont" ,size=20))
    plot_hyperspec1 
  ggsave(filename=paste(testpath,"//",name,"_CDR_boxplot.png", sep=""),plot=plot_hyperspec1, limitsize=T,width = width,height = height)  
}
 