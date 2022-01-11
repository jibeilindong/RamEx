Spec_all <- function(hyperspec,Meta_group)
{
  hyperspec_frame <- cbind(as.data.frame(hyperspec)[,1:5],as.data.frame(hyperspec$spc))
  hyperspec_frame$Number <- c(1:nrow(hyperspec_frame))
  hyperspec_melt <- spc_melt(hyperspec_frame,c(Meta_group,"Number"),450,3200)
  
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot_hyperspec <- ggplot(data = hyperspec_melt,aes(x = wavenumber, y = value, group = as.factor(Number))) +
    geom_line(aes(color = as.factor(Number)),size = 0.8) + theme_bw() + facet_wrap( . ~ group_A,nrow = 1) +#,
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
  plot_hyperspec
  return(plot_hyperspec)
}


Spec_Nor_all <- function(hyperspec,Meta_group)
{
  hyperspec_frame <- cbind(as.data.frame(hyperspec)[,1:5],as.data.frame(hyperspec$spc))
  hyperspec_frame$Number <- c(1:nrow(hyperspec_frame))
  hyperspec_melt <- spc_melt(hyperspec_frame,c(Meta_group,"Number"),450,3200)

  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot_hyperspec <- ggplot(data = hyperspec_melt,aes(x = wavenumber, y = value, group = as.factor(Number))) +
    geom_line(aes(color = as.factor(Number)),size = 0.8) + theme_bw() + facet_wrap( . ~ group_A,nrow = 1) +#,
    labs(y = "Normalized Intensity (a.u.)") + xlab(expression(paste('Wavenumber (cm'^{-1},')',sep = ""))) +
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
  plot_hyperspec
  return(plot_hyperspec)
}


tsne_plot <- function(tSNE_data)
{
  tSNE_data$label <- gsub("24433","*C. albicans*",tSNE_data$label)
  tSNE_data$label <- gsub("6258","*C. krusei*",tSNE_data$label)
  tSNE_data$label <- gsub("C auris","*C. auris*",tSNE_data$label)
  tSNE_data$label <- gsub("C glabrata","*C. glabrata*",tSNE_data$label)
  tSNE_data$label <- gsub("C parapsilosis","*C. parapsilosis*",tSNE_data$label)
  tSNE_data$label <- gsub("C tropicalis","*C. tropicalis*",tSNE_data$label)
  
  
#windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"
  plot <- ggplot(tSNE_data,aes(tSNE1,tSNE2,color = as.factor(label)))+ 
    geom_point() + theme_bw() + stat_ellipse(level = 0.8)+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      #legend.position = c(0.85,0.9),
      legend.text = element_markdown(size = 15,family = "myFont"),
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "myFont" ,size=15))
  return(plot)
}

###PCA
pca <- function(data,label)
{
  library("ggbiplot")
  data_pca <- prcomp(data,cor=T,scores=T,scale. = TRUE)
  return(data_pca)
}


###pca画图
pca_plot <- function(data.pca,label)
{
  #windowsFonts(myFont = windowsFont("Times New Roman"))
  myFont="Times New Roman"

  p <- ggbiplot(data.pca, obs.scale = 2, var.scale = 1,var.axes = FALSE,
                groups = label, ellipse = TRUE, circle = TRUE) +
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      legend.text = element_markdown(size = 15,family = "myFont"),
      #legend.position = "none",
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "myFont" ,size=15))
  return(p)
}
###

PLSDA <- function(data,label)
{
  library("mixOmics")
  plsda.res <- plsda(data, label, ncomp = 20) 
  plotIndiv(plsda.res, ind.names = F, legend = TRUE, ellipse = TRUE,title = 'PLS-DA')
  return(plsda.res) 
}
###plsda_plot
plsda_plot <- function(plsda_data,name,label)
{
  plsda_data_variates<- plsda_data$variates[[1]]
  plsda_data_variates <- data.frame(plsda_data_variates)
  plsda_data_variates <- cbind(label,plsda_data_variates)
  
#windowsFonts(myFont = windowsFont("Times New Roman"))a
  myFont="Times New Roman"

  plot1  <- ggplot(plsda_data_variates,aes(comp1,comp2,color=label)) +
    geom_point()+ stat_ellipse(level = 0.8)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_blank(),
      #legend.position = c(0.85,0.9),
      legend.text = element_markdown(size = 15,family = "myFont"),
      legend.background = element_blank(),
      text = element_text(color="black"),
      axis.text.x = element_text(size = 15,angle = 0,family = "myFont"),
      axis.text.y = element_text(size = 15,family = "myFont"),
      axis.ticks = element_line(size = 1),
      axis.ticks.length = unit(0.4,"lines"),
      axis.title= element_text(family = "myFont" ,size=15))
  return(plot1)
}

