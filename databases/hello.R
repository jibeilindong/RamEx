# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}
###library_all
library_all <- function()
{
  library("MASS")
  library("lattice")
  library("RColorBrewer")
  library("dplyr")
  library("mixOmics")
  library('optparse')
  library("pROC")
  library("R.utils")
  library("ggplot2")
  library("reshape2")
  library("tidyverse")
  library("hyperSpec")
  library("ggpubr")
  library("data.table")
  library("scales")
  library("Rtsne")
  library("pls")
  library("e1071")
  library("ggbiplot")
  library("ggtext")
  library("prospectr")
  library("Rmisc")
}
###baseline
baseline <- function(data,poly.order = 3)
{
  wavelength <- as.numeric(colnames(data))
  data_hyperSpec <- new ("hyperSpec", spc = data, wavelength = wavelength)
  data_baseline <- data_hyperSpec-
    spc.fit.poly.below (data_hyperSpec, data_hyperSpec, poly.order = poly.order)
  data_baseline.data.frame <- as.data.frame(data_baseline$spc)
  plotspc(data_baseline)#,spc.nmax = nrow(data_baseline)
  return(data_baseline.data.frame)
}


###CDR
CD_ratio <- function(hyperspec)
{
  meta_data <- select(as.data.frame(hyperspec),-spc,-.row)
  data_baseline_normalization_frame  <- cbind(meta_data,as.data.frame(hyperspec$spc))
  wave_nums <- as.numeric(as.character(colnames(data_baseline_normalization_frame)))

  CD_start <- which.min(abs(wave_nums - 2050))
  CD_end <- which.min(abs(wave_nums - 2300))
  CH_start <- which.min(abs(wave_nums - 3050))
  CH_end <- which.min(abs(wave_nums - 2800))
  #CD_ratio<-rowSums(CD)/(rowSums(CD)+rowSums(CH)) #参照
  data_baseline_normalization_frame_CD_ratio <- mutate(data_baseline_normalization_frame,
                                                       CD_ratio = rowSums(data_baseline_normalization_frame[,CD_start:CD_end])/(rowSums(data_baseline_normalization_frame[,CD_start:CD_end]) + rowSums(data_baseline_normalization_frame[,CH_start:CH_end])))
  CD_ratio_groups <- select(data_baseline_normalization_frame_CD_ratio,CD_ratio)
  CD_ratio_groups <- cbind(meta_data,CD_ratio_groups)
  return(CD_ratio_groups)
}

###folder reader
##批量读取文件夹数据####
folder_reader <- function(foldernames)
{
  meta_data <- c()
  setwd(foldernames)
  filefolder <- list.files()
  workingpath <- paste(foldernames,"/",filefolder[1],sep = "")
  setwd(workingpath)
  filenames <- list.files(pattern = "*.txt")
  final_data <- spectrum <- as.data.frame(fread(filenames[1], header = F, sep = "\t")[,1])
  for (foldername in filefolder)
  {
    workingpath <- paste(foldernames,"/",foldername,sep = "")
    setwd(workingpath)
    filenames <- list.files(pattern = "*.txt")
    everyfolder_data <- spectrum
    meta_data_n <- filenames
    meta_data <-  c(meta_data,meta_data_n)
    for (filename in filenames)
    {
      data <- as.data.frame(fread(filename,header = F, sep = "\t"))
      if (nrow(everyfolder_data) != nrow(data))
      {
        rowname_everyfolder_data <- as.numeric(everyfolder_data[,1])
        rowname_data <- as.numeric(data[,1])
        new_rowname <- c()
        for (i in rowname_everyfolder_data)
        {
          x <- which.min(abs(rowname_data - i))
          new_rowname <- c(new_rowname,x)
        }
        data <- data[new_rowname,]
      }
      everyfolder_data <- cbind(everyfolder_data,data[,2])
    }
    final_data <- cbind(final_data,everyfolder_data[,-1])
  }
  final_dataframe <- as.data.frame(t(final_data))
  meta_data <- as.data.frame(meta_data)
  final_dataframe <- data.frame(filename = meta_data,final_dataframe[-1,])
  colnames(final_dataframe) <- c("filename",as.character(spectrum[,1]))
  return(final_dataframe)
}
###getspc
###获取特定范围的光谱——dataframe
spc_waveumber <- function(spc,from,to)
{
  wave_nums <- as.numeric(as.character(colnames(spc)))
  start <- which.min(abs(wave_nums-from))
  end <- which.min(abs(wave_nums-to))
  final_spc <- spc[,start:end]
  return(final_spc)
}


###hcluster_hyperspec
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
###hyperspec_to_dataframe
hyperspec_to_dataframe <- function(hyperspec,from,to)
{
  data_frame <- cbind(as.data.frame(hyperspec)[,from:to],as.data.frame(hyperspec$spc))
  return(data_frame)
}
###initialization

initialization <- function(working_dir,output_folder = "output")
{
  outpath <- paste(working_dir,"\\",output_folder,sep = "")
  setwd(working_dir)
  dir.create(outpath)
  print(paste("当前的工作路径是: ",working_dir,sep = ""))
  print(paste("当前的输出路径是: ",outpath,sep = ""))
}

###install_packages_for_me
install_packages_for_me <- function(used_packages = c("MASS","lattice","RColorBrewer","dplyr","mixOmics","optparse","pROC","R.utils","ggplot2",
                                                      "reshape2","tidyverse","hyperSpec","ggpubr","data.table","scales","Rtsne","pls","e1071",
                                                      "ggbiplot","ggtext","prospectr","devtools"))
{
  already_installed_packages <- library()$results[,1]
  install_pack <- setdiff(used_packages, already_installed_packages)
  for ( i in install_pack)
  {
    install.packages(i)
  }

  for ( i in used_packages)
  {
    library(i,character.only = T)
  }
  install_github("vqv/ggbiplot")
  install_github("mixOmicsTeam/mixOmics")
  library(ggbiplot)
  library(mixOmics)
}
###Meanspc

Mean.spc <- function(pretreatment.spc,Group,save = TRUE,filelabel,outpath)
{
  ###计算处理后平均光谱
  clusters <- list(Group)
  pretreatment.spc <- cbind(clusters,pretreatment.spc)
  cluster.meansd <- aggregate(pretreatment.spc,clusters,mean)
  cluster.meansd.dataframe <- as.data.frame(cluster.meansd)

  if (save == TRUE)
  {
    ###保存单个文件-平均光谱
    setwd(outpath)
    dir.create(paste(filelabel,"Mean_by_group",sep = ""))
    for (i in 1:nrow(cluster.meansd.dataframe))
    {
      cell <- data.frame(t(cluster.meansd.dataframe[i,-1]))
      fwrite(cell,
             paste(outpath,"/",filelabel,"Mean_by_group/",cluster.meansd.dataframe$Group.1[i],"_.txt",sep=""),
             row.names=T,col.names=F,quote=F,sep = "\t")
    }
  }
  return(cluster.meansd.dataframe)
}

###newhyperspec
new_hyperSpec <- function(spc,Name_group)
{
  filename <- spc[,1]
  meta_data <-  as.data.frame(tstrsplit(gsub(".txt","",spc[,1]),"_"),col.names = Name_group)
  spectrum <- as.numeric(colnames(spc[,-1]))
  data_hyperSpec<-new ("hyperSpec", data = meta_data,spc = spc,wavelength=spectrum)
  return(data_hyperSpec)
}

###Normalization
normalization <- function(hyperSpec,from = 2500,to = 3000,cal = "max")
{
  if (cal == "max")
  { factors <- 1/apply(hyperSpec[,,from~to],1,max)}  # 以CH的均值 来归一化
  else if (cal == "sum")
  { factors <- 1/apply(hyperSpec,1,sum)  }
  else
  { print("Error: cal should be max or sum") }
  data_normalization <- sweep(hyperSpec,1,factors,"*")
  return(data_normalization)
}

###Renishaw单光谱

wavelength_summary <- function(hyperspec,Name_group)
{
  meta <- as.data.frame(tstrsplit(gsub(".txt","",hyperspec$filename),"_"),col.names = Name_group)
  meta <- mutate(meta,Group = apply(dplyr::select(meta,starts_with("group")),1,function(i) paste(i,collapse = "_")))
  hyperspec@data <- cbind(meta,hyperspec@data)
  wavelength_summary <- data.frame(Group = hyperspec$Group[1],
                                   first_wavnumber = hyperspec@wavelength[1],
                                   length = length(hyperspec@wavelength))
  return(wavelength_summary)
}

renishaw <- function(folderpath,Name_group,save = T)
{
  setwd(folderpath)
  filenames <- list.files(pattern = "*.txt")
  ###将所有Txt的第一行删除
  for (filename in filenames)
  {temp <- fread(filename,header = T,sep = "\t")
  fwrite(temp,file = filename,sep = "\t") }
  ###预读取所有数据并根据Wavenumber的统计结果分组
  wavelengths <- read.txt.Renishaw(filenames[1],data = "xyspc")
  wavelength_summarize <- wavelength_summary(wavelengths,Name_group)
  ####在这之前需要去除坐标都为0的情况
  for (filename in filenames[-1])
  {
    temp <- read.txt.Renishaw(filename,data = "xyspc")
    temp_summarize <- wavelength_summary(temp,Name_group)
    wavelength_summarize <- rbind(wavelength_summarize,temp_summarize)
  }
  wavelength_summarize$group <- paste(wavelength_summarize$first_wavnumber,wavelength_summarize$length)
  if (save == T)
  {
    ###根据分组结果读取数据
    for (i in levels(factor(wavelength_summarize$group)))
    {
      group_i <- which(wavelength_summarize$group == i)
      filenames_i <- filenames[group_i]
      final_dataframe <- read.txt.Renishaw(filenames_i[1],data = "xyspc")
      for (filename in filenames_i)
      {
        temp <- read.txt.Renishaw(filename,data = "xyspc")
        final_dataframe <- rbind(final_dataframe,temp)
      }
      data_hyperSpec <- final_dataframe
      data_hyperSpec$Number <- c(1:nrow(data_hyperSpec))
      meta <- as.data.frame(tstrsplit(gsub(".txt","",data_hyperSpec$filename),"_"),col.names = Name_group)
      meta <- mutate(meta,Group = apply(dplyr::select(meta,starts_with("group")),1,function(i) paste(i,collapse = "_")))
      data_hyperSpec@data <- cbind(meta,data_hyperSpec@data)
      single_cell("all_single_cell",data_hyperSpec)
    }
  }
  return(wavelength_summarize)
}



###tSNE
tSNE <- function(data,label,dims = 2,perplexity = 5,max_iter = 10000,outpath)
{
  tsne_data <- data
  tsne <- Rtsne(tsne_data, dims = dims, perplexity=perplexity,theta=0.5, verbose=TRUE, max_iter = max_iter)
  tsne_plot <- as.data.frame(tsne$Y)
  if (dims == 2)
  { colnames(tsne_plot) <- c("tSNE1","tSNE2") }
  tsne_plot$label <- label
  return(tsne_plot)
}
###Save_hyperspec
Save_hyperspec <- function(hyperspec,filename,path)
{
  hyperspec_dataframe <- cbind(select(as.data.frame(hyperspec),-spc,-.row),as.data.frame(hyperspec$spc))
  fwrite(hyperspec_dataframe,paste(path,filename,sep = "\\"),row.names = F,col.names = T,quote = F,sep = "\t")
}

### S-G smooth
S_G_smooth <- function(data,m=0,p = 3,w = 11)
{
  final_data <- data.frame(savitzkyGolay(data,m = m,p = p,
                                         w = w,delta.wav = 1))
  colnames(final_data) <- gsub("X","",colnames(final_data))
  return(final_data)
}
# 2020/08/22
# Rongze Chen
###输出单个光谱数据保存在文件夹####
###根据Group列的分组进行单细胞光谱保存

single_cell <- function(output_folder,hyperSpec)
{
  outpath <- getwd()
  outpath_all <- paste(outpath,"/",output_folder,"/", sep = "")
  factor_Meta <- hyperSpec$Group
  factor_Meta <- factor_Meta[!duplicated(factor_Meta)]
  dir.create(outpath_all)
  for (i in factor_Meta)
  {
    hyperSpec_i <- hyperSpec[hyperSpec$Group == i]
    spec_frame <- as.data.frame(hyperSpec_i$spc)
    meta_data <- data.frame(dplyr::select(as.data.frame(hyperSpec_i), -spc,-.row),metadata = gsub(".txt","",hyperSpec_i$filename))
    outpath_group <- paste(outpath_all,i,"/", sep = "")
    dir.create(outpath_group)
    for (j in (1:nrow(spec_frame)))
    {
      Cells <- data.frame(t(spec_frame[j,]))
      fwrite(Cells,paste(outpath_group,
                         meta_data$metadata[j],"_",meta_data$Number[j],"_",
                         meta_data$x[j],"_",meta_data$y[j],".txt",sep = ""),
             row.names = T,col.names = F,quote = F,sep = "\t")
    }
  }
}
###Spec_shadow
Spec_shadow <- function(hyperspec,Meta_group)
{
  hyperspec_frame <- cbind(as.data.frame(hyperspec)[,1:6],as.data.frame(hyperspec$spc))
  hyperspec_melt <- spc_melt(hyperspec_frame,Meta_group,450,3000)
  hyperspec_melt_summary <- summarySE(hyperspec_melt, measurevar = "value", groupvars = c(Meta_group,"wavenumber"))
  #hyperspec_melt_summary$group_C <- factor(hyperspec_melt_summary$group_C,levels = c("0h","2h","4h","8h","12h"))
  #hyperspec_melt_summary <- filter(hyperspec_melt_summary,group_B != "bg")

  windowsFonts(myFont = windowsFont("Times New Roman"))
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
###SNR####
SNR_function <- function(spect_frame)
{
  Noise_spect <- as.data.frame(spect_frame[,,c(1730~1958,2350~2650)]$spc)
  signal_spect <- as.data.frame(spect_frame[,,c(995~1010)]$spc)
  IBG <- apply(Noise_spect,1,mean)
  IRAM <- apply(signal_spect,1,max)-IBG
  SNR <- IRAM/sqrt(IBG+IRAM)
  return(SNR)
}
###function_all_source

function_all_source <- function(testpath = "G:\\Rcodes\\Function")
{
  setwd(testpath)
  filefolder <- list.files()
  for (i in filefolder)
  {
    workingpath <- paste(testpath,"/",i,sep = "")
    setwd(workingpath)
    filenames <- list.files(pattern = "*.R")
    for (j in filenames) {source(j,encoding = "UTF-8")}
  }
}
###Draw
spc_melt <- function(data,Meta_group,from,to)
{
  ###数据调整成适于画光谱的格式
  cluster.meansd.dataframe_melt <- melt(data,id.vars = Meta_group,variable.name = "wavenumber",value.name = "value")
  cluster.meansd.dataframe_melt$wavenumber <- as.numeric(as.character(cluster.meansd.dataframe_melt$wavenumber))
  cluster.meansd.dataframe_melt$value <- as.numeric(cluster.meansd.dataframe_melt$value)
  cluster.meansd.dataframe_melt <- filter(cluster.meansd.dataframe_melt,wavenumber < to & wavenumber > from)
  return(cluster.meansd.dataframe_melt)
}
###spc_quanlity_control
spc_quanlity_control <- function(hyperspec,abs_noise = 0.1, sd_noise = 0.03,save = T)
{
  hyperspec_1 <- hyperspec[,,c(600~1800)] -
    spc.fit.poly.below(hyperspec[,,c(600~1800)], hyperspec[,,c(600~1800)], poly.order = 1)#3151 Horiba
  hyperspec_2 <- data_hyperSpec[,,c(1800~3000)] -
    spc.fit.poly(data_hyperSpec[,,c(1800~2783,3000)], data_hyperSpec[,,c(1800~3000)], poly.order = 1)#3151 Horiba
  hyperspec_baseline <- cbind(hyperspec_1,hyperspec_2)
  factors <- 1/apply(hyperspec_baseline[,,2500~3000],1,max) #以CH的均值 来归一化
  hyperspec_baseline_normalization <- sweep(hyperspec_baseline,1,factors,"*")
  data <- hyperspec[apply(abs(hyperspec_baseline_normalization[,,1730~1958]),1,mean) < abs_noise & apply(abs(hyperspec_baseline_normalization[,,2350~2650]),1,mean) < abs_noise & apply(hyperspec_baseline_normalization[,,1730~1958],1,sd) < sd_noise & apply(hyperspec_baseline_normalization[,,2350~2650],1,sd) < sd_noise & apply(hyperspec_baseline_normalization[,,500~1730],1,max) < 0.6 & apply(hyperspec_baseline_normalization[,,1730~2830],1,max) < 0.1]
  
  data_min <- apply(data$spc,1,min)
  data <- data[data_min >= -0.5]
  
  
  if (save == T)
  {
    ###根据分组结果读取数据
    data$Number <- c(1:nrow(data))
    meta <- data@data
    meta <- mutate(meta,Group = apply(dplyr::select(meta,starts_with("group")),1,function(i) paste(i,collapse = "_")))
    data@data <- cbind(Group = meta$Group,data@data)
    single_cell("all_single_cell_spc_quanlity_control",data)
  }
  return(data)
}

###获取特定范围的光谱——dataframe
spc_waveumber <- function(spc,from,to)
{
  wave_nums <- as.numeric(as.character(colnames(spc)))
  start <- which.min(abs(wave_nums-from))
  end <- which.min(abs(wave_nums-to))
  final_spc <- spc[,start:end]
  return(final_spc)
}

###
Spec_all <- function(hyperspec,Meta_group)
{
  hyperspec_frame <- cbind(as.data.frame(hyperspec)[,1:5],as.data.frame(hyperspec$spc))
  hyperspec_frame$Number <- c(1:nrow(hyperspec_frame))
  hyperspec_melt <- spc_melt(hyperspec_frame,c(Meta_group,"Number"),450,3200)
  
  windowsFonts(myFont = windowsFont("Times New Roman"))
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
  
  windowsFonts(myFont = windowsFont("Times New Roman"))
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
  
  windowsFonts(myFont = windowsFont("Times New Roman"))
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
pca_plot <- function(data.pca)
{
  windowsFonts(myFont = windowsFont("Times New Roman"))
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
  windowsFonts(myFont = windowsFont("Times New Roman"))
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

####实际测量使用的方法
Sampling_depth2 <- function(CD_ratio)
{
  Time_sample <- 1000
  Pvalue <- Mean_sample <- c()
  if (length(CD_ratio) < 1000)
  { Number_sample <- c(1:(length(CD_ratio)-20)) }
  else { Number_sample <- c(1:1000) }
  
  for (n in Number_sample )#length(CD_ratio_groups$CD_ratio)
  {
    
    Sample_x <- sample(CD_ratio,(20+n))
    last_Mean_sample_x <- mean(Sample_x)
    for (j in 1:Time_sample)
    { 
      Sample_y <- c(Sample_x,sample(CD_ratio,10))
      factor_a <- mean(Sample_y)
      Mean_sample <- c(Mean_sample,factor_a)
    }
    Pvalue <- c(Pvalue,sum(abs(1-Mean_sample/last_Mean_sample_x) < 0.01)/Time_sample)
    Mean_sample <- c()
  }
  return(Pvalue)
}





