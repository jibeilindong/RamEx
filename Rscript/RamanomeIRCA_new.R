source(sprintf('%s/Rscript/util.R',"/home/gene/jinggc/RamEX"))
GlobalBandsAnnfile<- sprintf('%s/databases/Global_bands_annotation.txt',"/home/gene/jinggc/RamEX") #GlobalBandsAnnfile<-NULL
LocalBandsAnnfile<- sprintf('%s/databases/Local_bands_annotation.txt',"/home/gene/jinggc/RamEX")
category<-c("Celltype","Timepoint")#c("Species","Strain","Condition","Timepoint","User","Date")
cor_cutoff<-0.6

plot_local_IRCN<-FALSE
plot_global_IRCN<-FALSE
#device="pdf"

#---------------------
# output pathway
#---------------------
outpath <- "./Timepoint_neg0.6/"
dir.create(outpath)
outpath_png <- paste(outpath,"png/",sep="")
dir.create(outpath_png)

#-------------------------------
# Spec range trimming
#-------------------------------
#SpecRange_remove<-1800:2600

options(warn=-1)
#-------------------------------
# Add prefix to outpath
#-------------------------------
#prefix<-paste(paste(ifelse(Pos_Edge, "Pos", "NoPos"),ifelse(Neg_Edge, "Neg", "NoNeg"),sep="-"), cor_cutoff, "", sep="_")
#outpath<-paste(outpath, prefix, sep="")
#-------------------------------
# Spectral data input
#-------------------------------
matrixfile <- "CC124_new.txt"
mat <- read.table(matrixfile, header = T, row.names = 1, sep="\t"); 
#remove 3d_s1-s2outliers
#mat<-mat[-c(603,610,623,637,631),]
#mat<-mat[order(rownames(mat)),]

#-------------------------------
# scale by dividing sum area
#-------------------------------
scale_sum<-function(data){
  sum<-sum(data)
  nor<-data/sum
}
mat<-apply(mat, 1, scale_sum)
mat<-data.frame(t(mat))
#write.table(mat,paste(outpath,"Nor_",matrixfile,sep=""),row.names=T,quote=F,sep="\t")
write.csv(mat,paste(outpath,"Nor_",matrixfile,".csv",sep=""),row.names=T,quote=F)
#mat_show<-mat[1:5,1:5]
#-------------------------------
# Clean data (Negative, NA, Inf, or -Inf)
#-------------------------------
mat<-CleanData(mat)
cat("The number of Raman shifts: ", ncol(mat) , "\n")
#-------------------------------
# Wave number trimming
#-------------------------------
raw_wn<-as.numeric(gsub("[A-Z]", "", colnames(mat)))
wn0<-formatC(raw_wn, digits=0, format="f")
if(!is.null(GlobalBandsAnnfile)){ global_bands_ann<-read.table(GlobalBandsAnnfile,header=T, sep="\t")
Wave_num <- data.frame(Wave_num = wn0)
global_bands_ann <- merge(Wave_num, global_bands_ann, by="Wave_num", all.x=TRUE, sort = TRUE)
global_bands_ann[is.na(global_bands_ann)] <- "Unknown"
global_bands_ann$Wave_num<-as.numeric(global_bands_ann$Wave_num)
global_bands_ann<-global_bands_ann[with(global_bands_ann, order(Wave_num)),]

}else{
  stop("Please provide global Bands annotation file!")
}
if(!is.null(LocalBandsAnnfile)){ local_bands_ann<-read.table(LocalBandsAnnfile,header=T, sep="\t") }
# Local IRCNs by chordDiagram
# Generate local networks by selected bands
#-------------------------------
#Peak_sub<-c("B907","B959","B1005","B1129",
#         "B1156",  "B1186", "B1260","B1489")
#Peak_sub<-paste("B",local_bands_ann$Wave_num,sep = "")
Peak_sub<-as.character(local_bands_ann$Wave_num)
bands<-Peak_sub
#aaa<-cor_mats_rpvalue[["Yeast__day 0"]]$r
#which(bands%in%colnames(aaa))
Peaks_fixed_cor_mats_rpvalue <- lapply(cor_mats_rpvalue, 
                                   function(x) submat <- list(x$r[bands, bands],x$P[bands, bands],x$n[bands, bands]))
names(Peaks_fixed_cor_mats_rpvalue)<-names(cor_mats_rpvalue)
save(Peaks_fixed_cor_mats_rpvalue,
     file=paste(outpath,
     #file=paste("/mnt/data6/heyh/20180831_Analysis_RWAS/Group_neg0.6/",
                "Peaks_fixed_mannulselected_cor_mats_rpvalue.RData",sep=""))

#------------------------------------------------------------------

#------------------------------------------------------------------
if(!is.null(LocalBandsAnnfile) & plot_local_IRCN){ 
  
  #Peak_sub<-paste("B",local_bands_ann$Wave_num,sep = "")
  Peak_sub<-local_bands_ann$Wave_num
  #-------------------------------
  # Generate local networks by selected bands
  #-------------------------------
  bands<-as.character(Peak_sub)
  bands_cor_mats_rpvalue <- lapply(cor_mats_rpvalue, 
                                   function(x) submat <- list(x$r[bands, bands],x$P[bands, bands],x$n[bands, bands]))
  require(Cairo)
  #---------------------------------------------------
  # png plot_local_chorddiagram
  #---------------------------------------------------
  # 1.# Pos_Edge=F, Neg_Edge=T, Threshold=0.6,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_local_neg6_png_rpvalue <-paste(outpath_png,"local_neg6_png_rpvalue/",sep = "")
  dir.create(outpath_local_neg6_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20151101"
    #group_name<-"Chlamydomonas reinhardtii__CC124__N+__072h__HeYH__20160817"
	#group_name<-"CC124__072h"
    #Corr_mat<-bands_cor_mats[[group_name]]
    #a_degree<-data.frame(factors = rep("a",nrow(Corr_mat)),
    #                     x=as.numeric(gsub("[A-Z]", "", row.names(Degree))),
    #                     y=Degree[,group_name]/max(Degree))
    #group_name_degree<-paste("X",group_name,sep = "")
    CairoPNG(file = paste(outpath_local_neg6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""), width=3, height=3, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=F, Neg_Edge=T, Threshold=0.6, local_bands_ann)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 2.# Pos_Edge=T, Neg_Edge=F, Threshold=0.6,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_local_pos6_png_rpvalue <- paste(outpath_png,"local_pos6_png_rpvalue/",sep = "")
  dir.create(outpath_local_pos6_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    CairoPNG(file = paste(outpath_local_pos6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""),width=3, height=3, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=F, Threshold=0.6, local_bands_ann)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 3.# Pos_Edge=T, Neg_Edge=T, Threshold=0.6,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__144h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__024h__WangTT__20120710","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219")
  device="png"
  outpath_local_neg6pos6_png_rpvalue <- paste(outpath_png,"local_neg6pos6_png_rpvalue/",sep = "")
  dir.create(outpath_local_neg6pos6_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    CairoPNG(file = paste(outpath_local_neg6pos6_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""), width=3, height=3, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=T, Threshold=0.6, local_bands_ann)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 4.# Pos_Edge=F, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_local_neg0_png_rpvalue <-  paste(outpath_png,"local_neg0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_neg0_png_rpvalue/"
  dir.create(outpath_local_neg0_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    CairoPNG(file = paste(outpath_local_neg0_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png" , sep=""), width=3, height=3, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=F, Neg_Edge=T, Threshold=0, local_bands_ann)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 5.# Pos_Edge=T, Neg_Edge=F, Threshold=0,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_local_pos0_png_rpvalue <- paste(outpath_png,"local_pos0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_pos0_png_rpvalue/"
  dir.create(outpath_local_pos0_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    CairoPNG(file = paste(outpath_local_pos0_png_rpvalue, "Lobal_chordDiagram_", group_name, sep=""), width=3, height=3, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=F, Threshold=0, local_bands_ann)
    dev.off()  
  }
  #-----------------------------------------------------------------------------------------------
  
  #---------------------------------------------------
  # 6.# Pos_Edge=T, Neg_Edge=T, Threshold=0,
  #---------------------------------------------------
  #names_sub<-c("Chlamydomonas reinhardtii__CC124__N+__000h__HeYH__20160817","Chlamydomonas reinhardtii__CC124__N-__012h__HeYH__20160817",
  #             "Chlamydomonas reinhardtii__CC124__N-__096h__HeYH__20160817","Nannochloropsis oceanica__IMET1__N-__096h__WangTT__20120710",
  #             "Saccharomyces cerevisiae__Y50049__30mM__096h__ZhangP__20161219","Escherichia coli__DH5a__kan__003h__TengL__20160519")
  device="png"
  outpath_local_neg0pos0_png_rpvalue <- paste(outpath_png,"local_neg0pos0_png_rpvalue/",sep = "")#"/mnt/data6/heyh/Plot_Fig8BC_chordDiagram_20200623/local_neg0pos0_png_rpvalue/"
  dir.create(outpath_local_neg0pos0_png_rpvalue)
  for (group_name in names(bands_cor_mats_rpvalue)) {
  #for (group_name in names_sub) {
    CairoPNG(file = paste(outpath_local_neg0pos0_png_rpvalue, "Lobal_chordDiagram_", group_name, ".png", sep=""), width=3, height=3, bg="white")
    Plot_local_chordDiagram_rpvalue_new(group_name, Pos_Edge=T, Neg_Edge=T, Threshold=0, local_bands_ann)
    dev.off()  
  }
  #----------------------------------------------------------------------------------------------
}
#---------------------------------------------------------------------------------------------------------------------


#--------------------------------
# save data
#--------------------------------
#save(cor_mats,v_cor_mats,v_cor_mats_df,v_cor_mats_trimmed,g_array,g,Degree,mean_spectra,dist_cor_mats,
#     file = paste(outpath,"mydata_maindataset_by_", category,".RData",sep=""))
save.image(file = paste(outpath,"mydata_by_", category,".RData",sep="")) #save all data into one .RData file

#-------------------------------
# Local IRCNs by chordDiagram
# Desk top
#-------------------------------
#---------------------
# loading the dataset
#---------------------
#load(paste("/mnt/data6/heyh/IRCN_CC124_20200624/Timepoint_neg0.6/",
#           "mydata_by_Celltype__Timepoint.RData",sep=""))
load(paste("./Timepoint_neg0.6/",
           "Peaks_fixed_mannulselected_cor_mats_rpvalue.RData",sep=""))

bands_cor_mats_rpvalue<-Peaks_fixed_cor_mats_rpvalue

# global_bands_annotation 鎷撳睍(+-3)
#Global_bands_annotation<-read.table(GlobalBandsAnnfile,
#                                    header = T,sep="\t")
Global_bands_annotation<-global_bands_ann
i<-4
#i<-79
while (i<=(nrow(Global_bands_annotation)-3)){
  if(Global_bands_annotation$Group_new[i]!="Unknown"){
    Global_bands_annotation$Group_new[c(i-3,i-2,i-1,i+1,i+2,i+3)]<-rep(as.character(Global_bands_annotation$Group_new[i]),6)
    Global_bands_annotation$Assignment_new[c(i-3,i-2,i-1,i+1,i+2,i+3)]<-rep(as.character(Global_bands_annotation$Assignment_new[i]),6)
    i<-i+4
  }else{
      i<-i+1
    }
}
# Libraries
library(dplyr)
library(tidyverse)
library(circlize)
options(knitr.table.format = "html")
library(viridis)
library(igraph)

library(colormap)
library(hrbrthemes)
library(kableExtra)
library(ggraph)#graphlayouts #tidygraph


################
#1.Neg 0 & p.value<0.05 all peaks
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"
outpath <- "./Timepoint_neg0.6/1.Arcdiagram_Peaks_Neg0/"
dir.create(outpath)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  #Corr_mat[Corr_mat>(-0.6)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(0)]<-NA # p.value<0.05
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  
  #--------------------
  # 2. Make the graph
  # mannul grouped
  #--------------------
  #order_Peaks<-c("B957","B966","B999","B1007",
  #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  #order_Peaks<-c("B1526","B1158","B1007",
  #               "B1504","B1152","B999",
  #               "B1193","B1180","B957","B966")
  order_Peaks<-c(#"B16581441",
                 "B1007","B1602",
                 "B478","B865","B938","B1125",
                 "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
                 "B526",  "B577", 
                 "B612","B718",
                 "B2911")
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  
  #dim(coauth)
  new_order<-match(order_Peaks,coauth$name)
  coauth<-coauth[new_order,]
  
  # grouped
  #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  #                                    header = T,sep="\t")
  Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  
  coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  coauth$grp_mannul<-c(#"DU",
                       rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
                       rep("Lipids",2),"Lipids; Carbohydrates")
  coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(##"DU",
                                                       "Proteins","Starch","TAG", "Memanbrane Lipids",
                                                       "Lipids","Lipids; Carbohydrates"))
  #coauth$grp_mannul<-c(#"DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU","Starch","TAG", "Memanbrane Lipids"))
  #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,
  #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  
  #coauth$grp_mannul<-coauth$Group_new
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  #coauth<-coauth[order(coauth$grp_mannul),]
  
  # keep only this people in edges
  connect <- connect %>%
    filter(from %in% coauth$name) %>%
    filter(to %in% coauth$name)
  
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  # prepare a vector of n color in the viridis scale
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  #mycolor <- sample(mycolor, length(mycolor))
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  #mycolor<-c("#440154ff", "#5cc863ff")
  #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  
  #----------------------
  #no group and assignment
  p21<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  
  #----------------------
  #add assignment
  p22<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
         p22,height = 3,width = 3)
  
  #----------------------
  #add group
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  
  }


################
#2.Neg 0 & p.value<0.05 DU only
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"
outpath <- "./Timepoint_neg0.6/2.Arcdiagram_DU_Neg0/"
dir.create(outpath)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(0)]<-NA # r<0
  #Corr_mat[which(colnames(Corr_mat)!="B16581441"),which(colnames(Corr_mat)!="B16581441")]<-NA # DU only
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  
  #--------------------
  # 2. Make the graph
  # mannul grouped
  #--------------------
  #order_Peaks<-c("B957","B966","B999","B1007",
  #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  #order_Peaks<-c("B1526","B1158","B1007",
  #               "B1504","B1152","B999",
  #               "B1193","B1180","B957","B966")
  order_Peaks<-c(#"B16581441",
                 "B1007","B1602",
                 "B478","B865","B938","B1125",
                 "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
                 "B526",  "B577", 
                 "B612","B718",
                 "B2911")
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  
  #dim(coauth)
  new_order<-match(order_Peaks,coauth$name)
  coauth<-coauth[new_order,]
  
  # grouped
  #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  #                                    header = T,sep="\t")
  Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  
  coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  coauth$grp_mannul<-c(#"DU",
    rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
                       rep("Lipids",2),"Lipids; Carbohydrates")
  coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU",
                                                       "Proteins","Starch","TAG", "Memanbrane Lipids",
                                                       "Lipids","Lipids; Carbohydrates"))
  #coauth$grp_mannul<-c("DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("DU","Starch","TAG", "Memanbrane Lipids"))
  #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,
  #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  
  #coauth$grp_mannul<-coauth$Group_new
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  #coauth<-coauth[order(coauth$grp_mannul),]
  
  # keep only this people in edges
  connect <- connect %>%
    filter(from %in% coauth$name) %>%
    filter(to %in% coauth$name)
  
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  # prepare a vector of n color in the viridis scale
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  #mycolor <- sample(mycolor, length(mycolor))
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  #mycolor<-c("#440154ff", "#5cc863ff")
  #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  
  #----------------------
  #no group and assignment
  p21<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  
  #----------------------
  #add assignment
  p22<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
         p22,height = 3,width = 3)
  
  #----------------------
  #add group
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  
}

################
#3.Neg 6 & p.value<0.05 all peaks
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"
outpath <- "./Timepoint_neg0.6/3.Arcdiagram_Peaks_Neg6/"
dir.create(outpath)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(-0.6)]<-NA # p.value<0.05
  #Corr_mat[Corr_mat>(0)]<-NA # p.value<0.05
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  
  #--------------------
  # 2. Make the graph
  # mannul grouped
  #--------------------
  #order_Peaks<-c("B957","B966","B999","B1007",
  #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  #order_Peaks<-c("B1526","B1158","B1007",
  #               "B1504","B1152","B999",
  #               "B1193","B1180","B957","B966")
  order_Peaks<-c(#"B16581441",
                 "B1007","B1602",
                 "B478","B865","B938","B1125",
                 "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
                 "B526",  "B577", 
                 "B612","B718",
                 "B2911")
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  
  #dim(coauth)
  new_order<-match(order_Peaks,coauth$name)
  coauth<-coauth[new_order,]
  
  # grouped
  #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  #                                    header = T,sep="\t")
  Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  
  coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  coauth$grp_mannul<-c(#"DU",
                       rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
                       rep("Lipids",2),"Lipids; Carbohydrates")
  coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU",
                                                       "Proteins","Starch","TAG", "Memanbrane Lipids",
                                                       "Lipids","Lipids; Carbohydrates"))
  #coauth$grp_mannul<-c("DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("DU","Starch","TAG", "Memanbrane Lipids"))
  #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,
  #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  
  #coauth$grp_mannul<-coauth$Group_new
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  #coauth<-coauth[order(coauth$grp_mannul),]
  
  # keep only this people in edges
  connect <- connect %>%
    filter(from %in% coauth$name) %>%
    filter(to %in% coauth$name)
  
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  # prepare a vector of n color in the viridis scale
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  #mycolor <- sample(mycolor, length(mycolor))
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  #mycolor<-c("#440154ff", "#5cc863ff")
  #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  
  #----------------------
  #no group and assignment
  p21<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  
  #----------------------
  #add assignment
  p22<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
         p22,height = 3,width = 3)
  
  #----------------------
  #add group
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  
}



################
#4.Neg 0.6 & p.value<0.05 DU only
################
#---------------------
# output pathway
#---------------------
#outpath <- "E:/RWAS/IRCN_CC124_20200624/M0628Fig7-CC124alltimepoints-fixedbands/"
outpath <- "./Timepoint_neg0.6/4.Arcdiagram_DU_Neg6/"
dir.create(outpath)

for (x in names(bands_cor_mats_rpvalue)){
  #--------------------
  # 0.loading dataset
  #--------------------
  #x<-"00h"
  #x<-"CC124__144h"
  Corr_mat<-bands_cor_mats_rpvalue[[x]][[1]] # r
  Corr_mat_pvalue<-bands_cor_mats_rpvalue[[x]][[2]] # p.value
  Corr_mat[which(seq(1:length(Corr_mat))%in%which(Corr_mat_pvalue<0.05)==FALSE)]<-NA # p.value<0.05
  Corr_mat[Corr_mat>(-0.6)]<-NA # r<0
  #Corr_mat[which(colnames(Corr_mat)!="B16581441"),which(colnames(Corr_mat)!="B16581441")]<-NA # DU only
  dataUU<-data.frame(from=row.names(Corr_mat),Corr_mat)
  
  #--------------------
  # 1. Make the graph
  # no grouped
  #--------------------
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  #dim(coauth)
  new_order<-match(row.names(dataUU),coauth$name)
  coauth<-coauth[new_order,]
  
  # Create a graph object with igraph
  require(igraph)
  mygraph1 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  p1<-ggraph(mygraph1, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    #geom_node_point(aes(size=n, color=as.factor(grp), fill=grp), alpha=0.5) +
    geom_node_point(aes(size=(n+1)),colour="grey20", alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    #scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"1.local-IRCN_arcdiagram_ungrouped_",x,".pdf",sep = ""),p1,height = 1.5,width = 2)
  
  
  #--------------------
  # 2. Make the graph
  # mannul grouped
  #--------------------
  #order_Peaks<-c("B957","B966","B999","B1007",
  #               "B1152",  "B1158", "B1180","B1193","B1504","B1526")
  #order_Peaks<-c("B1526","B1158","B1007",
  #               "B1504","B1152","B999",
  #               "B1193","B1180","B957","B966")
  order_Peaks<-c(#"B16581441",
                 "B1007","B1602",
                 "B478","B865","B938","B1125",
                 "B968", "B1263","B1302","B1443","B1656","B1746","B2855","B3011",
                 "B526",  "B577", 
                 "B612","B718",
                 "B2911")
  # Transform the adjacency matrix in a long format
  connect <- dataUU %>% 
    gather(key="to", value="value", -1) %>%
    mutate(to = gsub("\\.", " ",to)) %>%
    na.omit() 
  
  # Number of connection per person
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> coauth
  
  colnames(coauth) <- c("name", "n")
  bands_more<-row.names(dataUU)[!row.names(dataUU)%in%coauth$name]
  bands_more<-data.frame(name=bands_more,n=rep(0,length(bands_more)))
  coauth<-rbind(coauth,bands_more)
  
  #dim(coauth)
  new_order<-match(order_Peaks,coauth$name)
  coauth<-coauth[new_order,]
  
  # grouped
  #Global_bands_annotation<-read.table(paste("E:/RWAS/IRCN_CC124_20200624/","Global_bands_annotation_forCC1247d.txt",sep = ""),
  #                                    header = T,sep="\t")
  Global_bands_annotation_tem<-Global_bands_annotation[which(paste("B",Global_bands_annotation$Wave_num,sep="")%in%coauth$name),]
  new_order<-match(order_Peaks,paste("B",Global_bands_annotation_tem$Wave_num,sep=""))
  Global_bands_annotation_tem<-Global_bands_annotation_tem[new_order,]
  
  coauth<-cbind(coauth,Global_bands_annotation_tem[,4:5])
  coauth$grp_mannul<-c(#"DU",
                       rep("Proteins",2),rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",2),
                       rep("Lipids",2),"Lipids; Carbohydrates")
  coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c(#"DU",
                                                       "Proteins","Starch","TAG", "Memanbrane Lipids",
                                                       "Lipids","Lipids; Carbohydrates"))
  #coauth$grp_mannul<-c("DU",rep("Starch",4),rep("TAG",8),rep("Memanbrane Lipids",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("DU","Starch","TAG", "Memanbrane Lipids"))
  #coauth$grp_mannul<-c(rep("C13-Carotenoids",3),rep("C12-Carotenoids",3),rep("Unknown",4))
  #coauth$grp_mannul<-factor(coauth$grp_mannul,
  #                          levels=c("C13-Carotenoids","C12-Carotenoids", "Unknown"))
  
  #coauth$grp_mannul<-coauth$Group_new
  #coauth$grp_mannul<-factor(coauth$grp_mannul,levels=c("carbohydrates","lipids" ))
  #coauth<-coauth[order(coauth$grp_mannul),]
  
  # keep only this people in edges
  connect <- connect %>%
    filter(from %in% coauth$name) %>%
    filter(to %in% coauth$name)
  
  # Create a graph object with igraph
  mygraph2 <- graph_from_data_frame( connect, vertices = coauth, directed = FALSE )
  
  # prepare a vector of n color in the viridis scale
  mycolor <- colormap(colormap=colormaps$viridis, nshades=nlevels(coauth$grp_mannul))
  #mycolor <- sample(mycolor, length(mycolor))
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff","grey20")
  #mycolor<-c("#440154ff", "#5cc863ff")
  #mycolor<-c("#440154ff", "#5cc863ff", "grey20")
  #mycolor<-c("#3b518bff","#440154ff", "#5cc863ff", "#21908dff",)
  
  #----------------------
  #no group and assignment
  p21<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_node_point(aes(size=n, color=as.factor(grp_mannul), fill=grp_mannul), alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=name), angle=65, hjust=1.5, nudge_y = -1.1, size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0.4,0), "null"),
      panel.spacing=unit(c(0,0,3.4,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-5.6, 1.2)) 
  ggsave(paste(outpath,"2.1 local-IRCN_arcdiagram_grouped_",x,".pdf",sep = ""),p21,height = 1.5,width = 2)
  
  #----------------------
  #add assignment
  p22<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", Assignment_new,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.2 local-IRCN_arcdiagram_grouped_assignment_",x,".pdf",sep = ""),
         p22,height = 3,width = 3)
  
  #----------------------
  #add group
  p23<-ggraph(mygraph2, layout="linear") + 
    #geom_edge_arc(edge_colour="black", edge_alpha=0.2, edge_width=0.3, fold=TRUE) +
    geom_edge_arc(edge_colour="#66CCFF", edge_alpha=0.2, edge_width=0.3, fold=TRUE,check_overlap = T) +
    geom_node_point(aes(size=n),color=mycolor[unclass=coauth$grp_mannul], fill=mycolor[unclass=coauth$grp_mannul],alpha=0.5) +
    #geom_node_point(aes(size=n+1), alpha=0.5) +
    scale_size_continuous(range=c(0.5,3.5)) +
    scale_color_manual(values=mycolor) +
    #geom_node_text(aes(label=name), angle=65, hjust=1, nudge_y = -1.1, size=2.3) +
    geom_node_text(aes(label=paste(name," [", grp_mannul,"]",sep="")), 
                   colour=mycolor[unclass=coauth$grp_mannul],
                   angle=65, 
                   hjust=1.1, 
                   nudge_y = -1.1, 
                   size=1.5) +
    theme_void() +
    theme(
      legend.position="none",
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null")
    ) +
    expand_limits(x = c(-1.2, 1.2), y = c(-30, 10)) 
  ggsave(paste(outpath,"2.3 local-IRCN_arcdiagram_grouped_group_",x,".pdf",sep = ""),
         p23,height = 4,width = 2)
  
  
}
#----------finished-------------
