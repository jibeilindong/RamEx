
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
      cell <- data.frame(t(cluster.meansd.dataframe[i,-(1:2)]))
      fwrite(cell,
             paste(outpath,"/",filelabel,"Mean_by_group/",cluster.meansd.dataframe$Group.1[i],"_.txt",sep=""),
             row.names=T,col.names=F,quote=F,sep = "\t")
    }
  }
  return(cluster.meansd.dataframe)
}



