# 2020/08/22 
# Rongze Chen
###输出单个光谱数据保存在文件夹####
###根据Group列的分组进行单细胞光谱保存

single_cell <- function(outpath,output_folder,hyperSpec)
{
  outpath_all <- paste(outpath,"/",output_folder,"/", sep = "")
  factor_Meta <- hyperSpec$Group
  factor_Meta <- factor_Meta[!duplicated(factor_Meta)]
  dir.create(outpath_all)
  for (i in factor_Meta)
  {
    hyperSpec_i <- hyperSpec[hyperSpec$Group == i]
    spec_frame <- as.data.frame(hyperSpec_i$spc)
    meta_data <- data.frame(dplyr::select(as.data.frame(hyperSpec_i), -spc,-.row),metadata = gsub(".txt","",basename(hyperSpec_i$filename)))
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

