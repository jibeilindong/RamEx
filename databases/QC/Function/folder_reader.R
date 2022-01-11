##批量读取文件夹数据####
folder_reader <- function(folderpath)
{
  meta_data <- c()
  filefolder <- list.files(folderpath)
  workingpath <- paste(folderpath,"/",filefolder[1],sep = "")
  filenames <- list.files(workingpath,pattern = "*.txt")
  final_data <- spectrum <- as.data.frame(fread(paste(workingpath,"/",filenames[1],sep = ""), header = F, sep = "\t")[,1])
  for (foldername in filefolder)
  {
    workingpath <- paste(folderpath,"/",foldername,sep = "")
    filenames <- list.files(workingpath,pattern = "*.txt")
    everyfolder_data <- spectrum
    meta_data_n <- filenames
    meta_data <-  c(meta_data,meta_data_n)
    for (filename in filenames)
    {
      data <- as.data.frame(fread(paste(workingpath,"/",filename,sep = ""),header = F, sep = "\t"))
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
