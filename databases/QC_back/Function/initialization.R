
initialization <- function(working_dir,output_folder = "output")
{
  outpath <- paste(working_dir,"\\",output_folder,sep = "")
  setwd(working_dir)
  dir.create(outpath)
  print(paste("当前的工作路径是: ",working_dir,sep = ""))
  print(paste("当前的输出路径是: ",outpath,sep = ""))
}

