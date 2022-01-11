Save_hyperspec <- function(hyperspec,filename,path)
{
  hyperspec_dataframe <- cbind(select(as.data.frame(hyperspec),-spc,-.row),as.data.frame(hyperspec$spc))
  fwrite(hyperspec_dataframe,paste(path,filename,sep = "\\"),row.names = F,col.names = T,quote = F,sep = "\t")
}
