new_hyperSpec <- function(spc,Name_group)
{
  filename <- spc[,1]
  meta_data <-  as.data.frame(tstrsplit(gsub(".txt","",spc[,1]),"_"),col.names = Name_group)
  spectrum <- as.numeric(colnames(spc[,-1]))
  data_hyperSpec<-new ("hyperSpec", data = meta_data,spc = spc,wavelength=spectrum)
  return(data_hyperSpec)
}
