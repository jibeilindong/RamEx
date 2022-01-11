new_hyperSpec <- function(spc,Name_group)
{
  filename <- spc[,1]
  meta_data <-  as.data.frame(tstrsplit(gsub(".txt","",spc[,1]),"_"),col.names = Name_group)
  meta_data <- mutate(meta_data,Group = apply(dplyr::select(meta_data,starts_with("group")),1,function(i) paste(i,collapse = "_")))
  meta_data <- cbind(filename,meta_data)
  spectrum <- as.numeric(colnames(spc[,-1]))
  data_hyperSpec<-new ("hyperSpec", data = meta_data,spc = spc[,-1],wavelength=spectrum)
  return(data_hyperSpec)
}
