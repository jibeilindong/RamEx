get_Meta <- function(all_spc)
{
  filename <- all_spc$filename
  meta_data <-  as.data.frame(tstrsplit(gsub(".txt","",filename),"_"),col.names = Name_group)
  meta_data <- cbind(filename,meta_data)
  return(meta_data)
}
