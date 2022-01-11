
###Normalization
Normalization <- function(my_hyperSpec,from = 2500,to = 3000,nor_cal = "max")
{
  my_hyperspec_spc <- my_hyperSpec$spc
  if (nor_cal == "max")
  {   my_hyperspec_spc_factor <- apply(my_hyperspec_spc,1,max)}  # 以CH的均值 来归一化 
  else if (nor_cal == "sum")
  {  my_hyperspec_spc_factor <- apply(my_hyperspec_spc,1,sum)}
  else
  { print("Error: nor_cal should be max or sum") }
  data_nor <- my_hyperspec_spc/my_hyperspec_spc_factor
  my_hyperSpec$spc <- data_nor
  my_hyperSpec <- my_hyperSpec[!duplicated(my_hyperSpec@data[,ncol(my_hyperSpec)])]
  return(my_hyperSpec)
}
