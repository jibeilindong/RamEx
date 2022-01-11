hyperspec_to_dataframe <- function(hyperspec,from,to)
{
  data_frame <- cbind(as.data.frame(hyperspec)[,from:to],as.data.frame(hyperspec$spc))
  return(data_frame)
}
