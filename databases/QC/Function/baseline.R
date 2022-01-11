
###baseline
baseline <- function(hyperspec,poly.order = 3)
{
  data_baseline <- hyperspec-
    spc.fit.poly.below (hyperspec, hyperspec, poly.order = poly.order)
  data_baseline <- data_baseline[!duplicated(data_baseline@data[,ncol(data_baseline)])]
  return(data_baseline)
}
