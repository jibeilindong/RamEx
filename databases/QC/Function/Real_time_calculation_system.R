RTCS_group <- function(CD_ratio_groups,Original_value = 20,initial_value = 5,by = 5)
{
  alldata  <-as.data.table(CD_ratio_groups)
  factor_alldata <- levels(factor(alldata$Group))
  setkey(alldata,Group)
  
  ###以mean为主的取样深度--实际测量使用的代码###
  alldata_N <- alldata_n <- c()
  for (i in 1:length(factor_alldata))
  {
    all_data_i <- alldata[factor_alldata[i]]
    for (j in seq(initial_value,(nrow(all_data_i)-Original_value),by = by))
    {
      alldata_i <- sample(all_data_i,(j+Original_value))
      alldata_i <- data.frame(alldata_i,N = j)
      alldata_n <- rbind(alldata_n,alldata_i)
    }
    alldata_N <- rbind(alldata_N,alldata_n)
    alldata_n <- c()
  }
  
  return(alldata_N)
}
