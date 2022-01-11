RTCS_Pvalue <- function(alldata_N,Original_value = 20)
{
  alldata_N  <-as.data.table(alldata_N)
  factor_alldata <- select(alldata_N,Group,N)
  factor_alldata <- factor_alldata[!duplicated(factor_alldata),]
  setkey(alldata_N,Group,N)
  
  Pvaluedata_n <- c()
  for (i in 1:nrow(factor_alldata))
  {
    Pvalue <- as.numeric(Sampling_depth2(alldata_N[factor_alldata[i,]]$CD_ratio))
    if (length(Pvalue) == 280)
    {
      Number_sample <- 1:280 
    }else
    {
      Number_sample = c(1:(length(alldata_N[factor_alldata[i,]]$CD_ratio)-Original_value))
    }
    Pvaluedata <- data.frame(Pvalue = Pvalue,Group = factor_alldata[i,1],N = factor_alldata[i,2],Number_sample,Number_sample_max_Pvalue = Number_sample[max(which(Pvalue < 0.95))])
    Pvaluedata_n <- rbind(Pvaluedata_n,Pvaluedata)
  }
  return(Pvaluedata_n)
}
