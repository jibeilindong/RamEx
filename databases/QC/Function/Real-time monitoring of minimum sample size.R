
####实际测量使用的方法
Minimum_Samplesize_Realtime <- function(CD_ratio,iterations = 1000,initial_number = 20, Max_sample_number = 300,increase_number = 5,error = 0.05,seq = 5)
{
  Pvalue <- P_colname <- Mean_sample <- c()
  if (length(CD_ratio) < Max_sample_number)
  { Number_sample <- length(CD_ratio)-initial_number }
  else { Number_sample <- Max_sample_number }
  
  for (n in seq(1,Number_sample,seq))#length(CD_ratio_groups$CD_ratio)
  {
    Sample_x <- sample(CD_ratio,(initial_number+n))
    last_Mean_sample_x <- mean(Sample_x)
    for (j in 1:iterations)
    { 
      Sample_y <- c(Sample_x,sample(CD_ratio,increase_number))
      factor_a <- mean(Sample_y)
      Mean_sample <- c(Mean_sample,factor_a)
    }
    Pvalue <- c(Pvalue,sum(abs(1-Mean_sample/last_Mean_sample_x) < error)/iterations)
    P_colname <- c(P_colname,(initial_number+n))
    Mean_sample <- c()
  }
  Pframe <- data.frame(Sample_number = P_colname,Pvalue = Pvalue)
  return(Pframe)
}
