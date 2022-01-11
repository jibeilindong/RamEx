Sampling_depth2 <- function(CD_ratio_group)
{
  Time_sample <- 1000
  Pvalue <- Mean_sample <- c()
  if (length(CD_ratio_group) < 300)
  { Number_sample <- c(1:(length(CD_ratio_group) - 20)) }
  else {Number_sample <- c(1:280) }
  
  for (n in Number_sample )#length(CD_ratio_groups$CD_ratio)
  {
    
    Sample_x <- sample(CD_ratio_group,(20 + n))
    last_Mean_sample_x <- mean(Sample_x)
    for (j in 1:Time_sample)
    { 
      Sample_y <- c(Sample_x,sample(CD_ratio_group,10))
      factor_a <- mean(Sample_y)
      Mean_sample <- c(Mean_sample,factor_a)
    }
    Pvalue <- c(Pvalue,sum(abs(1-Mean_sample/last_Mean_sample_x) < 0.01)/Time_sample)
    Mean_sample <- c()
  }
  return(Pvalue)
}
