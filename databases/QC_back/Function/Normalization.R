
###Normalization
normalization <- function(hyperSpec,from = 2500,to = 3000,cal = "max")
{
  if (cal == "max")
  { factors <- 1/apply(hyperSpec[,,from~to],1,max)}  # 以CH的均值 来归一化 
  else if (cal == "sum")
  { factors <- 1/apply(hyperSpec,1,sum)  }
  else
  { print("Error: cal should be max or sum") }
  data_normalization <- sweep(hyperSpec,1,factors,"*")
  return(data_normalization)
}
