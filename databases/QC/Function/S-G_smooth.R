
### S-G smooth
S_G_smooth <- function(data,m=0,p = 3,w = 11)
{
  final_data <- data.frame(savitzkyGolay(data[,-1],m = m,p = p,
                                         w = w,delta.wav = 1))
  colnames(final_data) <- gsub("X","",colnames(final_data))
  final_data <- cbind(filename = data[,1],final_data)
  return(final_data)
}
