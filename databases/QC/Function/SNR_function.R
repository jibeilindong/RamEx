###SNR####
SNR_function <- function(spect_frame)
{
  Noise_spect <- as.data.frame(spect_frame[,,c(1730~1958,2350~2650)]$spc)
  signal_spect <- as.data.frame(spect_frame[,,c(995~1010)]$spc)
  IBG <- apply(Noise_spect,1,mean)
  IRAM <- apply(signal_spect,1,max)-IBG
  SNR <- IRAM/sqrt(IBG+IRAM)
  return(SNR)
}
