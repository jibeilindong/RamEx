###获取特定范围的光谱——dataframe
spc_waveumber <- function(spc,from,to)
{
  wave_nums <- as.numeric(as.character(colnames(spc)))
  start <- which.min(abs(wave_nums-from))
  end <- which.min(abs(wave_nums-to))
  final_spc <- spc[,start:end]
  return(final_spc)
}


###获取特定位置的光谱——dataframe
spc_at_wavenumber <- function(spc,at)
{
  wave_nums <- as.numeric(as.character(colnames(spc)))
  if (length(at) == 1)
  {
    start <- which.min(abs(wave_nums-at))
  }else
  {
    for (i in 1:length(at))
    {
      start_i <- which.min(abs(wave_nums-at[i]))
      if (i == 1)
      {
        start <- start_i
      }else
      {
        start <- c(start,start_i)
      }
    }
  }
  final_spc <- spc[,start]
  return(final_spc)
}
