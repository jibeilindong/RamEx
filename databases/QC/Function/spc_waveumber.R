###获取特定范围的光谱——dataframe
spc_waveumber <- function(spc,from,to)
{
  wave_nums <- as.numeric(as.character(colnames(spc)))
  start <- which.min(abs(wave_nums-from))
  end <- which.min(abs(wave_nums-to))
  final_spc <- spc[,start:end]
  return(final_spc)
}