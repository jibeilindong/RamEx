CD_ratio <- function(hyperspec)
{
  meta_data <- select(as.data.frame(hyperspec),-spc,-.row)
  data_baseline_normalization_frame  <- cbind(meta_data,as.data.frame(hyperspec$spc))
  wave_nums <- as.numeric(as.character(colnames(data_baseline_normalization_frame)))
  
  CD_start <- which.min(abs(wave_nums - 2050))
  CD_end <- which.min(abs(wave_nums - 2300))
  CH_start <- which.min(abs(wave_nums - 3050))
  CH_end <- which.min(abs(wave_nums - 2800))
  #CD_ratio<-rowSums(CD)/(rowSums(CD)+rowSums(CH)) #参照
  data_baseline_normalization_frame_CD_ratio <- mutate(data_baseline_normalization_frame,
                                                       CD_ratio = rowSums(data_baseline_normalization_frame[,CD_start:CD_end])/(rowSums(data_baseline_normalization_frame[,CD_start:CD_end]) + rowSums(data_baseline_normalization_frame[,CH_start:CH_end])))
  CD_ratio_groups <- select(data_baseline_normalization_frame_CD_ratio,CD_ratio)  #需要设置组别，在实验中哪几个参数作为主要变量，需要在这里更换
  CD_ratio_groups <- cbind(meta_data,CD_ratio_groups)
  return(CD_ratio_groups)
}
