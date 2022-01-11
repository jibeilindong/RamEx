###删除异常光谱
spc_quanlity_control <- function(data,metadata,draw = FALSE,level = "hard")
{
  wavelength <- as.numeric(colnames(data))
  data_hyperSpec <- new ("hyperSpec", data= metadata,spc = data, wavelength = wavelength)
  data_baseline <- data_hyperSpec[,,c(500~3200)]-
    spc.fit.poly.below (data_hyperSpec[,,c(500~3200)], data_hyperSpec[,,c(500~3200)], poly.order = 3)
  factors <- 1/apply(data_baseline[,,c(2925~2935)],1,max)
  data_normalization <- sweep(data_baseline,1,factors,"*")
  if (level == "hard")
  {
    data_hyperSpec_good <- data_hyperSpec[apply(abs(data_normalization[,,1730~1958]),1,mean) < 0.05
                                          & apply(abs(data_normalization[,,2350~2650]),1,mean) < 0.05
                                          & apply(abs(data_normalization[,,1730~1958]),1,max) < 0.05
                                          & apply(abs(data_normalization[,,2350~2650]),1,max) < 0.05
                                          & apply(data_normalization[,,1730~1958],1,sd) < 0.03
                                          & apply(data_normalization[,,2350~2650],1,sd) < 0.03
                                          & SNR_function(data_normalization) > 0.1 ]
  }
  if (level == "easy")
  {
    data_hyperSpec_good <- data_hyperSpec[apply(abs(data_normalization[,,1730~1958]),1,mean) < 0.3
                                          & apply(abs(data_normalization[,,2350~2650]),1,mean) < 0.3
                                          & apply(abs(data_normalization[,,1730~1958]),1,max) < 0.3
                                          & apply(abs(data_normalization[,,2350~2650]),1,max) < 0.3
                                          & apply(data_normalization[,,1730~1958],1,sd) < 0.3
                                          & apply(data_normalization[,,2350~2650],1,sd) < 0.3
                                          & SNR_function(data_normalization) > 0.1 ]
  }
  data_hyperSpec_good.data.frame <-  data.frame(select(as.data.frame(data_hyperSpec_good),-spc,-.row),data_hyperSpec_good$spc)
  data_hyperSpec_good.data.frame <- na.omit(data_hyperSpec_good.data.frame)
  colnames(data_hyperSpec_good.data.frame) <- gsub("X","",colnames(data_hyperSpec_good.data.frame))
  Meta_group <- colnames(metadata)
  if (draw == TRUE)
  {
    draw_spc(data_hyperSpec_good.data.frame,Meta_group,name="good",weight = 12,height = 4,outpath) 
  }
  plot(data_hyperSpec_good,"spcprctl5")
  return(data_hyperSpec_good.data.frame)
}
