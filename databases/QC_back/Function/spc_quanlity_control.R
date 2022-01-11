spc_quanlity_control <- function(hyperspec,abs_noise = 0.1, sd_noise = 0.03,save = F)
{
  hyperspec_1 <- hyperspec[,,c(600~1800)] - 
    spc.fit.poly.below(hyperspec[,,c(600~1800)], hyperspec[,,c(600~1800)], poly.order = 1)#3151 Horiba
  hyperspec_2 <- data_hyperSpec[,,c(1800~3000)] - 
    spc.fit.poly(data_hyperSpec[,,c(1800~2783,3000)], data_hyperSpec[,,c(1800~3000)], poly.order = 1)#3151 Horiba
  hyperspec_baseline <- cbind(hyperspec_1,hyperspec_2)
  factors <- 1/apply(hyperspec_baseline[,,2500~3000],1,max) #以CH的均值 来归一化
  hyperspec_baseline_normalization <- sweep(hyperspec_baseline,1,factors,"*")
  data <- hyperspec[apply(abs(hyperspec_baseline_normalization[,,1730~1958]),1,mean) < abs_noise & apply(abs(hyperspec_baseline_normalization[,,2350~2650]),1,mean) < abs_noise & apply(hyperspec_baseline_normalization[,,1730~1958],1,sd) < sd_noise & apply(hyperspec_baseline_normalization[,,2350~2650],1,sd) < sd_noise & apply(hyperspec_baseline_normalization[,,500~1730],1,max) < 0.6 & apply(hyperspec_baseline_normalization[,,1730~2830],1,max) < 0.1] 
  data_min <- apply(data$spc,1,min)
  data <- data[data_min >= -0.5]
  data$Number <- c(1:nrow(data))
  meta <- data@data
  meta <- mutate(meta,Group = apply(dplyr::select(meta,starts_with("group")),1,function(i) paste(i,collapse = "_")))
  data@data <- cbind(Group = meta$Group,data@data)
  if (save == T)
  {
    ###根据分组结果读取数据
    single_cell("all_single_cell_spc_quanlity_control",data)
  }
  return(data)
}
