spc_quanlity_control <- function(my_hyperspec_1,abs_noise = 0.1, sd_noise = 0.03)
{
  hyperspec_1 <- my_hyperspec_1[,,c(600~1800)] - 
    spc.fit.poly.below(my_hyperspec_1[,,c(600~1800)], my_hyperspec_1[,,c(600~1800)], poly.order = 1)#3151 Horiba
  hyperspec_2 <- my_hyperspec_1[,,c(1800~3000)] - 
    spc.fit.poly(my_hyperspec_1[,,c(1800~2783,3000)], my_hyperspec_1[,,c(1800~3000)], poly.order = 1)#3151 Horiba
  hyperspec_baseline <- cbind(hyperspec_1,hyperspec_2)
  hyperspec_baseline_normalization <- Normalization(hyperspec_baseline,2500,3000,"max")
  my_hyperspec_spc <- hyperspec_baseline_normalization$spc
  good_vector <- apply(abs(spc_waveumber(my_hyperspec_spc,1730,1958)),1,mean) < abs_noise & 
    apply(abs(spc_waveumber(my_hyperspec_spc,2350,2650)),1,mean) < abs_noise & 
    apply(spc_waveumber(my_hyperspec_spc,1730,1958),1,sd) < sd_noise & 
    apply(spc_waveumber(my_hyperspec_spc,2350,2650),1,sd) < sd_noise
  hyperspec_good <- my_hyperspec_1[good_vector] 
  hyperspec_good_min <- apply(hyperspec_good$spc,1,min)
  hyperspec_good <- hyperspec_good[hyperspec_good_min >= -0.5]
  hyperspec_good$Number <- c(1:nrow(hyperspec_good))
  hyperspec_good <- hyperspec_good[!duplicated(hyperspec_good@data[,ncol(hyperspec_good)])]
  return(hyperspec_good)
}
