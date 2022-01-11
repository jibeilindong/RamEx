Mean_spc <- function(data_baseline_normalization,Name_group)
{
  hyperspec_frame <- data_baseline_normalization
  #hyperspec_frame <- cbind(as.data.frame(data_baseline_normalization@data),as.data.frame(data_baseline_normalization_spc))
  hyperspec_melt <- reshape2::melt(hyperspec_frame, id.vars = Name_group, variable.name="wavenumber")
  #hyperspec_melt_summary <- summarySE(hyperspec_melt, measurevar = "value", groupvars = c(Name_group,"wavenumber"))
  hyperspec_melt_summary <- select(hyperspec_melt,c(Name_group,"wavenumber","value"))
  return(hyperspec_melt_summary)
}

RBCS_selectpeak <- function(mean_all,Pvalue = 0.0001)
{
  diff_cta <- sub_mean <- NULL
  for (i in levels(factor(mean_all$Group)))
  { 
    sub_mean_i <- NULL
    sub_mean_j <- NULL
    mean_all_i <- filter(mean_all,Group == i)
    for (j in levels(factor(mean_all_i$wavenumber))) {
      mean_all_i_j <- filter(mean_all_i,wavenumber == j)
      control <- filter(mean_all_i_j,RBCS_group == "control")$value
      test <- filter(mean_all_i_j,RBCS_group == "test")$value
      p <- t.test(control,test,paired = FALSE)$p.value
      if (p > Pvalue) {
        sub_mean_i <- mean(test) - mean(control)
      }else{sub_mean[i] <- 0}
      sub_mean_j <- c(sub_mean_j,sub_mean_i)
    }
    if(length(sub_mean)==0){sub_mean <- sub_mean_j
	
    }
    else{sub_mean <- data.frame(sub_mean, sub_mean_j)}
	
    #sub_mean <- t(sub_mean)
    
    #rownames(sub_mean) <- i
  }
	
  diff_cta <- rbind(diff_cta,sub_mean)
  rownames(diff_cta) <- levels(factor(mean_all$wavenumber))
  colnames(diff_cta) <- levels(factor(mean_all$Group))
  return(diff_cta)
}

Spc_at_wavenumber <- function(spc,at)
{
  wave_nums <- as.numeric(as.character(colnames(spc)))
  #at <- as.numeric(at)
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

