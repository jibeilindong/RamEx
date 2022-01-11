
install_packages_for_me <- function(used_packages = c("MASS","lattice","RColorBrewer","dplyr","mixOmics","optparse","pROC","R.utils","ggplot2",
                                                      "reshape2","tidyverse","hyperSpec","ggpubr","data.table","scales","Rtsne","pls","e1071",
                                                      "ggbiplot","ggtext","prospectr","devtools","optparse2"))
{
  already_installed_packages <- library()$results[,1]
  install_pack <- setdiff(used_packages, already_installed_packages)
  for ( i in install_pack)
    {
      install.packages(i)
    }

  for ( i in used_packages)
    {
      library(i,character.only = T)
    }
  install_github("vqv/ggbiplot")
  install_github("mixOmicsTeam/mixOmics")
  library(ggbiplot)
  library(mixOmics)
  library(optparse)
}

