
install_packages_for_me <- function(used_packages = c("MASS","lattice","RColorBrewer","dplyr","optparse","pROC","R.utils","ggplot2",
                                                      "reshape2","tidyverse","hyperSpec","ggpubr","data.table","scales","Rtsne","pls","e1071",
                                                      "ggtext","prospectr","devtools","Rmisc"))
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
  library(ggbiplot)
  install_github("mixOmicsTeam/mixOmics")
  library(mixOmics)
}

