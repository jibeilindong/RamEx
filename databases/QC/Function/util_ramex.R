#################################################################
# Function:  RamEX util  
# Last update: 2021-11-02, Jing Gongchao;
#################################################################
# install necessary libraries
p <- c("optparse","MASS","lattice", "RColorBrewer", "dplyr", "pROC", "R.utils", "ggplot2", "reshape2", "tidyverse", "hyperSpec", "ggpubr","data.table","scales","Rtsne","pls","e1071","prospectr","Rmisc","ggtext")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://cloud.r-project.org/")
    suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
