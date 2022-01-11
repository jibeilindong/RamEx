
function_all_source <- function(testpath = "G:\\Rcodes\\Function")
{
  setwd(testpath)
  filefolder <- list.files()
  for (i in filefolder)
  {
    workingpath <- paste(testpath,"/",i,sep = "")
    setwd(workingpath)
    filenames <- list.files(pattern = "*.R")
    for (j in filenames) {source(j,encoding = "UTF-8")}
  }
}
