
###PCA
pca <- function(data,label,name = "",outpath)
{
  data_pca <- prcomp(data,cor=T,scores=T,scale. = TRUE)
  pca_data_plot <- data_pca$x
  pca_data_plot <- data.frame(pca_data_plot)
  pca_data_plot <- cbind(label,pca_data_plot)
  return(pca_data_plot)
}