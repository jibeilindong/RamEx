

###tSNE
tSNE <- function(data,label,dims = 2,perplexity = 5,max_iter = 10000,outpath)
{
  tsne_data <- data
  tsne <- Rtsne(tsne_data, dims = dims, perplexity=perplexity,theta=0.5, verbose=TRUE, max_iter = max_iter)
  tsne_plot <- as.data.frame(tsne$Y)
  if (dims == 2)
  { colnames(tsne_plot) <- c("tSNE1","tSNE2") }
  tsne_plot$label <- label
  return(tsne_plot)
}