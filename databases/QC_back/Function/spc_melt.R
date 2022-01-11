###Draw
spc_melt <- function(data,Meta_group,from,to)
{
  ###数据调整成适于画光谱的格式
  cluster.meansd.dataframe_melt <- melt(data,id.vars = Meta_group,variable.name = "wavenumber",value.name = "value")
  cluster.meansd.dataframe_melt$wavenumber <- as.numeric(as.character(cluster.meansd.dataframe_melt$wavenumber))
  cluster.meansd.dataframe_melt$value <- as.numeric(cluster.meansd.dataframe_melt$value)
  cluster.meansd.dataframe_melt <- filter(cluster.meansd.dataframe_melt,wavenumber < to & wavenumber > from)
  return(cluster.meansd.dataframe_melt)
}
