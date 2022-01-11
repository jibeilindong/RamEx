
###baseline
baseline <- function(data,poly.order = 3)
{
  wavelength <- as.numeric(colnames(data))
  data_hyperSpec <- new ("hyperSpec", spc = data, wavelength = wavelength)
  data_baseline <- data_hyperSpec-
    spc.fit.poly.below (data_hyperSpec, data_hyperSpec, poly.order = poly.order)
  data_baseline.data.frame <- as.data.frame(data_baseline$spc)
  plotspc(data_baseline)#,spc.nmax = nrow(data_baseline)
  return(data_baseline.data.frame)
}