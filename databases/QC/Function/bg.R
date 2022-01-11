
###bg
bg <- function(data,metadata,K = FALSE)
{
  library("hyperSpec")
  wavelength <- as.numeric(colnames(data))
  data_hyperSpec <- new ("hyperSpec", data= metadata,spc = data, wavelength = wavelength)
  bg_hyperSpec <- data_hyperSpec[ data_hyperSpec$C== "bg" ]
  bgsub_hyperSpec <- data_hyperSpec[ data_hyperSpec$C!= "bg" ]
  
  for(i in levels(factor(metadata$A)))
  {
    for ( j in levels(factor(metadata$C)))
    {
      if (j == 1 )
      {
        bg_hyperSpec_ij <- apply(bg_hyperSpec[bg_hyperSpec$A == i], 2, mean)
        bgsub_hyperSpec_ij <- bgsub_hyperSpec[bgsub_hyperSpec$A == i&bgsub_hyperSpec$C == j] 
        for (N in bgsub_hyperSpec_ij$Number)
        {
          bgsub_hyperSpec_ijN <- bgsub_hyperSpec[bgsub_hyperSpec$A == i&bgsub_hyperSpec$C == j&bgsub_hyperSpec$Number == N]
          k <- mean(bgsub_hyperSpec_ijN$spc)/mean(bg_hyperSpec_ij$spc)
          if (K == FALSE) { k <- 1  }
          bgsub_hyperSpec[bgsub_hyperSpec$A == i&bgsub_hyperSpec$C == j&bgsub_hyperSpec$Number == N]  <- bgsub_hyperSpec_ijN - k*bg_hyperSpec_ij
        }
        plot(bgsub_hyperSpec[bgsub_hyperSpec$A == i&bgsub_hyperSpec$C == j][,,400~3050] )
        plot(bg_hyperSpec_ij )
      }
    }
  }

  data_hyperSpec_good.data.frame <-  data.frame(select(as.data.frame(bgsub_hyperSpec),-spc,-.row),bgsub_hyperSpec$spc)
  data_hyperSpec_good.data.frame <- na.omit(data_hyperSpec_good.data.frame)
  colnames(data_hyperSpec_good.data.frame) <- gsub("X","",colnames(data_hyperSpec_good.data.frame))
  return(data_hyperSpec_good.data.frame)
}
