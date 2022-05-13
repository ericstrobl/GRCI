earth_wrap <- function(x,y,z,suffStat){
  
  x1=suffStat$data[,x];
  y1=suffStat$data[,y];
  z1=suffStat$data[,z,drop=FALSE];
  
  out = earth_test(x1,y1,z1)
  
  return(out$p)
  
}