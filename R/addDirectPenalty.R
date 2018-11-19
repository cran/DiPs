addDirectPenalty<-function(dist,z,dx,positive = TRUE,penalty=1){

  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))

  m <- sum(z)
  dx0 <- dx[z==0]
  dx1 <- dx[z==1]

  dif<-dx1[dist$start]-dx0[dist$end-m]
  # if (positive) d<-dist$d+dif*as.numeric(dif>0)*penalty
  # else d<-dist$d+dif*as.numeric(dif<0)*penalty
  if (positive) d<-dist$d+as.numeric(dif>0)*penalty
  else d<-dist$d+as.numeric(dif<0)*penalty

  list(d=d,start=dist$start,end=dist$end)
}
