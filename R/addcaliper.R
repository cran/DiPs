addcaliper<-function(dist,z,dx,rg, stdev = FALSE, penalty = 1000){

  stopifnot(is.vector(rg)&(length(rg)==2))
  stopifnot((rg[1]<=0)&(rg[2]>=0))


  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))

  dx1<-dx[z==1]
  dx0<-dx[z==0]
  # Standardize p using an equally weighted pooled variance
  v1<-stats::var(dx1)
  v2<-stats::var(dx0)
  sp<-sqrt((v1+v2)/2)
  stopifnot(sp>0)

  if (stdev) rg <- rg *sp

  m=sum(z)
  dif<-dx1[dist$start]-dx0[dist$end-m]
  d <- dist$d + as.numeric(dif>rg[2])*penalty
  d <- d + as.numeric(dif<rg[1])*penalty

  list(d=d,start=dist$start,end=dist$end)
}
