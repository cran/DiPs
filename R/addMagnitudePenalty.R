addMagnitudePenalty<-function(dist,z,dx,positive=TRUE,hstick=0,multiplier=2){

  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))

  if (hstick<0) hstick<-(-hstick)

  m <- sum(z)
  dx0 <- dx[z==0]
  dx1 <- dx[z==1]
  # Standardize p using an equally weighted pooled variance
  v1<-stats::var(dx1)
  v2<-stats::var(dx0)
  sp<-sqrt((v1+v2)/2)
  stopifnot(sp>0)
  dx <- dx/sp
  dx0 <- dx[z==0]
  dx1 <- dx[z==1]

  dif<-dx1[dist$start]-dx0[dist$end-m]
  if (!positive) dif <- (-dif)
  if (hstick!=0) dif<-dif-hstick
  d<-dist$d+dif*as.numeric(dif>0)*multiplier

  list(d=d,start=dist$start,end=dist$end)
}
