net<-function(z,dist,ncontrol=1,fine=rep(1,length(z)),penalty=round(max(dist$d)*10000),s.cost=100){

  #check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=(ncontrol*ntreat))
  stopifnot(nobs==length(fine))

  #create basic treated-vs-control bipartite graph
  fine1=fine[z==1]
  fine0=fine[z==0]
  startn<-dist$start
  endn<-dist$end
  cost<-dist$d
  tcarcs<-length(startn) # number of treatment-control arcs
  ucap<-rep(1,tcarcs)
  b<-rep(ncontrol,ntreat) #supply for treated nodes
  b<-c(b,rep(0,ncontr)) #flow conservation at control nodes
  #Make costs integer
  if (any(cost<0)) cost<-cost-min(cost)
  #cost<-round(100*cost)
  cost<-round(cost*s.cost)

  #create a duplicate for each control to make sure each control is only used once
  startn=c(startn,(ntreat+1):nobs)
  endn=c(endn,(nobs+1):(nobs+ncontr))
  cost=c(cost,rep(0,ncontr))
  ucap=c(ucap,rep(1,ncontr))
  b<-c(b,rep(0,ncontr))

  #Add structure to the bipartite graph for near fine balance

  tb<-table(z,fine)
  nc<-as.vector(tb[1,])
  nt<-as.vector(tb[2,])
  nwant<-nt*ncontrol #desired number
  nlow<-pmin(nc,nwant) #available number
  nextra<-sum(nwant-nlow) #gap between desired and available
  finelevels<-as.vector(as.numeric(colnames(tb)))

  #Add a node for fine balance category k
  sinks<-NULL
  for (k in 1:length(nlow)){
    if (nlow[k]>0){
      sinkk<-length(b)+1
      sinks<-c(sinks,sinkk)
      who0<-fine0==finelevels[k]
      b<-c(b,0)
      if (sum(who0)>0){
        startn<-c(startn,rep(nobs,sum(who0))+which(who0))
        endn<-c(endn,rep(sinkk,sum(who0)))
        ucap<-c(ucap,rep(1,sum(who0)))
        cost<-c(cost,rep(0,sum(who0)))
      }
    }
  }

  #Add a node to take the extras
  sinkex<-length(b)+1
  b<-c(b,0)
  startn<-c(startn,(nobs+1):(nobs+ncontr))
  endn<-c(endn,rep(sinkex,ncontr))
  ucap<-c(ucap,rep(1,ncontr))
  cost<-c(cost,rep(0,ncontr))

  #Add a sink
  finalsink<-length(b)+1
  b<-c(b,-ntreat*ncontrol) #finalsink absorbs all flow
  #Connect balance nodes to finalsink
  startn<-c(startn,sinks)
  endn<-c(endn,rep(finalsink,length(sinks)))
  ucap<-c(ucap,nlow[nlow>0])
  cost<-c(cost,rep(0,length(sinks)))

  #Connect sinkex to finalsink
  startn<-c(startn,sinkex)
  endn<-c(endn,finalsink)
  #ucap<-c(ucap,ntreat*ncontrol)
  ucap<-c(ucap,nextra)
  cost<-c(cost,penalty)

  #Make costs integer
  cost<-round(cost)
  net<-list(startn=startn,endn=endn,ucap=ucap,b=b,cost=cost,tcarcs=tcarcs)
  net
}
