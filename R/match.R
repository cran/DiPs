match<-function(z,fine=rep(1,length(z)),dist,dat,ncontrol=1,penalty=round(max(dist$d)*10000),s.cost=100){
  #Check input
  stopifnot(is.data.frame(dat))
  stopifnot(is.vector(z))
  if (is.factor(fine)){
    levels(fine)<-1:nlevels(fine)
    fine<-as.integer(fine)
  }
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  n<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=(ncontrol*ntreat))
  stopifnot(length(z)==length(fine))

  stopifnot(length(z)==(dim(dat)[1]))


  #Must have treated first
  n<-length(z)
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    warning("Your input data has been sorted to place treated (z=1) subjects first.  This was done in an effort
            to maintain compatibility between dist and rownames of dat,
            assuming the distance dist was built using maha_dense or maha_sparse.")
    o<-order(1-z)
    z<-z[o]
    dat<-dat[o,]
    fine<-fine[o]
  }

  #do match
  if (!requireNamespace("optmatch", quietly=TRUE)) {
    stop("Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.")
  }

  net<-net(z,dist,ncontrol,fine,penalty,s.cost)
  if (any(net$cost==Inf)) net$cost[net$cost==Inf]<-2*max(net$cost[net$cost!=Inf])
  output<-rcbalance::callrelax(net)
  if (output$feasible!=1){
    warning("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
    m<-list(feasible=output$feasible,d=NULL)
  }else{
    x<-output$x[1:net$tcarcs]
    treated<-net$startn[1:net$tcarcs]
    control<-net$endn[1:net$tcarcs]
    match.df<-data.frame('treat' = treated, 'x' = x, 'control' = control)
    matched.or.not<-plyr::daply(match.df, plyr::.(match.df$treat),
                                  function(treat.edges) c(as.numeric(as.character(treat.edges$treat[1])), sum(treat.edges$x)), .drop_o = FALSE)
    if(any(matched.or.not[,2] == 0)){
      match.df<-match.df[-which(match.df$treat %in% matched.or.not[which(matched.or.not[,2] == 0),1]),]
    }
    match.df$treat<-as.factor(as.character(match.df$treat))
    matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
                                     function(treat.edges) treat.edges$control[treat.edges$x == 1], .drop_o = FALSE))
    id1<-(1:n)[z==1]
    id0<-(1:n)[z==0]
    matchid<-matrix(c(id1[as.numeric(row.names(matches))], id0[as.vector((matches-sum(z)))]),ncol=ncontrol+1)
    matchid<-as.vector(t(matchid))
    dat1<-dat[matchid,]
    zm<-z[matchid]
    mset<-rep(1:ntreat,each=ncontrol+1)
    dat1<-cbind(mset,dat1)
    m<-list(feasible=output$feasible,data=dat1,x=x)
  }
  #m<-nearfine(z,dist,dat,X,ncontrol,fine,penalty,max.cost,nearexPenalty,dx,disto)
  if(m[[1]]==0) {
    warning("The match you requested is infeasible, reconsider caliper or ncontrol or exact for distance.")
  }
  m
}
