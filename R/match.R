match<-function(z,p,dist,dat,exact=NULL,fine=rep(1,length(z)),ncontrol=1,penalty=round(max(dist$d)*1000),s.cost=100,subX=NULL){
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

  if (!is.null(subX)){
    if (is.factor(subX)){
      levels(subX)<-1:nlevels(subX)
      subX<-as.integer(subX)
    }
    stopifnot(is.vector(subX))
  }

  #sort input
  if (is.null(exact)){
    o<-order(1-p)
  }else{
    o<-order(exact,1-p)
    exact<-exact[o]
  }

  z<-z[o]
  p<-p[o]
  fine<-fine[o]
  dat<-dat[o,]

  #Must have treated first
  if(!(min(z[1:(n-1)]-z[2:n])>=0)){
    o<-order(1-z)
    z<-z[o]
    dat<-dat[o,]
    fine<-fine[o]
    if (!is.null(subX)) subX<-subX[o]
  }

  #do match
  if (!requireNamespace("optmatch", quietly=TRUE)) {
    stop("Error: package optmatch (>= 0.9-1) not loaded.  To run match command, you must install optmatch first and agree to the terms of its license.")
  }

  net<-net(z,dist,ncontrol,fine,penalty,s.cost,subX)
  if (any(net$cost==Inf)) net$cost[net$cost==Inf]<-2*max(net$cost[net$cost!=Inf])

  callrelax <- function (net) {
    if (!requireNamespace("optmatch", quietly = TRUE)) {
      stop('Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.')
    }
    startn <- net$startn
    endn <- net$endn
    ucap <- net$ucap
    b <- net$b
    cost <- net$cost
    nnodes <- length(b)
    my.expr <- parse(text = '.Fortran("relaxalg", nnodes, as.integer(length(startn)),
                     as.integer(startn), as.integer(endn), as.integer(cost),
                     as.integer(ucap), as.integer(b), x1 = integer(length(startn)),
                     crash1 = as.integer(0), large1 = as.integer(.Machine$integer.max/4),
                     feasible1 = integer(1), NAOK = FALSE, DUP = TRUE, PACKAGE = "optmatch")')
    fop <- eval(my.expr)
    x <- fop$x1
    feasible <- fop$feasible1
    crash <- fop$crash1
    list(crash = crash, feasible = feasible, x = x)
  }

  #output<-rcbalance::callrelax(net)
  output<-callrelax(net)

  if (output$feasible!=1){
    warning("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
    m<-list(feasible=output$feasible,d=NULL)
  }else{
    x<-output$x[1:net$tcarcs]
    treated<-net$startn[1:net$tcarcs]
    control<-net$endn[1:net$tcarcs]
    treated=treated[which(x==1)]
    control=control[which(x==1)]
    match.df=data.frame('treat'=treated,'control'=control)
    match.df$treat<-as.factor(as.character(match.df$treat))
    matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
                                   function(treat.edges) treat.edges$control,.drop_o=FALSE))
    id1<-(1:n)[z==1]
    id0<-(1:n)[z==0]
    matchid<-matrix(c(id1[as.numeric(row.names(matches))], id0[as.vector((matches-sum(z)))]),ncol=ncontrol+1)
    matchid<-as.vector(t(matchid))
    dat1<-dat[matchid,]
    zm<-z[matchid]
    mset<-rep(1:nrow(matches),each=ncontrol+1)
    dat1<-cbind(mset,dat1)
    m<-list(feasible=output$feasible,data=dat1,x=x)
  }
  #m<-nearfine(z,dist,dat,X,ncontrol,fine,penalty,max.cost,nearexPenalty,dx,disto)
  if(m[[1]]==0) {
    warning("The match you requested is infeasible, reconsider caliper or ncontrol or exact for distance.")
  }
  m
}
