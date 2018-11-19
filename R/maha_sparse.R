maha_sparse<-function(z,X,p=rep(1,length(z)),caliper=1,exact=NULL,nearexact=NULL,penalty=100){
  Xmatrix<-function(x){
    if (is.vector(x) || is.factor(x)) x<-matrix(x,nrow=length(z))

    if(is.data.frame(x) || is.character(x)){
      if(!is.data.frame(x)) x <- as.data.frame(x)
      X.chars <- which(plyr::laply(x, function(y) 'character' %in% class(y)))
      if(length(X.chars) > 0){
        for(i in X.chars){
          x[,i] <- factor(x[,i])

        }
      }
      #if some variables are factors convert to dummies
      X.factors <-  which(plyr::laply(x, function(y) 'factor' %in% class(y)))

      #handle missing data
      for(i in which(plyr::laply(x, function(y) any(is.na(y))))){
        if(i %in% X.factors){
          #for factors, make NA a new factor level
          x[,i] <- addNA(x[,i])
        }else{
          #for numeric/logical, impute means and add a new indicator for missingness
          x[[paste(colnames(x)[i],'NA', sep = '')]] <- is.na(x[,i])
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }
      for(i in rev(X.factors)){
        dummyXi <- model.matrix(as.formula(
          paste('~',colnames(x)[i], '-1')),data=x)
        x <- cbind(x[,-i], dummyXi)
      }

    }else{
      #handle missing data
      for(i in c(1:ncol(x))){
        if(any(is.na(x[,i]))){
          x <- cbind(x,is.na(X[,i]))
          colnames(x)[ncol(x)] <- paste(colnames(X)[i],'NA', sep = '')
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }

    }

    #get rid of columns that do not vary
    varying <- apply(x,2, function(y) length(unique(y)) > 1)
    x <- x[,which(varying),drop = FALSE]

    as.matrix(x)
  }

  #Check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(all((z==1)|(z==0)))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(length(z)==length(p))
  X<-Xmatrix(X)

  if (is.factor(nearexact)){
    levels(nearexact)=1:nlevels(nearexact)
    nearexact<-as.integer(nearexact)
  }

  stopifnot(length(z)==(dim(X)[1]))

  if (!is.null(exact)){
    stopifnot(length(exact)==length(z))
    tb<-table(z,exact)
    if (!all(tb[2,]<=tb[1,])){
      warning("An exact match for exact is infeasible for every caliper.")
      return(NULL)
    }
  }

  stopifnot(caliper>=0)

  # Standardize p using an equally weighted pooled variance
  v1<-stats::var(p[z==1])
  v2<-stats::var(p[z==0])
  sp<-sqrt((v1+v2)/2)

  n<-dim(X)[1]

  # #Must have treated first
  # if(!(min(z[1:(n-1)]-z[2:n])>=0)){
  #   o<-order(1-z)
  #   z<-z[o]
  #   p<-p[o]
  #   X<-X[o,]
  #   if (!is.null(exact)) exact<-exact[o]
  #   if (!is.null(nearexact)) nearexact<-nearexact[o]
  # }

  # if (is.vector(X)) X<-matrix(X,ncol=1)

  ids<-1:n
  k<-dim(X)[2]
  m<-sum(z)

  for (j in 1:k) X[, j]<-rank(X[, j])
  cv<-stats::cov(X)
  vuntied<-stats::var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  LL<-chol(cv)

  if (is.vector(X)) X<-matrix(X,ncol=1)
  X0<-X[z==0,]
  X1<-X[z==1,]
  if (is.vector(X0)) X0<-matrix(X0,ncol=1)
  if (is.vector(X1)) X1<-matrix(X1,ncol=1)
  cids <-ids[z==0]
  p0<-p[z==0]
  p1<-p[z==1]
  if (!is.null(exact)){
    exact0<-exact[z==0]
    exact1<-exact[z==1]
  }
  if (!is.null(nearexact)){
    nearex0<-nearexact[z==0]
    nearex1<-nearexact[z==1]
  }

  distance<-c()
  start_node<-c()
  end_node<-c()
  nearex<-c()

  for (i in 1:m){
    #use caliper
    d<-abs(p1[i]-p0)
    who<-(d<=caliper)
    #use exact
    if (!is.null(exact)) who<-who&(exact1[i]==exact0)
    if (!any(who)) return(NULL)
    num<-sum(who)
    cc<-X0[who,]
    if (num==1) cc<-matrix(cc,nrow=1)
    else if (is.vector(cc)) cc<-matrix(cc,nrow=num)
    tt<-t(as.matrix(X1[i,]))
    distancei<-mvnfast::maha(cc,tt,LL,isChol=TRUE)

    distance<-c(distance, distancei)
    start_node<-c(start_node, rep(i, num))
    end_node<-c(end_node,cids[who])

    #use nearexact
    if (!is.null(nearexact)){
      nearex<-c(nearex, nearex1[i]!=nearex0[who])
      out<-list(d=distance+nearex*penalty,start=start_node,end=end_node)
    }else{
      out<-list(d=distance,start=start_node,end=end_node)
    }
  }
  out
}
