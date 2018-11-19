\name{maha_sparse}
\alias{maha_sparse}
\title{Creates a robust Mahalanobis distance for matching based on a sparse network.}
\usage{
maha_sparse(z, X, p=rep(1,length(z)), caliper=1, exact=NULL,
nearexact=NULL, penalty=100)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{p}{ A vector of length(z)=length(p) giving the variable used to define the caliper.  Typically, p is the propensity score or its rank.}
\item{caliper}{ If two individuals differ on p by more than caliper, we will not calculate the distance for this pair. }
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{penalty}{ The penalty for a mismatch on nearexact.}
}

\value{
\item{d}{ A distance list for each pair within the caliper distance and constant constraint.}
\item{start}{ The treated subject for each distance.}
\item{end}{ The control subject for each distance.}
}
\description{
Computes a robust Mahalanobis distance list for use in sparse matching. In this case, we will only calculate the distance for pairs within the caliper on p. If the caliper is too small, the matching may be infeasible. For the smallest caliper that keeps feasiblibity, refer to optcal() in package 'bigmatch'.

This function and its use are discussed in Rosenbaum (2010). It is preferred when the dataset is large. The robust
Mahalanobis distance in described in Chapter 8 of Rosenbaum (2010).

}
\details{
The usual Mahalanobis distance works well for multivariate Normal covariates, but can exhibit odd behavior with typical covariates. Long tails or an outlier in a covariate can yield a large estimated
variance, so the usual Mahalanobis distance pays little attention to large differences in this covariate. Rare binary covariates have a small variance, so a mismatch on a rare binary covariate is
viewed by the usual Mahalanobis distance as extremely important. If you were matching for binary covariates indicating US state of residence, the usual Mahalanobis distance would regard a
mismatch for Wyoming as much worse than a mismatch for California.

The robust Mahalanobis distance uses ranks of covariates rather than the covariates themselves, but the variances of the ranks are not adjusted for ties, so ties do not make a variable more important. Binary covariates are, of course, heavily tied.
}
\examples{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
\donttest{
#apply symmetric caliper 0.15 on propensity score
dist1<-maha_sparse(d$z,X,p,0.15)
length(dist1$d)
}
#apply symmetric caliper 0.1 on propensity score
dist2<-maha_sparse(d$z,X,p,0.1)
length(dist2$d)
}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
