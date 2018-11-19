\name{maha_dense}
\alias{maha_dense}
\title{Creates a robust Mahalanobis distance for matching based on a dense network.}
\usage{
maha_dense(z, X,  exact=NULL,
nearexact=NULL, penalty=100)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{penalty}{ The penalty for a mismatch on nearexact.}
}

\value{
\item{d}{ A distance object for each pair of treated and control.}
\item{start}{ The treated subject for each distance.}
\item{end}{ The control subject for each distance.}
}
\description{
Computes a robust Mahalanobis distance list for use in dense matching. In this case, we compute the distance for all possible pairs of treated and control.

This function and its use are discussed in Rosenbaum (2010).  The robust
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
dist<-maha_dense(z,X)
head(dist$d)
detach(nh0506Homocysteine)

}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
