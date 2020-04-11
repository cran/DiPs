\name{check}
\alias{check}
\title{Check standardized mean differences (SMDs) of the matched data set.}
\usage{
check(fdata, mdata, fz, mz)
}
\arguments{
\item{fdata}{ A full data frame with length(fz) rows and columns being variables that need to check SMDs. fdata and mdata must have the same variables with same column names in the same order.}
\item{mdata}{ A matched data frame with length(mz) rows and columns being variables that need to check SMDs. fdata and mdata must have the same variables with same column names in the same order.}
\item{fz}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control for subjects in the full data set.}
\item{mz}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control for subjects in the matched data set.}
}
\value{ A matrix with one row for each variable and five columns being the mean of treated group, mean of matched control group, mean of full control group, SMD of matched control group and SMD of full control group.
}
\description{
The function is used to create a table of mean and SMDs to check the balance before and after matching.
}
\examples{
# To run this example, you must load the optmatch package.
\donttest{
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
X<-cbind(female,age,black,education,povertyr,bmi)
p<-glm(z~female+age+black+education+povertyr+bmi,family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist<-maha_dense(d$z,X)
o<-match(d$z, dist, d)
matcheddata<-o$data
Xm<-subset(matcheddata, select=c('female','age','black','education','povertyr','bmi','p'))
check(cbind(X,p),Xm,d$z,matcheddata$z)
}
}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
