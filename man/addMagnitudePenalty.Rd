\name{addMagnitudePenalty}
\alias{addMagnitudePenalty}
\title{Add a directional magnitude penalty to a distance matrix}
\usage{
addMagnitudePenalty(dist, z, dx, positive = TRUE, hstick = 0, multiplier = 2)
}
\arguments{
\item{dist}{  A distance object with three components: d, start, end, typically created by maha_dense or maha_sparse. d[i] gives the distance between the (start[i])th treated and the (end[i]-sum(z))th control.}
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.  Must have
treated subjects (z=1) before controls (z=0).}
\item{dx}{ A vector of with length(z)=length(dx) giving the variable used to define the caliper.  For instance,
dx might be the propensity score.}
\item{positive}{ If positive = TRUE, a treated-minus-control difference di in dx that is positive is
increased by a multiple of |di|, but if positive = FALSE a a treated-minus-control difference in dx that is negative is increased by a multiple of |di|. }
\item{hstick}{ Hockey-stick value, which is a nonnegative number.  See the description.}
\item{multiplier}{ The magnitide added is multiplier*|di|/s where s is an equally weighted, pooled
within group standard deviation of dx.}
}
\value{ Returns a new distance object whose distance component d is updated by the sum of dmat and the penalties.}
\description{
Adds a penalty to the distance component d of the distance object dist depending upon value of dx.  The distance object dist has three components: d, start, end. d[i] gives the distance between the tth treated and the cth control, with t=start[i] and c=end[i]-sum(z). The value of dx for treated
unit t, say dxt, is dx[z==1][t] and the value of dx for control c, say dxc, is dx[z==0][c].  Then,
d[i] is adjusted using dxt-dxc.  If positive=TRUE, the adjustment changes d[i] to
d[i] + multiplier*(max(0, (dxt-dxc)-hstick)).  That is, a penalty is imposed if dxt exceeds dxc
by more than hstick.  If positive=FALSE, the adjustment changes d[i] to
d[i] + multiplier*(max(0, (dxc-dxt)-hstick)).
}
\examples{
\dontrun{
library(MASS)
#
data("nh0506Homocysteine")
attach(nh0506Homocysteine)
# Select covariates
X<-cbind(female,age,black,education,povertyr,bmi)
#Propensity score
p<-glm(z~female+age+black+education+povertyr+bmi,
    family=binomial)$fitted.values
d<-cbind(nh0506Homocysteine,p)
detach(nh0506Homocysteine)
dist<-maha_dense(d$z,X)
head(dist$d)
# Impose a penalty when a treated individual has a higher propensity
# score than a control
dist<-addMagnitudePenalty(dist, d$z, d$p, positive=TRUE, multiplier = 20)
head(dist$d)
}
}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
