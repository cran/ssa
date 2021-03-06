% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prs.R
\name{prs.predict.cv}
\alias{prs.predict.cv}
\title{Polygenic risk score; prediction for classifier trained with CV}
\usage{
prs.predict.cv(newX, prs, P = NULL)
}
\arguments{
\item{newX}{n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls}

\item{prs}{output of prs.train.cv()}

\item{P}{prevalence of cases in the testing set; if NULL, P is taken from the train object}
}
\value{
\item{score}{risk score}
\item{class}{predicted class, 0=control, 1=case}
}
\description{
Polygenic risk score; prediction for classifier trained with CV
}
\examples{
p <- 1000; ## number of snps
I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
pi1 <- pi0;
pi1[I==1] <- expit(ors+logit(pi0[I==1]));
## training data
n0 <- 100; ## number of controls
X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
n1 <- 50; ## number of cases
X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
prs <- prs.train.cv(X0,X1,K=3,w=1,nlambda=100,verbose=TRUE);
## testing data
newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
              t(replicate(n1,rbinom(p,2,pi1))));
newY <- c(rep(0,n0),rep(1,n1));
Yhat <- prs.predict.cv(newX,prs);
mean(abs(newY-Yhat$class));

}
