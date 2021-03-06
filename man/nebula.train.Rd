% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nebula.R
\name{nebula.train}
\alias{nebula.train}
\title{Nonparametric empirical Bayes classifier using latent annotations: wrapper function; training}
\usage{
nebula.train(pi0, pi1, n0, n1, T = NULL, I = NULL, d = 25,
  maxit = 200, tol = 1e-04, verbose = FALSE)
}
\arguments{
\item{pi0, pi1}{p x 1 vectors of control and case minor allele frequencies, respectively; IMPORTANT: must be relative to the same allele in both cases and controls}

\item{n0, n1}{number of controls and number of cases, respectively}

\item{T}{p x 1 vector of chi-square test statistics}

\item{I}{p x 1 vector of binary indicators}

\item{d}{if a single number, G0 and G1 are estimated on d x d x d grids; if a three-component vector (d0,d1,dt), G0 and G1 are estimated on d0 x d1 x dt grids}

\item{maxit}{maximum number of EM iterations}

\item{tol}{error tolerance}

\item{verbose}{TRUE to print the error attained by each EM iteration}
}
\value{
\item{type}{1=given neither T nor I; 2=given T but not I; 3=not given T but given I; 4=given both T and I}
\item{nebula}{trained classifier}
}
\description{
Nonparametric empirical Bayes classifier using latent annotations: wrapper function; training
}
\examples{
p <- 1000; ## number of snps
I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
pi1 <- pi0;
pi1[I==1] <- expit(ors+logit(pi0[I==1]));
set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,50); ## ncps
## training data
n0 <- 100; ## number of controls
X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
n1 <- 50; ## number of cases
X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
T <- rchisq(p,1,lam); ## chi-square statistics
nebula1 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,d=c(10,15));
nebula2 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T=T,d=c(10,15,20));
nebula3 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,I=I,d=c(10,15));
nebula4 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T=T,I=I,d=c(10,15,20));

}
