#' Polygenic risk score; training
#'
#' Equivalent to maximum likelihood naive Bayes classifier.
#' 
#' @param X0,X1 n x p matrix of additively coded control and case genotypes, respectively; IMPORTANT: must be coded relative to the same allele in both cases and controls
#'
#' @return
#' \item{pi0}{minor allele frequencies in controls}
#' \item{pi1}{minor allele frequencies in cases}
#' \item{P}{proportion of cases}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' prs.train(X0,X1);
#'
#' @import stats
#' @export

prs.train <- function(X0,X1){
    n0 <- nrow(X0); n1 <- nrow(X1); p <- ncol(X0);
    pi0 <- .colMeans(X0,na.rm=TRUE,n0,p)/2;
    pi1 <- .colMeans(X1,na.rm=TRUE,n1,p)/2;
    return(list(pi0=pi0,pi1=pi1,P=n1/(n0+n1)));
}

#' Polygenic risk score; prediction
#'
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param prs output of prs.train()
#' @param lambda threshold parameter; variants with log ORs with magnitudes less than lambda are not used for classification
#' @param P prevalence of cases in the testing set; if NULL, P is taken from the train object
#'
#' @return
#' \item{score}{risk score}
#' \item{class}{predicted class, 0=control, 1=case}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' prs <- prs.train(X0,X1);
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat <- prs.predict(newX,prs);
#' mean(abs(newY-Yhat$class));
#' Yhat.thresh <- prs.predict(newX,prs,lambda=0.2);
#' mean(abs(newY-Yhat.thresh$class));
#' 
#' @import stats
#' @export

prs.predict <- function(newX,prs,lambda=0,P=NULL){
    if(min(c(prs$pi0,prs$pi1))<=0){
        stop("One or more MAFs equals zero.");
    }
    if(is.null(P)){
        P <- prs$P;
    }
    b <- logit(prs$pi1)-logit(prs$pi0);
    thresh <- which(abs(b)>=lambda);
    score <- 2*sum(log((1-prs$pi1[thresh])/(1-prs$pi0[thresh])))+
        newX[,thresh]%*%b[thresh];
    return(list(score=score,class=as.numeric(score>=-logit(P))));
}
