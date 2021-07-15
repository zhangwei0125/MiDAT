#' An adaptive microbiome-based direction-assisted test (aMiDAT) for microbiome differential
#'  abundance analyses??
#'
#' 'aMiDAT()' returns the test statistic value of aMiDAT and its correpsonding p-value.
#'
#'  This function is used to calculate the adaptive microbiome-based direction-assisted test (aMiDAT)
#'  statistic and its p-value for microbiome differential abundance analyses. aMiDAT combines multiple
#'  MiDAT statistics by using their minimal p-value as the test statistic. The p-value of aMiDAT is
#'  calcuated based on an one-layer permutation procedure.
#'
#' @param X1 relative abundance observations of OTUs for the first sample; the rows represent individuals and
#'  the columns represent OTUs.
#' @param X2 relative abundance observations of OTUs for the second sample; the rows represent individuals and
#'  the columns represent OTUs.
#' @param lambdaSet a vector containing candidate values for lambda.
#' @param B total number of permutations
#'
#' @return A list with the first argument being the aMiDAT statistic value and the second argument being
#' the p-value of aMiDAT based on permutations.
#'
#' @examples
#' p = 20
#' n1 = 100
#' n2 = 100
#' mu1 = runif(p, min=0, max=10)
#' mu2 = mu1 + c(rep(0.5,5), rep(0,p-5))
#' sigMat1 = diag(p)
#' rho = 0.5
#' sigMat2 = rho^(abs(matrix(1:p, nrow=p, ncol=p, byrow=TRUE)-matrix(1:p, nrow=p, ncol=p, byrow=FALSE)))
#' Z1 = MASS::mvrnorm(n1, mu=mu1, Sigma=sigMat1); W1 = exp(Z1)
#' Z2 = MASS::mvrnorm(n2, mu=mu2, Sigma=sigMat2); W2 = exp(Z2)
#' X1 = W1/matrix(rowSums(W1), nrow=n1, ncol=p, byrow=FALSE)
#' X2 = W2/matrix(rowSums(W2), nrow=n2, ncol=p, byrow=FALSE)
#' aMiDAT(X1=X1, X2=X2, lambdaSet=c(1,2,3,4,Inf), B=1000)
#'
#' @export
#'
aMiDAT <- function(X1, X2, lambdaSet, B){
  n1 = nrow(X1)
  n2 = nrow(X2)
  ## calcuate the value of MiDAT
  Tstat =  MiDAT_stat(X1=X1, X2=X2, lambdaSet=lambdaSet)
  ### calculate the p-value using one-layer permutation procedure
  lamNum = length(lambdaSet)
  Tstat.perm = matrix(NA, nrow=B, ncol=2*lamNum)
  X = rbind(X1, X2)
  for(j in 1:B){
    ## permute the data
    samLoc = sample(1:(n1+n2), replace=F)
    X1.perm = X[samLoc[1:n1], ]
    X2.perm = X[samLoc[(n1+1):(n1+n2)], ]
    Tstat.perm[j,] = MiDAT_stat(X1=X1.perm, X2=X2.perm, lambdaSet=lambdaSet)
  }
  ### calculate the p-values for the proposed adaptive tests based on permutation
  InnerP.right = colSums(matrix(abs(Tstat), nrow=B, ncol=2*lamNum, byrow=T) <= Tstat.perm)
  InnerP.left = colSums(matrix(-abs(Tstat), nrow=B, ncol=2*lamNum, byrow=T) >= Tstat.perm)
  aTstat = min((InnerP.right+InnerP.left+1)/(B+1))  ## adaptive MiDAT (minimal p-value)
  ### calculate the p-value of the minP test statistic for one replication
  aTstat.perm = rep(NA, B)
  for(j in 1:B){
    aTstat.InnerP.right = colSums(matrix(abs(Tstat.perm[j,]), nrow=B, ncol=2*lamNum, byrow=T) <= Tstat.perm)
    aTstat.InnerP.left = colSums(matrix(-abs(Tstat.perm[j,]), nrow=B, ncol=2*lamNum, byrow=T) >= Tstat.perm)
    aTstat.perm[j] = min((aTstat.InnerP.right+aTstat.InnerP.left)/B)
  }
  aTstat.pval = mean(aTstat > aTstat.perm)

  ### return the results (adaptive MiDAT, and the aMiDAT's p-value)
  return(list(stat=aTstat, pval=aTstat.pval))
}

