#'  A microbiome-based direction-assisted test (MiDAT) for microbiome differential abundance analyses.
#'
#'  'MiDAT_stat()' returns the test statistic values of MiDAT for a set of cadidates lambda.
#'
#'  This function is used to calculate the microbiome-based direction-assisted test (MiDAT) statistics
#'  for microbiome differential abundance analyses. MiDAT is developed to detect the overall
#'  difference of micorbiome relative abundances between two or more biological/clinical conditions.
#'
#' @param X1 relative abundance observations of OTUs for the first sample; the rows represent individuals and
#'  the columns represent OTUs.
#' @param X2 relative abundance observations of OTUs for the second sample; the rows represent individuals and
#'  the columns represent OTUs.
#' @param lambdaSet a vector containing candidate values for lambda.

#' @return A vector containing the MiDAT statistics for all candidate values of lambda.
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
#' MiDAT_stat(X1=X1, X2=X2, lambdaSet=c(1,2,3,4,Inf))
#'
#' @export
#'
MiDAT_stat <- function(X1, X2, lambdaSet){
  n1 = nrow(X1)
  n2 = nrow(X2)
  ## check the dimension of two samples
  if(ncol(X1) != ncol(X2)){
    stop("The two samples have different numbers of OTUs!")
  }
  p = ncol(X1)
  X1.bar = colMeans(X1)
  X2.bar = colMeans(X2)
  var1.jj = apply(X1, MARGIN=2, var)
  var2.jj = apply(X2, MARGIN=2, var)
  ## caclulate the sample variance for each individual OTU
  var.bar = var1.jj/n1+var2.jj/n2
  ## check whether all the observations for a single OTU are zeros
  if(min(var.bar)==0){
    stop("The individuals of two samples have equal observation for at least one OTU!")
  }
  U = (X1.bar-X2.bar)/sqrt(var.bar)
  ## divide the differences across OTUs by their directions
  Uneg = U[U<0]
  Upos = U[U>0]
  Tneg = rep(NA, length(lambdaSet))
  Tpos = rep(NA, length(lambdaSet))
  for(k in 1:length(lambdaSet)){
    if(lambdaSet[k]==Inf){
      Tneg[k] = max(abs(Uneg))
      Tpos[k] = max(abs(Upos))
    }else{
      Tneg[k] = sum(Uneg^lambdaSet[k])
      Tpos[k] = sum(Upos^lambdaSet[k])
    }
  }
  statVal = c(Tneg, Tpos)
  ## return the value of MiDAT for the candidate value of lambda
  return(statVal)
}

