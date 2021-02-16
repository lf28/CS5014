logisGrad <- function(y, X, theta){
  score <- X%*%theta
  mu<-pracma::sigmoid(score)
  return(t(X)%*%(y-mu))
}

logisticDescent <- function(alpha, iter, X, y, theta0){
  thetas <- matrix(0,nrow = ncol(X), ncol= iter+2)
  thetas[,1] <- theta <- theta0
  for(i in 1:iter){
    g <- logisGrad(y, X, theta)
    thetas[,i+1] <- theta <- theta + alpha*g
  }
  return(thetas)
}