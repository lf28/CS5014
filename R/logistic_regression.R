sigmoid <- function(x){
  1/(1+exp(-x))
}


logisGrad <- function(y, X, theta){
  score <- X%*%theta
  mu<-sigmoid(score)
  return(t(X)%*%(y-mu))
}


logisHessian <- function(y, X, theta){
  score <- X%*%theta
  mu <- sigmoid(score)
  dd <- mu*(1-mu)
  return(-1*(t(X)%*%diag(as.vector(dd))%*%X))
}


logisNll <- function(y, X, theta){
  score <- X%*%theta
  p <- sigmoid(score)
  -1*sum(y*log(p) + (1-y)*(log(1-p)), na.rm=T)
}

logisticDescent <- function(alpha, iter, X, y, theta0){
  thetas <- matrix(0,nrow = ncol(X), ncol= iter+1)
  thetas[,1] <- theta <- theta0
  nll <- vector(mode = "numeric", length= iter+1)
  nll[1] <- logisNll(y,X, theta)
  for(i in 1:iter){
    g <- logisGrad(y, X, theta)
    thetas[,i+1] <- theta <- theta + alpha*g
    nll[i+1] <- logisNll(y,X, theta)
  }
  return(list(thetas, nll))
}


logisNewton <- function(iter, X, y, theta0){
  thetas <- matrix(0,nrow = ncol(X), ncol= iter+1)
  thetas[,1] <- theta <- theta0
  nll <- vector(mode = "numeric", length= iter+1)
  nll[1] <- logisNll(y,X, theta)
  for(i in 1:iter){
    g <- logisGrad(y, X, theta)
    H <- logisHessian(y, X, theta)
    d <- solve(H, g)
    thetas[,i+1] <- theta <- (theta - d)
    nll[i+1] <- logisNll(y,X, theta)
  }
  return(list(thetas, nll))
}

# x_1 <- mvrnorm(100, c(-2,-2),  matrix(c(10,3,3,5), ncol = 2, byrow=T))
# x_2 <- mvrnorm(100, c(3, 3),  matrix(c(10,2,2,5), ncol = 2, byrow=T))
# X <-  cbind(1, rbind(x_1, x_2))
# y<- c(rep(0, nrow(x_1)), rep(1, nrow(x_2)))
# alpha <- 0.01
# iter <- 1000
# 
# 
# 
# # par(mfrow=c(1,3)) 
# par(mar=c(5,5,5,0.1)+.1)
# theta0 <- c(0,1,1)
# rstGD2 <- logisticDescent(alpha, iter, X, y, theta0)
# rstND2 <- logisNewton(iter, X, y, theta0)
# plot(rstGD2[[2]][1:100], type="b", xlab="Steepest GD", ylab= "NLL", col="red" ,cex.lab= 2.5)
# plot(rstND2[[2]][1:100], type="b", xlab="Newton's method", ylab= "NLL", col="red" ,cex.lab= 2.5)
# 
# library(MASS)
# x_1 <- mvrnorm(100, c(-2,-2),  matrix(c(10,3,3,5), ncol = 2, byrow=T))
# x_2 <- mvrnorm(100, c(3, 3),  matrix(c(10,4,4,5), ncol = 2, byrow=T))
# X <-  cbind(1, rbind(x_1, x_2))
# y<- c(rep(0, nrow(x_1)), rep(1, nrow(x_2)))
# alpha <- 0.01
# iter <- 1000
# 
# 
# 
# # par(mfrow=c(1,3))
# par(mar=c(5,5,5,0.1)+.1)
# theta0 <- c(0,0.1,0.1)
# rstGD2 <- logisticDescent(alpha, iter, X, y, theta0)
# plot(rstGD[[2]][1:20], type="b", xlab="Steepest GD", ylab= "NLL", col="red" ,cex.lab= 2.5)