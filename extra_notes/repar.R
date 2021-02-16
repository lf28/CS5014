trueMu <- 0
trueSig2 <- 25
m <- 1000
y<-rnorm(m, trueMu, sd= sqrt(trueSig2))
# sum of squared errors sum_{i}(yi-mu)^2
sse <- sum((y-trueMu)^2)
# closed form solution
sigm2ML <- sse/m
  
# calculate the gradient of L (gaussian loglikelihood) w.r.t reparameterised beta, where e^beta = sigma2
# m is total observation count
# sse is sum of squared error, here we assume mu is known
# beta is reparameterised parameter, and beta \in R
sigm2grad<- function(m, sse, beta){
  g<-(-m/(2*exp(beta))+sse/(2*(exp(beta)^2)))*exp(beta)
  return(g)
}

logLik <- function(m, sigm2, sse){
  return(-m/2*log(2*pi*sigm2) - 1/(2*sigm2)*sse)
}

# alpha is learning rate
# beta0 is the starting value
# m and sse are defined as before
gradientAscentForVariance <- function(iter, alpha, beta0, m, sse){
  beta <- beta0
  beta_history <- vector(mode="numeric", length = iter+1)
  beta_history[1] <- beta
  for(i in 1:iter){
    grad <- sigm2grad(m, sse, beta)
    # gradient ascent for logLik rather than descent
    beta <- beta + alpha*grad
    beta_history[i+1] <- beta
    cat("The log likelihood at ", i, "iteration is: ", logLik(m, exp(beta), sse), "\n", sep = "")
  }
  # transform back to sigma2
  sigm2_history <- exp(beta_history)
  return(list(beta_history, sigm2_history))
}

