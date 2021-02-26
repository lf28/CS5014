
xx <- seq(-6, 4, by=0.05)
ff<- expression(0.5*x^2+2*x+exp((-16)*x^2)+2*exp((-5)*x))
fx <- function(x, deriv=0){
  if(deriv == 0){
    return(0.5*x^2+2*x+exp((-16)*x^2)+2*exp((-5)*x))
  }else if(deriv==1){
    return(x+2+exp(-16*x^2)*(-16*x^2)*(-16*2*x)+2*exp(-5*x)*(-5))
  }
}
ffx<-deriv(ff, "x", hessian = T, func=T)
qf <- function(x0, x, FF=fx, FFp=ffx){
  evals <- FFp(x0)
  return(evals[1] + as.vector(attr(evals, "gradient"))*(x-x0) +  attr(evals, "hessian")[1]/2*(x-x0)^2)
}

plot(xx, fx(xx), type = "l", ylim=c(-1, 40), xlim=c(-7,4), lwd=4, ylab="f(x)", xlab="x")
lines(xx, qf(1,xx), type="l", lty=2, lwd=3, col="red")