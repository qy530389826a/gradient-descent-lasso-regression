rm(list=ls(all=T))
install.packages('glmnet')
library(glmnet)
#simulate a data p>n only first 5 columns have infulence to y the other can be treated as noise
n <- 100
p <- 300
sigma <- 0.5
beta <- c(rep(1, 5), rep(0, p - 5))

X <- matrix(nrow = n, ncol = p, rnorm(n * p))
e <- rnorm(n, mean = 0, sd = sigma)
y <- X %*% beta + e

# fit the model with different lambda
length=101
lambda = exp(seq(log(1), log(0.005), length = 101))
beta_estimate<- matrix(0,length,p)

# our original model
ptm_mymodel<-system.time(for(l in 1:length){
  b <- rep(0,p)
  r=y-X%*%b
  error=1
  if(error>0.0001){
    temp<-b
    for(j in 1:p){
      #update r in place
      r <- r+X[,j]*b[j]
      Xj<-(sum(X[,j]*X[,j]))^(-1)
      z<-sum(X[,j]*r)
      lam<-lambda[l]*n
      b[j] <- (abs(z)-lam)*Xj
      b[j] <- sign(z)*ifelse(b[j]>0,b[j],0)
      r <- r-X[,j]*b[j]
    }
    error<-sum(temp-b)^2
  }
  beta_estimate[l,]<-b
}
)



## lasso path generate
pct <- rowSums(abs(beta_estimate))/sum(abs(beta_estimate[2,]))
matplot(pct,beta_estimate,type="l",lty=1,
        xlab="L 1 norm",ylab="Coefficients")
text(1.02,beta_estimate[1,],1:p,cex=1,col=1:p)

## result from package
ptm_package<- system.time(fit1 <- glmnet(
  X,
  y,
  intercept = FALSE,
  standardize = FALSE,
  lambda = lambda
))
plot(fit1)


## compare computation time with package glmnet
print(cbind(ptm_package,ptm_mymodel))

