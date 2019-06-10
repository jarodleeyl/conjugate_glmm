################################################################################
## Name:         fn_CGLMM.R
## Article:      Analysis of grouped data using conjugate generalized linear
##               mixed models
## Author:       Jarod Y. L. Lee 
## Last updated: 1 MARCH 2018
## Purpose:      Functions for fitting Poisson CGLMM, as described in Section
##               4.2 of the article.
## Input data:   NA                                                      
## Output data:  NA
## Packages:     broom - extract parameter estimates and their standard errors 
##                       from the fitted model object
################################################################################

# Calculate the log-likelihood of Poisson CGLMM.
loglik.PoisCGLMM <- function(theta)
{
  beta = theta[1:ncol(X)];  # fixed effects
  A = theta[ncol(X)+1];     # shape par
  B = 1/A;                  # scale par 
  
  y.iPLUS = tapply(y,grp,sum) 
  inter = tapply(exp(X%*%beta),grp,sum)
  ngrp = length(unique(grp))
  
  comp1 = sum((X%*%beta)*y)
  comp2 = -ngrp*A*log(B)
  comp3 = -ngrp*lgamma(A)
  comp4 = sum(lgamma(y.iPLUS+A))
  comp5 = -sum((y.iPLUS+A)*log(inter+1/B))
  
  lik = comp1+comp2+comp3+comp4+comp5
  return(lik)
}

# Calculate the gradient of the log-likelihood function of Poisson CGLMM.
grad.PoisCGLMM <- function(theta)
{
  beta = theta[1:ncol(X)];  # fixed effects
  A = theta[ncol(X)+1];     # shape par
  
  ngrp = length(unique(grp))
  y.iPLUS = tapply(y,grp,sum)
  inter = as.vector(exp(X%*%beta))
  # inter.i is sum_j e(x_{ij}\beta) for each i
  inter.i = tapply(inter,grp,sum)  
  lambda.tilde = as.vector((A+y.iPLUS)/(inter.i+A))
  # inter2 is sum_j x_{ij} e(x_{ij}\beta) for each i
  inter2 = cbind(tapply((X*inter)[,1],grp,sum), tapply((X*inter)[,2],grp,sum))
  
  comp.b1 = -apply(inter2*lambda.tilde,2,sum)
  comp.b2 = apply(X*y,2,sum)
  comp.b = comp.b1 + comp.b2
  
  comp.a = ngrp*(1+log(A)-digamma(A)) + sum(digamma(y.iPLUS+A)) - 
           sum((y.iPLUS+A)/(inter.i+A)) - sum(log(inter.i+A))
  
  grad = c(comp.b,comp.a)
  return(grad)  
}

# Fit Poisson CGLMM via maximum likelihood.
fit.PoisCGLMM <- function()
{
  theta.start = c(tidy(glm(y~x1, family=poisson, data=data))$estimate, 1)
  
  result = optim(par=theta.start, 
                 fn=loglik.PoisCGLMM, gr=grad.PoisCGLMM,
                 method="L-BFGS-B",
                 lower=c(rep(-Inf,ncol(X)),0.001),
                 upper=rep(Inf,(ncol(X)+1)),
                 control = list(fnscale = -1,maxit = 1000,trace = 0),
                 hessian=T);
  
  vc = -solve(result$hessian);
  
  return(list(
    converge = result$convergence,
    est = result$par,
    se = sqrt(diag(vc))
  ))
}

# Best predictor for estimating random effects - minimizes mean squared
# error of prediction.
bestpred.PoisCGLMM <- function(y,X,beta,A,grp)
{
  B = 1/A;
  
  num = tapply(y,grp,sum) + A 
  denom = tapply(exp(X%*%beta),grp,sum) + 1/B
    
  out = num/denom
  return(out)
}

# Calculate the deviance of Poisson regression models.
deviance.Pois <- function(mu)
{
  temp1 = y*log(y/mu)
  temp1[y==0] = 0  # to avoid 0*-Inf=NaN
  temp2 = (y-mu)
  
  temp = temp1-temp2
  deviance = 2*sum(temp)
  return(deviance)
}



