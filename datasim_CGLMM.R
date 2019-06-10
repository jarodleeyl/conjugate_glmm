################################################################################
## Name:         datasim_CGLMM.R
## Article:      Analysis of grouped data using conjugate generalized linear
##               mixed models
## Author:       Jarod Y. L. Lee 
## Last updated: 1 MARCH 2018
## Purpose:      Function for generating data for simulation study, as described 
##               in Section 4.2 of the article.
## Input data:   NA                                                      
## Output data:  NA
## Packages:     NA
################################################################################

datasim <- function(ngrp = 50000,
                    nwithin = 2,
                    b = c(0.5,1),
                    sd.re = 1,
                    distre = c("gaussian","gamma"))
{
  nwithin = rep(nwithin,ngrp)
  N = sum(nwithin)

  x1 = rep(0:1, ngrp)
  X = cbind(1,x1)
  
  if(distre == "gaussian")
  {
    re = rnorm(ngrp,0,sd.re)
    re.expand = rep(re,nwithin)
    rate = exp(X%*%b+re.expand)
  }  
  
  if(distre == "gamma")
  {
    A = 1/sd.re^2
    B = 1/A
    re = rgamma(ngrp,shape=A,scale=B)
    re.expand = rep(re,nwithin)
    rate = re.expand*exp(X%*%b)
  }
  
  y = rpois(N,rate)
  
  grp = rep(1:ngrp,nwithin)
  obs = seq(1:N)
  
  data = as.data.frame(cbind(grp,obs,X[,2],y))
  colnames(data)[3] = "x1"
  head(data)
  
  return(list(data,re))
}

