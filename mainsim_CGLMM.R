################################################################################
## Name:         mainsim_CGLMM.R
## Article:      Analysis of grouped data using conjugate generalized linear
##               mixed models
## Author:       Jarod Y. L. Lee 
## Last updated: 1 MARCH 2018
## Purpose:      Script for conducting the simulation study, as described in
##               Section 4.2 of the article.
## Input data:   NA                                                      
## Output data:  RData file consisting of the following objects:
##               result.GLM   - simulation ...
##               result.LAP   - results ...   
##               result.AGQ2  - of ...
##               result.AGQ5  - each ...
##               result.AGQ10 - iteration ...
##               result.CGLMM - ...
##               result.summary.mean - averaged results across simulations
## Packages:     lme4 - fit GLMM (Laplace, Quadrature) 
##               rbenchmark - compare the running time of various methods
##               broom - extract parameter estimates and their standard errors 
##                       from the fitted model object
################################################################################

# == Clear R memory ============================================================
rm(list=ls()) 

# == Load the required packages ================================================
libname <- c("lme4","rbenchmark","broom")
lapply(libname, require, character.only=T)

# == Source the required functions =============================================
source("datasim_CGLMM.R")
source("fn_CGLMM.R")

# == Set the simulation parameters =============================================
seed = 1234
nsim = 1000
ngrp = 50000
nwithin = 2
b = c(0.5,1)
sd.re = 1
distre = "gamma"
#distre = "gaussian"

if(distre == "gamma")     filename = "resultGamma_1000sim.RData"
if(distre == "gaussian")  filename = "resultGaussian_1000sim.RData"

# == Initialize the results matrices ===========================================
result.GLM = matrix(rep(0,5*nsim), nrow=5, ncol=nsim)
result.LAP = result.GLM 
result.AGQ2 = result.GLM 
result.AGQ5 = result.GLM 
result.AGQ10 = result.GLM 
result.CGLMM = result.GLM

# == Simulation ================================================================
set.seed(seed)
for(i in 1:nsim)
{
  cat("Simulation "); cat(i); cat("; "); cat("\n");
  
  # Simulate data
  datsim = datasim(ngrp, nwithin, b, sd.re, distre)
  
  data = datsim[[1]]
  re = datsim[[2]]
  true.int = b[1] + re
  
  grp = data$grp
  X = cbind(1,data$x1)
  y = data$y
  N = nrow(data)
  
  # Fit the various Poisson models
  timebenchmark = benchmark(
    # GLM
    "GLM" = {
      cat("Fitting GLM. \n")
      GLM = glm(y~x1, family=poisson, data=data)
    },
    
    # GLMM (Laplace Approximation)
    "LAP" = {
      cat("Fitting GLMM (Laplace Approximation). \n")
      LAP = glmer(y~x1+(1|grp), family=poisson, data=data,
                  control = glmerControl(check.conv.grad = .makeCC("warning", tol = 0.1, relTol = NULL)))
    },
    
    # GLMM (2 adaptive Gauss-Hermite quadrature points)
    "AGQ2" = {
      cat("Fitting GLMM (Quadrature - 2 points). \n")
      AGQ2 = glmer(y~x1+(1|grp), family=poisson, nAGQ = 2, data=data,
                   control = glmerControl(check.conv.grad = .makeCC("warning", tol = 0.1, relTol = NULL)))
    },
    
    # GLMM (5 adaptive Gauss-Hermite quadrature points)
    "AGQ5" = {
      cat("Fitting GLMM (Quadrature - 5 points). \n")
      AGQ5 = glmer(y~x1+(1|grp), family=poisson, nAGQ = 5, data=data,
                   control = glmerControl(check.conv.grad = .makeCC("warning", tol = 0.1, relTol = NULL)))
    },
    
    # GLMM (10 adaptive Gauss-Hermite quadrature points)
    "AGQ10" = {
      cat("Fitting GLMM (Quadrature - 10 points). \n")
      AGQ10 = glmer(y~x1+(1|grp), family=poisson, nAGQ = 10, data=data,
                    control = glmerControl(check.conv.grad = .makeCC("warning", tol = 0.1, relTol = NULL)))
    },
    
    # CGLMM
    "CGLMM" = {
      cat("Fitting CGLMM. \n")
      CGLMM = fit.PoisCGLMM();
    },
    
    columns = c("test", "replications", "elapsed",  "user.self", "sys.self",
                "relative"),
    order = NULL,
    replications = 1
  )
  
  # Extract summaries from the fitted models
  est.GLM = c(tidy(GLM)$estimate, NA)
  deviance.GLM = deviance(GLM)
  result.GLM[,i] = c(est.GLM, deviance.GLM, timebenchmark$elapsed[1])
  
  est.LAP = tidy(LAP)$estimate
  deviance.LAP = deviance(LAP)
  result.LAP[,i] = c(est.LAP, deviance.LAP, timebenchmark$elapsed[2])
  
  est.AGQ2 = tidy(AGQ2)$estimate
  deviance.AGQ2 = deviance(AGQ2)
  result.AGQ2[,i] = c(est.AGQ2, deviance.AGQ2, timebenchmark$elapsed[3])
  
  est.AGQ5 = tidy(AGQ5)$estimate
  deviance.AGQ5 = deviance(AGQ5)
  result.AGQ5[,i] = c(est.AGQ5, deviance.AGQ5, timebenchmark$elapsed[4])
  
  est.AGQ10 = tidy(AGQ10)$estimate
  deviance.AGQ10 = deviance(AGQ10)
  result.AGQ10[,i] = c(est.AGQ10, deviance.AGQ10, timebenchmark$elapsed[5])
  
  est.CGLMM = c(CGLMM$est[1:2],sqrt(1/CGLMM$est[3]))
  #se.CGLMM = CGLMM$se      # se = diag(sqrt(-solve(hess(CGLMM$est))))
  repred.CGLMM = bestpred.PoisCGLMM(y,X,
                                    beta=CGLMM$est[1:ncol(X)],
                                    A=CGLMM$est[ncol(X)+1],
                                    grp)
  repred.CGLMM.expand = rep(repred.CGLMM, rep(nwithin,ngrp))
  fitted.CGLMM = repred.CGLMM.expand*exp(X%*%CGLMM$est[1:ncol(X)])
  deviance.CGLMM = deviance.Pois(fitted.CGLMM)
  result.CGLMM[,i] = c(est.CGLMM, deviance.CGLMM, timebenchmark$elapsed[6])
}

rowname = c("b0","b1","sigma","deviance","time")
colname = c("GLM","LAP","AGQ2","AGQ5","AGQ10","CGLMM")

rownames(result.GLM) = rowname
rownames(result.LAP) = rowname
rownames(result.AGQ2) = rowname
rownames(result.AGQ5) = rowname
rownames(result.AGQ10) = rowname
rownames(result.CGLMM) = rowname

# == Summarizes the results ====================================================
summary.mean.GLM = apply(result.GLM,1,mean)
summary.mean.LAP = apply(result.LAP,1,mean)
summary.mean.AGQ2 = apply(result.AGQ2,1,mean)
summary.mean.AGQ5 = apply(result.AGQ5,1,mean)
summary.mean.AGQ10 = apply(result.AGQ10,1,mean)
summary.mean.CGLMM = apply(result.CGLMM,1,mean)

result.summary.mean = cbind(summary.mean.GLM, summary.mean.LAP, 
                            summary.mean.AGQ2, summary.mean.AGQ5, 
                            summary.mean.AGQ10, summary.mean.CGLMM)
result.summary.mean[4,] = result.summary.mean[4,]/(ngrp*nwithin)
rownames(result.summary.mean) = rowname
colnames(result.summary.mean) = colname


# == Print the summarized results ==============================================
cat("Estimates: "); cat("\n");
print(round(result.summary.mean,2)); cat("\n");


# == Save the results ==========================================================
save(result.GLM, result.LAP, result.AGQ2, result.AGQ5, result.AGQ10, 
     result.CGLMM, result.summary.mean,
     file = filename)


