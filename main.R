## Main function running the estimation procedure
rm( list = ls())
options(digits=3)
## libraries
library(sn)
library(MASS)
library(expm)
library(mnormt)
###########################
### Initialization ########
###########################

##################################################
## Parameter set-up
## Uni variate fixed effect
B_o = B= matrix(3, nrow=1,ncol=1)       # Fixed effect 
units	= 200                           # no. of units
obs 	=5                              # observation per unit
sig_b  = 4                              # variance of random effect
xi =0.2                                 # retrogressive coef.
b = c(1)                                # a vector indicating the size of the random effect
X 	= if (NROW(B)>1) cbind(1,matrix(rnorm((NROW(B)-1)*obs*units,0,0.2),ncol=NROW(B)-1  )) else matrix(1,nrow=obs*units ) # Design matrix X
D	= matrix(X[,1:length(sig_b)],nrow=nrow(X),ncol=length(sig_b) ) # random effect design matrix
lambda = rep(1, obs)                                               # skewness parameter
response.family = 'exp'                 # possible types 'exp' or 'gamma'
source('cglmm.R')                       # cglmm object
source('random.data.R')                 # generating a random response Y and the Skew vector Z in a variable called lmm.in
source('functions.R')       # supporting function
LLK_o = best_cglmm(lmm.in, lmm.in$b_o, lambda)  # calculating the approximate best LLK, 

##################################################
## Multivariate fixed effect
## Arelano-Valle (2005) Set-up, see paper
B_o= B= c(1, 2,1)                       # Fixed effect
units	= 200                           # no. of units
obs 	=5                              # observations per unit
X 		= cbind(alpha = 1, t= rep(obs:1 -3, units), w =rep(c(1,0),each =obs*units/2))
b= c(1)
sig_b  = if(NROW(b)>1) diag(c(1,1)) else c(4) # variance of random effect
xi =0.2                                 # For some values the algorithm is not stable.
D	= matrix(1, ncol=1, nrow= nrow(X))  # random effect design matrix
lambda = rep(1, obs)                    # universal skewness vector
response.family = 'exp'         # 'exp' or 'gamma'
source('cglmm.R')
source('random.data.R')                 # creates a random dataset called lmm.in
source('functions.R')
LLK_o = best_cglmm(lmm.in, lmm.in$b_o, lambda)

##################################################
##################################################
## MC-EM estimation for a single chain
source('mcem.cglmm.R')                   # MC-EM loop
print(sprintf('Best LLK %0.3f', LLK_o)) # best LLK
mcmc.sim = mcem.cglmm(obj=lmm.in,itra = 40,verbose=TRUE, sink.to.file=FALSE, update.single.param = TRUE)

mcmc.sim$B
mcmc.sim$B[1] + mean(unlist(mcmc.sim$b))
mcmc.sim$G
mcmc.sim$xi

plot(lmm.in, mcmc.sim) 

## To run for multiple chains
## Settings for server
## options(echo=TRUE) # if you want see commands in output file
## args <- commandArgs(trailingOnly = TRUE)
## DATAFILE = args[1]
## load(DATAFILE)

## To run for multiple chains
no.chains = 50
no.cores = 13
library(parallel)
mcmc.sim = mclapply(1:no.chains,function(x,ob){
    source('mcem.cglmm.R', local=TRUE)
    source('functions.R', local=TRUE)
    tol.err <<- 1e-5
    LLK_o = best_cglmm(lmm.in, lmm.in$b_o, lambda)  # calculating the approximate best LLK, 
    mcem.cglmm(obj= ob,itra=74,verbose=FALSE, sink.to.file=FALSE,
               update.single.param = TRUE)}
   ,ob=lmm.in,mc.preschedule = TRUE, mc.cores =no.cores) 

## Saving for analyzing
fname   = paste0('SIM_',format(Sys.time(), "-%d-%m-%Hh%M_"), DATAFILE)
save.image(file =fname)


## to receive and email notification
subject = '"End of processing"'
#######################################################
body = paste("Your processing is done for" , fname)
## Sending notification
email = paste("echo",body," . | mailx -s " , subject ," myname@mymail.com")
system(email)
q('no')

