## A script made to generate a random dataset to fit using the model.

## storing all environment object names before the script
obj.names = ls()

tol.err=1e-5      # tolerance error 
xpnd <-	function(x){
    l.x = length(x)
    aux = matrix(0, nrow = l.x, ncol=l.x)
    for(t in 1:l.x)  	aux[t,t:l.x]  =	x[1:(l.x-t+1)]
    t(aux)+ aux
}

auto 		<- 	function(i)  exp(i)
sig.auto	<-	function(p=xi) round(auto(xpnd(T[F==1]))^(-p),3)
glm.link	<- 	function(x)  exp(x)

## Time variable
T 	<- 1*unlist(lapply(1:units, function(i){ 
    x = sample(1:obs, obs, replace = FALSE, prob = NULL)
    x = x[sort.list(x)]
    (x - min(x)) }))

F 	<- factor(rep(1:units, each=obs))    # unit factors
G 	<- sig_b                             # covariance(variance) of random effect variable
Sig <- sig.auto(xi)                      # covariance of Z

k	<- 3                                # the gamma parameter

b_o = lapply(1:units,function(r) mvrnorm(1,mu = rep(0,NROW(G)),Sigma=G))# original used b
Db	= lapply(1:units, function(i) if(NROW(b)>1) c(D[F==i,]%*%b_o[[i]]) else  D[F==i]*b_o[[i]])
## @@ Inverting to get Y from Z @@Y
alpha = lambda/sqrt(1+ t(lambda)%*%lambda)
delta = alpha
lambda_unit = delta/sqrt(1-delta^2)
Z_o = rmsn(units, xi = 0, Omega = Sig , alpha = lambda)
pz =  apply(Z_o, 1, function(r) psn(r, xi =0, omega=1, alpha =lambda_unit))

if(response.family=='exp')
    Y = lapply(1:units, function(i){
        qexp(pz[,i], rate = 1/glm.link(X[F==i,]%*%B +Db[[i]]), lower.tail = TRUE)
    })

if(response.family=='gamma')
    Y = lapply(1:units, function(i){
        qgamma(pz[,i], shape =k, scale= glm.link(X[F==i,]%*%B +Db[[i]])/k)
    })

Z_o = lapply(1:units, function(r) Z_o[r,] + b_o[[r]])
Y= unlist(Y)
# Default initialization
lmm.in <- cglmm(X, Y, D,F,T,tol.err=tol.err, response.family=response.family)
# Saving original set-up initialization 
lmm.in$B = B
lmm.in$G = G
lmm.in$xi = xi
lmm.in$delta = alpha
lmm.in$Sigdelta = sqrtm(Sig)%*%alpha
lmm.in$Z <-unlist(Z_o)
lmm.in$b_o <- b_o

## removing undesired variables.
rm(list= grep('lmm.in', setdiff(ls(), obj.names), invert=TRUE, value=TRUE))
