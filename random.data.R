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

## @@ MASS LIBRARY @@: mvrnorm(n = 1, mu, Sigma, tol = 1e-6, empirical = FALSE)
b_o = lapply(1:units,function(r) mvrnorm(1,mu = rep(0,NROW(G)),Sigma=G))# original used b
Db	= lapply(1:units, function(i) if(NROW(b)>1) c(D[F==i,]%*%b_o[[i]]) else  D[F==i]*b_o[[i]])

## @@ Inverting to get Y from Z @@Y
## alpha = lambda/sqrt(1+ t(lambda)%*%lambda)
## delta = sqrtm(Sig)%*%alpha
## lambda_unit = delta/sqrt(1-delta^2)
## Z_o = rmsn(units, xi = 0, Omega = Sig , alpha = lambda)
## Using a normal
delta = lambda/sqrt(1+ t(lambda)%*%lambda)
Sigdelta = c(sqrtm(Sig)%*%delta)
lambda_unit = delta/sqrt(1-delta^2)
v= abs(rnorm(units,0,1))
Psi = Sig - Sigdelta%*%t(Sigdelta)
Z_o = do.call('rbind',lapply(1:units, function(r)rmnorm(1, mean = Sigdelta*v[r], varcov = Psi) ))
##Z_o= rmnorm(units, mean = Sigdelta*v, varcov = Psi)
pz =  apply(Z_o, 1, function(r) psn(r, xi =0, omega=1, alpha =lambda_unit))

if(response.family=='exp')
    Y = lapply(1:units, function(i){
        qexp(pz[,i], rate = 1/glm.link(X[F==i,]%*%B +Db[[i]]), lower.tail = TRUE)
    })

if(response.family=='gamma')
    Y = lapply(1:units, function(i){
        qgamma(pz[,i], shape =k, scale= glm.link(X[F==i,]%*%B +Db[[i]])/k)
    })
Y= unlist(Y)

Z_o = lapply(1:units, function(r) Z_o[r,] + b_o[[r]])

# Default initialization
lmm.in <- cglmm(X, Y, D,F,T,tol.err=tol.err, response.family=response.family)
# Saving original set-up initialization 
lmm.in$B = B
lmm.in$G = G
lmm.in$xi = xi
lmm.in$delta = delta
lmm.in$Sigdelta = Sigdelta
lmm.in$Z <-unlist(Z_o)
lmm.in$b_o <- b_o

## removing undesired variables.
rm(list= grep('lmm.in', setdiff(ls(), obj.names), invert=TRUE, value=TRUE))

