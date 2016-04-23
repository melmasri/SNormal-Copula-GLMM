# A class of copula GLMM to structure the 
cglmm	<-function(X,...)
	UseMethod("cglmm")
	require('sn')
	require('MASS')
		
cglmm.default <- function(X, Y, D,F,T, response.family = 'exp',tol.err=1e-5,...){
	call<-match.call()
## Checking appropriate formating
	if(!all(is.matrix(X),is.matrix(D))) 	     stop("Variable X and D must be of matrix object")
	if(!all(is.vector(Y),is.vector(T)))		     stop("T and Y must be of vector type")
	if(!is.factor(F)) 				     stop("F must be a list of factors")
	if(length(F)!=nrow(X)) 				     stop("Factors and data do not match in size") 

    xpnd <-function(x){
        l.x = length(x)
        aux = matrix(0, nrow = l.x, ncol=l.x)
        for(t in 1:l.x)  	aux[t,t:l.x]  =	x[1:(l.x-t+1)]
        t(aux)+ aux
    }
    t.normalize <-function(i)  T[F==i]
    structure( list(tol.err = 1e-5,
                    X	=	X,
                    levels 	= 	F,
                    Y	=	Y,
                    response.family = response.family,
                    D	=	D,
                    T 	=	T,
                    p  	=	ncol(X), 	# dimension of fixed effect
                    q 	= 	ncol(D),		# dimension of the random effect
                    n  	= 	length(Y),
                    units 	=	nlevels(F),
                    lambda =rep(0, NROW(X[F==1,])),
                    obs	=	 NROW(X[F==1,]),
                    xi	=	0.5,		# autocorrelation function coefficient
                    G	=	diag(ncol(D)),
                    B	=	rep(0,ncol(X)),
                    delta	=rep(0, NROW(X[F==1,])),
                    xpnd	= 	xpnd,
                    auto 	= 	function(i)  exp(i),
                    glm.link =	function(x) exp(x),
                    glm.pfun = if(response.family=='exp')
                                   function(x,y){ pexp(y ,rate = 1/(x +tol.err))} else
                                                                                    function(x,y) pgamma(y, shape =3, scale= x/3),
                    glm.fun =  if(response.family=='exp')
                                   function(x,y) sum(log(dexp(y, rate = 1/(x+tol.err), log = FALSE)+tol.err)) else
                    function(x,y) sum(log(dgamma(y, shape =3, scale= x/3, log = FALSE)+tol.err))
                    ),
              class = "cglmm"
              )
}

summary.cglmm<-function(obj){
    print(sprintf('%d observation for each of %d units for a total of %d values.',obj$obs,obj$units, obj$n))

    print(sprintf('%d fixed effect and %d random effect parameters', obj$p, obj$q))
    print(sprintf('Fixed effect design matrix of dimension %dx%d', NROW(obj$X), ncol(obj$X)))
    print(sprintf('Random effect design matrix of dimension %dx%d', NROW(obj$D), ncol(obj$D)))
    print(sprintf('Response family is %s', obj$response.family))
    print('Summary of response')
    summary(obj$Y)
    print('see print for more details.')
}

print.cglmm<-function(obj){
    obj[c('units', 'obs', 'n', 'p', 'q', 'response.family','auto', 'glm.link', 'glm.pfun', 'glm.fun', 'tol.err')]
    
}

plot.cglmm<-function(obj, mcmc, replications=50){
    org.fixed = obj$X%*%obj$B
    org.line = org.fixed + obj$D*rep(unlist(obj$b_o), each= obj$obs)
    glty=c(1)
    gcol=c('black')
    par(mfrow=c(1,2))
    hist(log(obj$Y),freq=FALSE, col = 'grey', xlab = expression(log(E(y))), main='', breaks=20)
    lines(density(org.line), main='',lty=glty, lwd = 1, xlab='')
    
    if(!missing(mcmc)){
        b.in = rep(sapply(mcmc$b, mean),each = obj$obs)
        est.line = obj$X%*%mcmc$B + obj$D*b.in
        glty=c(glty, 3)
        gcol = c(gcol, 'red')
        lines(density(est.line), main='',lty=glty[2], lwd = 1, xlab='', col=gcol[2])
    }
    legend(x='topleft', legend=c('Original', 'fit'), lty=glty, col=gcol)

    s= replications
    glty=c(1,3)
    gcol=c('black', 'red')
    plot(density(org.line), main='',lty=glty[1], lwd = 2, xlab=paste(s, 'replications'))
    for(i in 1:s){
        lines(density(org.fixed + obj$D*(rep(rnorm(obj$units, 0, sqrt(obj$G)),each=obj$obs))), main='',lty=glty[2], lwd = 1, xlab='')
        if(!missing(mcmc))
            lines(density(obj$X%*%mcmc$B+obj$D*(rep(rnorm(obj$units, 0, sqrt(mcmc$G)),
                                                each=obj$obs))),
                  main='',lty=glty[2], lwd = 1, xlab='', col=gcol[2])
    }
    legend(x='topleft', legend=c('Original', 'fit'), lty=glty, col=gcol)
}
