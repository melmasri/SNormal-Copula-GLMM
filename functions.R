## All functions needed for parameter estimation procedure.
## General Log-likelihood
log.lik<-function(obj,b,Psi){
    ## obj := an object from class 'cglmm'
    ## b is an m x n matrix of random effects, for m replication for each unit n
    ## Psi is the covariance matrix
    inv.psi=  pd.solve(Psi, silent=TRUE, log.det=TRUE)
    det.psi =  attributes(inv.psi)$log.det
    G_in = obj$G
    dv = unlist(lapply(obj$v, function(r) r*obj$Sigdelta))
    z = obj$Z-  dv
    aux = apply(b,1, function(s){
        dd = obj$D*s
        s = matrix(s)
        sum(tapply(z-dd, obj$levels, function(zz){
           -.5*det.psi- .5*t(zz) %*% inv.psi %*%(zz)
        }))+
            -0.5*sum(((z -dd)^2)/(1-obj$Sigdelta^2) + log(1-obj$Sigdelta^2))+
                +obj$glm.fun(x = obj$glm.link(obj$X%*%obj$B+dd),obj$Y)+
          #          ## Random effect
                    + sum(sapply(s, function(r) -0.5*log(G_in) -.5 *t(r)%*%r/obj$G))/obj$obs +
                        -sum(obj$v^2/2)
    })
    mean(aux)
}

beta.hat.optim<-function(blist, all=FALSE)	{
	### Beta Parameter
	#####################################################
	#Preliminary functions
	# Log-likelihood for the beta parameter only. A simplified version of the general likelihood.
	LLK.beta.simple<- function(beta, i, b){
		dd =  if(obj$q>1) c(obj$D[obj$levels==i,] %*% b) else obj$D[obj$levels==i]*b
		if(obj$p==1) beta=t(beta)
        ## G_in = if(is.matrix(obj$G)) det(obj$G) else obj$G
		-(
            + obj$glm.fun(x = obj$glm.link(obj$X[obj$levels==i,]%*% beta +dd) , obj$Y[obj$levels==i])
        )
	}
	LLK.beta.simple.all <-function(beta,b){
		if(obj$q>1){
                                        #sum(sapply(1:length(b), function(i) sum(sapply(1:units, function(j) LLK.beta.simple(beta, j, b[[i]][j,])))))
			sum(sapply(1:nrow(b), function(i) LLK.beta.simple(beta, i, b[i,] )))
		}else{
			sum(sapply(1:obj$units, function(i) LLK.beta.simple(beta, i, b[i])))
		}
	}
#####################################################
	beta= obj$B
	if(all){
		b_l =  lapply(1:NROW(blist[[1]]), function(i) t(sapply(1:obj$units, function(j) blist[[j]][i,])))
		sapply(b_l,function(b)
			optim(beta, fn=LLK.beta.simple.all, b = b,
                  method = "BFGS", lower = -Inf, upper = Inf,  control = list(), hessian = FALSE))

	}else{
		optim(beta, fn=LLK.beta.simple.all, b = blist, method = "BFGS", lower = -Inf, upper = Inf,  control = list(), hessian = FALSE)
	}
}

B_analytic<-function(obj,link,b){
    ## estimating fixed effect analytically
    if (obj$q>1){
        DD = lapply(1:NROW(b[[1]]),
            function(i) unlist(lapply(1:obj$units,
                                      function(j)  obj$D[obj$levels==j,]%*%b[[j]][i,])))
    }
    if (obj$q==1)
        DD = apply(b, 1, function(r) obj$D*r)
    
    if(obj$p>1){
        if(link =="exp") {
            mu  = matrix(rowMeans(apply(DD, 2, function(a) obj$SXXlY -obj$SXX%*%a )))
            temp = c(obj$X[,-1]%*%mu[-1])
            mu[1] = mean(log(obj$Y%*%exp(-DD - temp)/obj$n))
        }
        if(link =="inv"){
            mu = rowMeans(sapply(DD, function(a)obj$SXX%*%(1/obj$Y) - obj$SXX%*%a ))
            mu = matrix(rowMeans(apply(DD, 2, function(a)obj$SXX%*%(1/obj$Y) - obj$SXX%*%a )))
        }
    }
    if(obj$p==1){
        if(link == 'exp') mu = matrix(mean(log(obj$Y%*%exp(-DD)/obj$n)))
        if(link =="inv")  mu = matrix(mean(sapply(DD, function(a)obj$SXX%*%(1/obj$Y) - obj$SXX%*%a )))
    }
    mu
}

sig.beta.optim<-function(obj,b, Psi, MLE = TRUE){
    ## THE MLE Maximum estimate of a multivariate distribution
    if(MLE){
        mean(apply(b, 1, var))
    }else{
        G.in<-function(i){
            spsi=pseudoinverse(Psi)
            x = t(obj$D[obj$levels== i,])%*%spsi%*%obj$D[obj$levels== i,]
            H = tauH[i]
            solve(solve(H) - x)
        }
        mean(sapply(1:obj$units, function(i) G.in(i)))
    }
}

sig.auto <-	function(obj,p = obj$xi)
    ## Calculates the auto-regressive co-variance matrix Sigma based on xi
    round(obj$auto(obj$xpnd(obj$T[obj$levels==1]))^(-p),3)

ZGenFromY<-function(obj,blist, mean = TRUE){
    ## Given Y, estimate Z by distribution transformation
    ## mean=TRUE does the transformation T as T(E(b)), rather than E(T(b))
    ## where b is the random effect variable. The former is faster but not exact as a result of convexity.
    auxf <- function(a,b) qsn(a, xi=0, omega=1, alpha=b,lower.tail=TRUE)
    lambda_unit = rep(obj$delta/sqrt(1-obj$delta^2), obj$units)

    if(mean){
        b = colMeans(blist)
        Db = obj$D*b
        x = obj$glm.link(obj$X %*% obj$B + Db)
        pp = obj$glm.pfun(x,  obj$Y)
        pp[which(pp>.999999)]=.9999
        pp[which(pp<.000001)]=.0001
        mapply(auxf, pp, lambda_unit) + b
    } else{
        b =  apply(blist, 1,function(r) obj$D*r)
        x  = obj$glm.link(c(obj$X %*% obj$B) + b)
        pp =  obj$glm.pfun(x,  obj$Y)
        pp[which(pp>.999999)]=.9999
        pp[which(pp<.000001)]=.0001
        ll = matrix(lambda_unit, nrow=length(lambda_unit), ncol=ncol(pp))
        aux = mapply(auxf,pp, ll)
        rowMeans(matrix(aux, dim(pp)) + b)
    }
}

bGen.single<-function(obj,n, Psi){
    ## Random effect generator for each level i. Generates n copies
    spsi=  pd.solve(Psi, silent=TRUE)
    x = t(obj$D[obj$levels== 1,])%*%spsi%*%obj$D[obj$levels== 1,]
    tau2 =solve( solve(obj$G)+ x)
    if(tau2<0)
        tau2=tol.err
    dz = obj$Z  - rep(obj$Sigdelta,obj$units) *rep(obj$v, each = obj$obs)
    dd = obj$D[1:obj$obs]
    tapply(dz, obj$levels, function(zz){
        mu = c(tau2%*%t(dd)%*%spsi%*%zz)
        mvrnorm(n=n,  mu=mu, Sigma=tau2)
    }, simplify = FALSE)
}


psi.llk<-function(obj,b.ln, sig.old){
    LLK.xi.all<-function(i,b){
        D_in = obj$D[obj$levels ==i]
        dv = obj$Sigdelta *sqrt(2/pi)
        z = obj$Z[obj$levels==i]
       aux =  matrix(rowMeans(sapply(b, function(bb){
           dd = D_in*bb
           ((z - dd - dv))%*%t((z - dd - dv))
       })), ncol = obj$obs, nrow = obj$obs)
        aux
    }
    aux = sapply(1:obj$units, function(i) LLK.xi.all(i, b.ln[[i]]))
    psi.temp =  matrix(rowMeans(aux), ncol= obj$obs)
    dv = obj$Sigdelta*sqrt(2/pi)
    sig.temp = psi.temp + dv%*%t(dv)

    sig.temp  = cov2cor(sig.temp)
    hxi =-log(sig.temp)/ obj$xpnd(obj$T[obj$levels ==1])
    hxi = mean(hxi[lower.tri(hxi)], na.rm =TRUE)
    hxi = max(0.1, hxi)
    sig.temp = sig.auto(obj,p=hxi)
    
    muZ = do.call('rbind',tapply(obj$Z, obj$levels, function(r)r)) - sapply(b.ln,mean)
    muZ = colMeans(muZ)*sqrt(2/pi)
    delta = c(solve(sqrtm(sig.old))%*%muZ)
    delta= sign(delta)*pmin(0.99, abs(delta))
    Sigdelta = c(sqrtm(sig.temp)%*%delta)
    psi  = sig.temp -Sigdelta%*%t(Sigdelta)
    list(Psi = psi , xi =hxi, delta = delta, Sigdelta=Sigdelta, sig = sig.temp)
}

best_cglmm<-function(obj, b_o, lambda){
    ##    source('functions.R', local=TRUE)
    ## calculates the log-likelihood given the original random effect variable (b_o)
    Sig <-sig.auto(obj,obj$xi)
	obj$v <- abs (rnorm(obj$units,0,1))
    delta = sqrtm(Sig)%*%(lambda/sqrt(1+t(lambda)%*%lambda))
	Psi<- Sig-delta%*%t(delta)
    
    if(is.null(obj$Z))
        stop('Cannot calculate best log-likelihood, since original SN variable Z is null')
    b.ln.m  = do.call('cbind', b_o)
    b.ln.m = t(apply(b.ln.m, 1, function(r) rep(r, each= obj$obs)))
    ## Best log-likelihood
    log.lik(obj, b.ln.m, Psi)
}
