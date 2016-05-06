# A function that holds the MC-EM loop and the estimation procedure
mcem.cglmm<- function(obj, itra =20,verbose=FALSE, delta, sink.to.file=FALSE, update.single.param = FALSE){
    ## Validation conditions

    ## global setup variable
    n <- 50                             # Initial number or replications
    n.max <- 3000                       # Maximum number of replication
    n.increase.freq <- 15                 # The frequency of increasing n
    n.increase.factor <- 3                # Increase factor of n as in n = x*n, where x is the increase factor
    start.single.param.update <- 20     # when to start updating single parameters
    link <- 'exp'                                   # The only supported link function for GLM so far

    ## local needed functions
    thresh<-function(x,n,t) min(t*(-x + n)/(n-1),0) #smoothing function.
    param.update<-function(name, i)     # to update or not the parameter in 'name' in iteration 'i'
        ((!update.single.param|i<start.single.param.update) | (update.single.param & param[(i-1)%%3 +1] ==name & i >=start.single.param.update))
    
    ## writing to file
    if(sink.to.file)  sink(paste0('est ',format(Sys.time(),"%Y%m%d %H-%M-%S"),'.csv'))

    ## Setting up the starting values of parameters.
    obj$SXX =solve(t(obj$X)%*%obj$X)%*%t(obj$X)
	obj$SXXlY= -obj$SXX%*%log(1/obj$Y)
	b.ln = lapply(1:obj$units, function(r) mvrnorm(n, mu=rep(0,obj$q),Sigma=1 ))
    b.ln.m  = do.call('cbind', b.ln)
    b.ln.m = t(apply(b.ln.m, 1, function(r) rep(r, each= obj$obs)))
    update.Z = if(is.null(obj$Z)) TRUE else FALSE
    obj$Z<- if(update.Z) ZGenFromY(obj,b.ln.m, mean=TRUE) else obj$Z

    obj$B[1:NROW(obj$B)]=1;B_best<-rep(1, obj$p)
	G_best= obj$G =1*diag(obj$q)
	Sig<- round(sig.auto(obj,p=obj$xi),2)
    sig.old = Sig
    xi_best<-obj$xi
    lam= rep(1, obj$obs)
    delta_best<-obj$delta<-lam/sqrt(1+t(lam)%*%lam)
    obj$Sigdelta<-sqrtm(Sig)%*%obj$delta
    Psi<- Sig-obj$Sigdelta%*%t(obj$Sigdelta)
    obj$v <- abs (rnorm(obj$units,0,1))
	v_best <-obj$v
    LLK_best = -1e10
    BIC_best = -1e10
    
    ## ME-MC loop    
	for(i in 1:itra){
        ## Increasing n
        if(i%% n.increase.freq==0) n=n.increase.factor*n
        ## Maximum allowed increase
        if(n>n.max) n=n.max
        
        ## Selecting a parameter to update
        if((i-1)%% 3==0){
            param = sample(c('b', 'g', 'xi'), 3)
            if(i==1) param = c('b', 'g', 'xi')
            if(update.single.param & verbose) print(paste('Updating order', param, sep= ' '))
        }

        ## updating parameters
        if(param.update('b', i)){
            ## Estimating Beta analytically for exponential or gamma marginals with exponential link function
            out = 	B_analytic(obj,link,b.ln.m);
            obj$B	<- 	out
        }  
        
        if(param.update('xi', i)){
            ## Estimating xi
            Psi_new <- psi.llk(obj,b.ln, sig.old)
            Psi <- Psi_new$Psi
            sig.old = Psi_new$sig
            obj$delta <-Psi_new$delta
            obj$Sigdelta<-Psi_new$Sigdelta
            obj$xi = Psi_new$xi
        }
        
        if(param.update('g', i)){
            ## Estimating G
            out_G = sig.beta.optim(obj,b.ln.m,Psi, MLE = TRUE); obj$G<- out_G
        }    

        ## Calculating likelihood and BIC
        auxlikelihood = log.lik(obj,b.ln.m,Psi)
        BIC=-2*log(abs(sum(auxlikelihood))) +(obj$p+2)*log(sum(obj$obs))

        ## Printing results on screen
        if(verbose | sink.to.file){
            tb.head<-c('itr', 'n','LLK', 'BIC',
                       paste('B', 1:obj$p-1, sep=''),
                       paste('G', 1:obj$q, sep=''), 'xi',
                       paste('d', 1:obj$obs,sep=''))

            tb.new = c(i, n,sum(auxlikelihood), BIC,
                obj$B, if(obj$q>1) diag(obj$G) else obj$G,obj$xi,
                obj$d)

            tb.best = c(i,n,sum(LLK_best), BIC_best,
                B_best, if(obj$q>1) diag(G_best) else G_best, xi_best,
                delta_best)

            tb = rbind(tb.new,tb.best)
            colnames(tb)<-tb.head 
            rownames(tb)<-c('update', 'best')
            print(tb)
        }

        ## Storing updates
        if((i<3)|(BIC - BIC_best > thresh(i,25, -0.05))){
            ## Recording best iteration
            Z_best = obj$Z
            B_best = obj$B
            b_best = b.ln
            G_best = obj$G
            xi_best = obj$xi
            Psi_best =Psi
            delta_best = obj$delta
            Sigdelta_best = obj$Sigdelta
            v_best <-obj$v
            LLK_best = auxlikelihood
            i_best = i
            BIC_best=BIC
        }

        # reverting to best updates
        obj$G = G_best
        obj$B = B_best
        obj$delta<-delta_best
        obj$xi = xi_best
        Sigdelta_best = obj$Sigdelta
        Psi =Psi_best
        obj$v <- abs (rnorm(obj$units,0,1))
        b.ln <-bGen.single(obj,n, Psi)
        b.ln.m = t(apply(do.call('cbind', b.ln), 1, function(r) rep(r, each= obj$obs)))
        if(update.Z)  obj$Z<-ZGenFromY(obj,b.ln.m, mean=FALSE) ## no need to update Z if the original is given 

        
    }
    if(sink.to.file)   sink()           # closing sink
    ## return estimated parameters.
	list(B=B_best,G = G_best,xi=xi_best, delta= delta_best,Sigdelta = Sigdelta_best,LLK = LLK_best ,Z = Z_best, v= v_best, b=b_best, Psi = Psi_best, itra = i_best)
}



