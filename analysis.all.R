### Script for testing only, might have many mistakes

library(expm)
library(ellipse)
##################################################
##################################################
## Functions
mse<-function(sim,obj, real=FALSE){
    ## Function that returns Mean, MC SD and MSE for some paramters
    n = length(sim)
    param.name= names(sim[[1]])
    param.name = grep('(b|LLK|v|Psi|Sigdelta)', param.name, invert=TRUE, value = TRUE)
    aux = lapply(param.name, function(param){
        aux = sapply(sim, function(r) r[[param]])
        if(param == 'Z'){
            if(!real) org = obj$Z
        }else{
            if(!real)  org = obj[[param]] 
        }
        if(!is.null(dim(aux))){
            m = rowMeans(aux)
            sd = sqrt((rowMeans(aux^2)-rowMeans(aux)^2)*n/(n-1))
            if(!real)
                mse =colMeans((org - t(aux))^2)
        }else{
            m = mean(aux)
            sd = sd(aux)
            if(!real)
                mse = mean((org - aux)^2)
        }
        cbind(m, sd, if(!real) mse else 1)
    })
    names(aux)<-param.name
    aux
}

## General coverage probability
coverage.prop<-function(x,v,fisher, a=.05){
    if(is.vector(fisher)){
        v <= x+	qnorm(1-a/2)*sqrt(1/fisher) &  v>= x +	qnorm(a/2)*sqrt(1/fisher)
    }else if(is.matrix(fisher)){
        e= ellipse(fisher, scale = c(1, 1), centre =x, level = 1-a)
        d = sqrt(sum((x-v)^2))
        dall = sqrt(rowSums((e - v)^2))
        dclose= e[which.min(dall),]
        dclose = sqrt(sum((dclose - x)^2))
        d<=dclose
    }
}


files = list.files(path = ".", pattern = 'SIM')
sink('results_table.txt')
for (f in files){
    load(f)
    param.mse = mse(mcmc.sim, lmm.in)
    ## Coverage probability for Beta
    if(length(lmm.in$B)>1){
        fisher.B = solve(t(lmm.in$X[,-1])%*%lmm.in$X[,-1])
        EC.B = mean(sapply(mcmc.sim, function(r){ coverage.prop(r[['B']][-1], lmm.in$B[-1], fisher.B )}))
        EC.B = rep(EC.B, length(lmm.in$B))
        fisher.B = c(t(lmm.in$X[,1])%*%lmm.in$X[,1])
        EC.B[1] =mean(sapply(mcmc.sim, function(r){ coverage.prop(r[['B']][1], lmm.in$B[1], fisher.B)}))
        alpha = mean(sapply(mcmc.sim, function(r){ coverage.prop(r[['B']][1] + mean(unlist(r[['b']])), lmm.in$B[1], fisher.B)}))
        EC.B = c(alpha, EC.B)
    }else {
        fisher.B = c(t(lmm.in$X)%*%lmm.in$X)
        EC.B = mean(sapply(mcmc.sim, function(r) coverage.prop(r[['B']], lmm.in$B, fisher.B)))
        EC.B = cbind(EC.B, mean(sapply(mcmc.sim, function(r) coverage.prop(r[['B']] + mean(unlist(r[['b']])), lmm.in$B, fisher.B))))
    }
    ## Coverage probability for G
    fisher.sig.b = 2.5*lmm.in$units/param.mse$G[1,'m']
    EC.G = mean(1*sapply(mcmc.sim, function(r){ coverage.prop(r[['G']], lmm.in$G, fisher.sig.b )}))
    ## Beta
    options(digits=5)
    ##B = cbind(lmm.in$B, param.mse$B + Eb['m'], EC.B)
    alpha = sapply(mcmc.sim, function(r) r[['B']][1] + mean(unlist(r[['b']])))
    alpha = c(mean(alpha), sd(alpha), mean((alpha -lmm.in$B[1])^2))
    xx = rbind(alpha, param.mse$B)
    B = cbind(c(lmm.in$B[1],lmm.in$B),xx, c(EC.B))
    rownames(B)<-c('alpha', paste('B',2:NROW(B)-1, sep=''))
    ## G
    G = cbind(lmm.in$G, param.mse$G, EC.G)
    rownames(G)<-'G'
    ## E[b]
    Eb  = sapply(mcmc.sim, function(r){ a = do.call(rbind, r[['b']]);  colMeans(a)       })
    Eb =  c(0, m = mean(Eb) ,sd = sd(Eb), NA, NA)
    ## xi
    Xi = c(xi =lmm.in$xi[1], param.mse$xi, NA)
    ## lambda
    Lam= rowMeans(sapply(mcmc.sim, function(r) r[['delta']]/(sqrt(1-r[['delta']]%*%r[['delta']]))))

    Lam = c(1, mean(Lam), sd(Lam) , mean((Lam - 1)^2), NA)
    ## All
    print(f)
    tb = rbind(B, Eb, G, Xi, Lam)
    print(tb)
    print('##################################')
    ## 1st
}
sink()
q('no')
