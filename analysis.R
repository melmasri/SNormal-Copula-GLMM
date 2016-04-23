### Script for testing only, might have many mistakes
library(xtable)
library(expm)
library(ellipse)
##################################################
##################################################
## Functions
mse<-function(sim, real=FALSE){
    ## Function that returns Mean, MC SD and MSE for some paramters
    ## sim = mcmc.sim
    n = length(sim)
    param.name= names(sim[[1]])
    param.name = grep('(b|LLK|v|Psi)', param.name, invert=TRUE, value = TRUE)
    aux = lapply(param.name, function(param){
        print(param)
        aux = sapply(sim, function(r) r[[param]])
        if(param =='xi'){
            param='xi'
            aux = aux[1,]
            if(!real)
                org = lmm.in[[param]][1]
        }else if( param == 'Z'){
            if(!real) org = unlist(Z_o)
        }else{
            if(!real)  org = lmm.in[[param]] 
               }
        if(!is.null(dim(aux))){
            m = rowMeans(aux)
            sd = apply(aux,1,sd)
            if(!real)
                mse =rowMeans((org - aux)^2)
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

##################################################
##################################################
### Analysis
PRINT = FALSE

load(FILENAME)
param.mse = mse(mcmc.sim)
aux = sapply(mcmc.sim, function(r) r[['itra']])
hist(aux)
summary(aux)

aux1 = sapply(mcmc.sim, function(r)sum(r[['LLK']]))
hist(aux1)
LLK_o
plot(aux,aux1)

## Single replication
if(PRINT) pdf('ExpExpSingleSingle.pdf')
i = sample(length(mcmc.sim), 1)
i = 46
org.line = lmm.in$X*c(lmm.in$B) + lmm.in$D*rep(unlist(b_o), times = lmm.in$obs)
plot(density(org.line), main='',lty=1, lwd = 2, xlab='')
b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
lines(density(lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in), lty=2, lwd=1)
if(PRINT) dev.off()


## multi replication
if(PRINT)pdf('ExpExpSingle100.pdf')
s = sample(length(mcmc.sim), 50)
plot(density(org.line), main='',lty=1, lwd = 2, xlab='')
for(i in s){
    b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
    lines(density(lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in), lty=3, lwd=1)
}
if(PRINT)    dev.off()


### Plotting Z versus a sample of Z
dz = density(unlist(Z_o))
plot(dz, col='red', xlab = 'Z', main ='')
#s = sample(100, 100)
for(i in s) lines(density(mcmc.sim[[i]]$Z), col='black')
plot(density(unlist(Z_o)), col='red', xlab = 'Z', main='')
lines(dz$x, rowMeans(sapply(mcmc.sim, function(r) density(r[['Z']])$y)))

### plot Y versus a sample of Y
plot(density(log(lmm.in$Y)), lwd=2, xlab = expression(log(y)), main='')
s = sample(length(mcmc.sim), 10)
for(i in s){
    b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
    ## lines(density(lmm.in$X%*%mcmc.sim[[i]]$B))
    lines(density(lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in ), lwd=1, lty =3)
}

plot(density(log(lmm.in$Y)), col = 'red', xlab = expression(log(y)), main='')
b.in= rep(rowMeans(sapply(mcmc.sim, function(r) sapply(r[['b']], mean ))), times= lmm.in$obs)
b.in = rep(sapply(mcmc.sim[[1]]$b, mean),times = lmm.in$obs)
aux = lmm.in$X%*%param.mse$B[,'m'] + lmm.in$D*b.in
aux = lmm.in$X%*%param.mse$B[,'m']
lines(density(aux))

##################################################
##################################################
##################################################

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
    EC.B = mean(sapply(mcmc.sim, function(r) coverage.prop(r[['B']]+ mean(unlist(r[['b']])), lmm.in$B, fisher.B)))
    EC.B = cbind(EC.B, mean(sapply(mcmc.sim, function(r) coverage.prop(r[['B']] , lmm.in$B, fisher.B))))
}
## Coverage probability for G
fisher.sig.b = 2.5*lmm.in$units/lmm.in$G
EC.G = mean(1*sapply(mcmc.sim, function(r){ coverage.prop(sqrt(r[['G']]), sqrt(lmm.in$G), fisher.sig.b )}))
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
obj=lmm.in
aux = sapply(mcmc.sim, function(r){
    x = r[[xi]][1]
    sig  = sig.auto(1, x)
    d = r[['delta']][1:obj$obs[1]]
    a = solve(sqrtm(sig))%*%d
    a = a/sqrt(1-a)
})
Lam = colMeans(aux)
Lam = c(1, mean(Lam), sd(Lam) , mean((Lam - 1)^2), NA)
## All
tb = rbind(B, Eb, G, Xi, Lam)
print(tb)

xtable(tb, digits=4)
print('##################################')





##################################################
##################################################
##################################################
##################################################


## Real Data
param.mse = mse(mcmc.sim, real=TRUE)

hist(log(lmm.in$Y), main ='', freq = FALSE)
s = sample(length(mcmc.sim), 10)
for(i in s){
    b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
    lines(density(log(3) +lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in))
}

Y = lmm.in$Y
par(mfrow=c(1,2))
dz = density(unlist(mcmc.sim[[1]]$Z))
plot(dz, col='red', xlab = 'Z', main ='')
#s = sample(100, 100)
for(i in s) lines(density(mcmc.sim[[i]]$Z), col='black')
plot(density(unlist(mcmc.sim[[1]]$Z)), xlab = 'Z', main='')
lines(dz$x, rowMeans(sapply(mcmc.sim, function(r) density(r[['Z']])$y)))

### plot Y versus a sample of Y
par(mfrow=c(1,2))
hist(lmm.in$Y,freq=FALSE, col = 'red', xlab = expression(log(y)), main='', xlim = c(0,10))
#s = sample(100, 100)

for(i in s){
    b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
    lines(density(lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in- log(3)))
}

plot(density(log(lmm.in$Y)), col = 'red', xlab = expression(log(y)), main='')
b.in= rep(rowMeans(sapply(mcmc.sim, function(r) sapply(r[['b']], mean ))), times= lmm.in$obs)
b.in = rep(sapply(mcmc.sim[[1]]$b, mean),times = lmm.in$obs)
aux = lmm.in$X%*%param.mse$B[,'m'] + lmm.in$D*b.in
lines(density(aux))


   ##lines(density(lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in))



plot(density(Y*100), main ="" , ylab= "", xlab= "")
for(i in s){
    b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
    lines(density(lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in- log(3)))
}




hist(Y, xlab = "cholesterol levels", ylab ="",main ="", col='grey', labels=FALSE,border='white',freq = FALSE, ylim = c(0,.012),yaxt='n')

for(i in s){
    b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
    lines(density((100)*exp(lmm.in$X%*%mcmc.sim[[i]]$B + lmm.in$D*b.in)))
}


for(i in 1:100) lines(density( 10*(XL%*%Bhat + bhatgen(sdbhat) )+256), lty=2, lwd=1)


par(mfrow=c(1,1))
hist(Y, xlab = "cholesterol levels", ylab ="",main ="", col='grey', labels=FALSE,border='white',freq = FALSE,yaxt='n')
for(i in s[1:10]){
    ##    b.in = rep(sapply(mcmc.sim[[i]]$b, mean),times = lmm.in$obs)
    b.in = rep(rnorm(lmm.in$units, m = 0, sd = sqrt(mcmc.sim[[i]]$G)), times = lmm.in$obs)
    ##    lines(density(lmm.in$X%*%mcmc.sim[[i]]$B +b.in))
    lines(density(exp(lmm.in$X%*%mcmc.sim[[i]]$B)))
}



## Plotting Real Data
LastQ = sapply(mcmc.sim, function(r) sum(r[['LLK']]))
est = which.max(LastQ)
## Based on Average
pdf('CholHist.pdf')
hist(100*Y, xlab = "cholesterol levels", ylab ="",main ="", col='grey', labels=FALSE,border='white',freq = FALSE,yaxt='n', ylim=c(0,0.012))
lines(density(100*exp(lmm.in$X%*%param.mse$B[,'m'])),lty=1, lwd=1)
dev.off()


## pdf('CholHist100.pdf')
## hist(100*Y, xlab = "cholesterol levels", ylab ="",main ="", col='grey', labels=FALSE,border='white',freq = FALSE,yaxt='n',ylim=c(0,0.012))

## org = exp(lmm.in$X%*%param.mse$B[,'m'])
## aux = sapply(mcmc.sim, function(r) mean(abs(org - lmm.in$X%*%r[['B']])^2))
## s = order(aux, decreasing = FALSE)
## for(i in s[1:10]){
##     lines(density(100*exp(lmm.in$X%*%mcmc.sim[[i]]$B)),lty=2, lwd=1)
## }


## Beta
options(digits=3)
B = param.mse$B
rownames(B)<-paste('B', 1:length(lmm.in$B), sep='')
## G
G = param.mse$G
rownames(G)<-'G'
## xi
Xi = param.mse$xi
rownames(Xi)<-'Xi'
## All
rbind(B, G, Xi)
xtable(rbind(B, G, Xi))

BIC=-2*log(-LastQ) +(lmm.in$p+2)*log(sum(lmm.in$obs))
AIC=2*(length(lmm.in$p)+2)-2*log(-LastQ)
mean(LastQ)
mean(AIC)
mean(BIC)


## Based on Max LLK
options(digits=3)
B = cbind(mcmc.sim[[est]]$B, 1/sqrt(colSums(lmm.in$X^2)))
rownames(B)<-paste('B', 1:length(lmm.in$B), sep='')
## G
G = cbind(mcmc.sim[[est]]$G, 1/sqrt(2.5*lmm.in$units* mcmc.sim[[est]]$G))
rownames(G)<-'G'
## xi
Xi = cbind(mcmc.sim[[est]]$xi[1], NA)
rownames(Xi)<-'Xi'
## All
rbind(B, G, Xi)
xtable(rbind(B, G, Xi))

BIC=-2*log(-LastQ) +(lmm.in$p+2)*log(sum(lmm.in$obs))
AIC=2*(length(lmm.in$p)+2)-2*log(-LastQ)
LastQ[est]
AIC[est]
BIC[est]


pdf('CholHist.pdf')
hist(100*Y, xlab = "cholesterol levels", ylab ="",main ="", col='grey', labels=FALSE,border='white',freq = FALSE,yaxt='n', ylim=c(0,0.012))
lines(density(100*exp(lmm.in$X%*%mcmc.sim[[est]]$B)),lty=1, lwd=1)
dev.off()

## hist(100*Y, xlab = "cholesterol levels", ylab ="",main ="", col='grey', labels=FALSE,border='white',freq = FALSE,yaxt='n', ylim=c(0,0.012))

## lines(density(100*exp(lmm.in$X%*%mcmc.sim[[est]]$B+ rep(colMeans(b), times = lmm.in$obs)) ),lty=1, lwd=1)

## b = do.call('cbind', mcmc.sim[[est]]$b)
## s = sample(1:nrow(b), 20)
## for(i in s){
##     lines(density(100*exp(lmm.in$X%*%mcmc.sim[[est]]$B)*rep(b[i,], times = lmm.in$obs) ),lty=1, lwd=1)
## }
