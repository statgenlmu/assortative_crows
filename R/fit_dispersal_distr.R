## fitting dispersal models to data in Tab 3 and its caption in
## Siefke (1994) Vogelwelt 115:83-89

## chisq(2) distribution corresponding to 2-dim normal dispersal model
logLh_chisq <- function(sig) {
    ## sig is the standard deviation of location change in each dimension
    25 * log(pchisq((5/sig)^2,2)) +
        17 * log(pchisq((10/sig)^2,2)-pchisq((5/sig)^2,2)) +
        8 * log(pchisq((15/sig)^2,2)-pchisq((10/sig)^2,2)) +
        10 * log(pchisq((20/sig)^2,2)-pchisq((15/sig)^2,2)) +
        1 * log(pchisq((30/sig)^2,2)-pchisq((20/sig)^2,2)) +
        4 * log(pchisq((40/sig)^2,2)-pchisq((30/sig)^2,2)) +
        log(dchisq((44/sig)^2,2)*2*44/sig/2) + log(dchisq((46/sig)^2,2)*2*46/sig/2) +
        log(dchisq((50/sig)^2,2)*2*50/sig/2) + log(dchisq((51/sig)^2,2)*2*51/sig/2) +
        log(dchisq((58/sig)^2,2)*2*58/sig/2) + log(dchisq((82/sig)^2,2)*2*81/sig/2)
}

op1 <- optimize(function(x) {-logLh_chisq(x)},c(10,30))

simuldata <- sqrt(apply(array(rnorm(2*71,0,sd=15.24)^2,dim=c(71,2)),1,sum))

## mix of two chisq(2) distributions, corresponding to mix of two 2-dim normal dispersal distris
logLh_chisq_mix <- function(sig1,sig2,x) {
    ## sig is the standard deviation of location change in each dimension
    25 * log(pchisq((5/sig1)^2,2)*x+pchisq((5/sig2)^2,2)*(1-x)) +
        17 * log(pchisq((10/sig1)^2,2)*x+pchisq((10/sig2)^2,2)*(1-x)-(pchisq((5/sig1)^2,2)*x+pchisq((5/sig2)^2,2)*(1-x))) +
        8 * log(pchisq((15/sig1)^2,2)*x+pchisq((15/sig2)^2,2)*(1-x)-(pchisq((10/sig1)^2,2)*x+pchisq((10/sig2)^2,2)*(1-x))) +
        10 * log(pchisq((20/sig1)^2,2)*x+pchisq((20/sig2)^2,2)*(1-x)-(pchisq((15/sig1)^2,2)*x+pchisq((15/sig2)^2,2)*(1-x))) +
        1 * log(pchisq((30/sig1)^2,2)*x+pchisq((30/sig2)^2,2)*(1-x)-(pchisq((20/sig1)^2,2)*x+pchisq((20/sig2)^2,2)*(1-x))) +
        4 * log(pchisq((40/sig1)^2,2)*x+pchisq((40/sig2)^2,2)*(1-x)-(pchisq((30/sig1)^2,2)*x+pchisq((30/sig2)^2,2)*(1-x))) +
        log(dchisq((44/sig1)^2,2)*2*44/sig1^2*x+dchisq((44/sig2)^2,2)*2*44/sig2^2*(1-x)) + log(dchisq((46/sig1)^2,2)*2*46/sig1^2*x+dchisq((46/sig2)^2,2)*2*46/sig2^2*(1-x)) +
        log(dchisq((50/sig1)^2,2)*2*50/sig1^2*x+dchisq((50/sig2)^2,2)*2*50/sig2^2*(1-x)) + log(dchisq((51/sig1)^2,2)*2*51/sig1^2*x+dchisq((51/sig2)^2,2)*2*51/sig2^2*(1-x)) +
        log(dchisq((58/sig1)^2,2)*2*58/sig1^2*x+dchisq((58/sig2)^2,2)*2*58/sig2^2*(1-x)) + log(dchisq((82/sig1)^2,2)*2*81/sig1^2*x+dchisq((82/sig2)^2,2)*2*82/sig2^2*(1-x))
}

op2 <- optim(c(10,20,0.5), function(par) -logLh_chisq_mix(par[1],par[2],par[3]),lower=c(2,15.24445,0),upper=c(15.24445,60,1),method="L-BFGS-B")
op2$par
## 6.0668316 28.1198216  0.7674638
(loglikelihoodratio <- - (op2$value - op1$objective))

par(mfcol=c(2,1))

## dispersal distances according to Siefke '94:

hist(c(rep(2.5,25),rep(7.5,17),rep(12.5,8),rep(17.5,10),25,rep(35,4),44,46,50,51,58,82),
     breaks=c(0,5,10,15,20,30,40:90),main="Fitting dispersal models to data of Siefke '94",xlab="Dispersal distance [km]",
     col="gray",ylim=c(0,0.08))

x <- 0:90
## distribution of dispersal distance in mixture model:
lines(x,dchisq((x/op2$par[1])^2,2)*2*x/op2$par[1]^2*op2$par[3]+dchisq((x/op2$par[2])^2,2)*2*x/op2$par[2]^2*(1-op2$par[3]),
      lwd=2)
abline(h=0)

lines(x,dchisq((x/op1$minimum)^2,2)*2*x/op1$minimum^2,lty=2,lwd=2)
abline(h=0)

legend("topright",lty=1:2,legend=c("resulting from mixture of two normal distr. dispersal",
                          "resulting from normal distr. dispersal"))

## distribution of east-west direction dispersal in mixture model:
x <- -60:60
plot(x,dnorm(x,0,op2$par[1])*op2$par[3]+dnorm(x,0,op2$par[2])*(1-op2$par[3]),t="l",main="Dispersal model: mixture of two normal distribution",ylab="Density",xlab="west <--- [km] ---> east")
abline(h=0)
abline(v=0)
x <- -10:10*5
b <- c(-Inf,(-10:-1)*5+2.5,(1:10)*5-2.5,Inf)
points(x,(pnorm(b[2:22],0,op2$par[1])*op2$par[3]+pnorm(b[2:22],0,op2$par[2])*(1-op2$par[3])-(pnorm(b[1:21],0,op2$par[1])*op2$par[3]+pnorm(b[1:21],0,op2$par[2])*(1-op2$par[3])))/5)
legend("topright",pch=c("-","o"),legend=c("density","bin probability"))
## dev.copy2pdf(file="dispersal_model_fit.pdf")

## Dispersal probabilities for bins:
(pnorm(b[2:22],0,op2$par[1])*op2$par[3]+pnorm(b[2:22],0,op2$par[2])*(1-op2$par[3])-(pnorm(b[1:21],0,op2$par[1])*op2$par[3]+pnorm(b[1:21],0,op2$par[2])*(1-op2$par[3])))

##  [1] 0.010601454 0.004593536 0.006005553 0.007607918 0.009340766 0.011184832
##  [7] 0.014224777 0.027894418 0.083391491 0.194232956 0.261844598 0.194232956
## [13] 0.083391491 0.027894418 0.014224777 0.011184832 0.009340766 0.007607918
## [19] 0.006005553 0.004593536 0.010601454


simulflight <- function(par) {
    if(runif(1)<par[3]) {
        return(rnorm(2,0,par[1]))
    }
    return(rnorm(2,0,par[2]))
}

ka <- array(NA,dim=c(710,2))
for(cr in 1:710) {
    ka[cr,] <- simulflight(op2$par)
}

## hist(ka[,1],60)
hist(apply(ka,1,function(x) sqrt(sum(x^2))),breaks=c(0,5,10,15,20,30,40:150))

hist(c(rep(2.5,25),rep(7.5,17),rep(12.5,8),rep(17.5,10),25,rep(35,4),44,46,50,51,58,82),
     breaks=c(0,5,10,15,20,30,40:150),main="Fitting dispersal models to data of Siefke '94",xlab="Dispersal distance [km]",
     col="gray",ylim=c(0,0.08))
