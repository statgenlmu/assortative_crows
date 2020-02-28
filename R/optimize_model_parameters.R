
funk <- function(par) {
    print(par)
    system(paste("../c++/hybridzone_simulator ",par[1]," ",par[2]," ",par[3]," ",par[1],par[2],par[3],"temp_4_ > /dev/null 2> /dev/null",sep=""))
    sout <- system(paste("../c++/genotype_llh ",par[1],par[2],par[3],"temp_4_finfreqs.csv ",par[1]," ",par[2]," ../Data/data_central_red.txt",sep=""),intern=TRUE)
    system(paste("rm ",par[1],par[2],par[3],"temp_4_finfreqs.csv ",par[1],par[2],par[3],"temp_4_simres.csv",sep=""))
    sp <- sapply(strsplit(sout,"\t"),"[",2)
    sp1 <- sapply(strsplit(sout,"\t"),"[",1)
    x <- as.numeric(sp)
    ma <- max(x,na.rm=TRUE)
    w <- which(x==ma)
    cat("Maximum",ma,"at km", sp1[w],"\n")
    write(c(par,ma),file="opp4.txt",append=TRUE)  
    -ma
}


## op <- optim(c(300,0.05,0.6),funk)
## opb <- optim(c(300,0.05,0.6),funk,method="L-BFGS-B",lower=c(0,0,0),upper=c(1e9,1,1))
## write(opb$par,file="op.txt",append=TRUE)

## op <- optim(c(6.263104,0.03678094,0.618768),funk)
## write(op$par,file="op3.txt",append=TRUE)

library(optimParallel)
cl <- makeCluster(12)
setDefaultCluster(cl=cl)
opp <- optimParallel(c(5.527374,0.02368335,0.6188977),funk,lower=c(0.01,0.01,0.01),upper=c(1000,0.99,0.99),hessian=TRUE,
                     control=list(parscale=c(100,1,1)),parallel=list(loginfo=TRUE))
write(opp$par,file="opp4.txt",append=TRUE)

funktr <- function(par) {
    print(c(par[1],exp(par[1]),par[2],1/(1+exp(-par[2])),par[3],1/(1+exp(-par[3]))))
    system(paste("../c++/hybridzone_simulator",exp(par[1]),1/(1+exp(-par[2])),1/(1+exp(-par[3])),
                 "temp_2_ > /dev/null 2> /dev/null"))
    sout <- system(paste("../c++/genotype_llh temp_2_finfreqs.csv",exp(par[1]),1/(1+exp(-par[2])),
                         "../Data/data_central_red.txt"),intern=TRUE)
    sp <- sapply(strsplit(sout,"\t"),"[",2)
    sp1 <- sapply(strsplit(sout,"\t"),"[",1)
    x <- as.numeric(sp)
    ma <- max(x,na.rm=TRUE)
    w <- which(x==ma)
    cat("Maximum",ma,"at km", sp1[w],"\n")
    write(c(exp(par[1]),1/(1+exp(-par[2])),1/(1+exp(-par[3])),ma),file="op2.txt",append=TRUE)  
    -ma
}
## op2 <- optim(c(log(300),log(0.05/(1-0.05)),log(0.6/(1-0.6))),funktr)
## write(c(op2$par,exp(par[1]),1/(1+exp(-par[2])),1/(1+exp(-par[3]))),file="op2.txt",append=TRUE)
