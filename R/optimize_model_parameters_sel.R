
funk <- function(par) {
    print(par)
    system(paste("../c++/hybridzone_simulator ",par[1]," ",par[2]," ",par[3]," ",par[1],par[2],par[3],"temp_sel9_ selection > /dev/null 2> /dev/null",sep=""))
    sout <- system(paste("../c++/genotype_llh ",par[1],par[2],par[3],"temp_sel9_finfreqs.csv ",par[1]," ",par[2]," ../Data/data_central_red.txt selection",sep=""),intern=TRUE)
    system(paste("rm ",par[1],par[2],par[3],"temp_sel9_finfreqs.csv ",par[1],par[2],par[3],"temp_sel9_simres.csv",sep=""))
    sp <- sapply(strsplit(sout,"\t"),"[",2)
    sp1 <- sapply(strsplit(sout,"\t"),"[",1)
    x <- as.numeric(sp)
    ma <- max(x,na.rm=TRUE)
    w <- which(x==ma)
    cat("Maximum",ma,"at km", sp1[w],"\n")
    write(c(par,ma),file="opp_sel9.txt",append=TRUE)  
    -ma
}

if (FALSE) { ## set to TRUE to use parallel optimization package
    library(optimParallel)
    cl <- makeCluster(12)
    setDefaultCluster(cl=cl)
    opp <- optimParallel(c(2.7,0.25,0.45),funk,lower=c(0.0001,0.0000001,0.0001),upper=c(1e6,0.9999999,0.9999999),hessian=TRUE,
                         control=list(parscale=c(100,1,1)),parallel=list(loginfo=TRUE))
    
    write(opp$par,file="opp_sel8.txt",append=TRUE)

    
} else {
    opp <- optim(c(2.171128,0.6172853,0.7766008),funk) ## starting values were results from previous optimization runs
    write(opp$par,file="opp_sel9.txt",append=TRUE)
}
