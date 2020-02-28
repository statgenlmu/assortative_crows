
## here, the programs hybridzone_simulator and genotype_llh must be compiled with
## initialize_w(argv[1],argv[2]) replaced by initialize_w_pure(atof(argv[1]),atof(argv[2]))

funk <- function(par) {
    print(par)
    rrr <- sample(.Machine$integer.max,1)
    system(paste("../c++/hybridzone_simulator ",par[1]," ",par[2]," ",par[3]," ",rrr,"_temp_sel_pur_ selection > /dev/null 2> /dev/null",sep=""))
    sout <- system(paste("../c++/genotype_llh ",rrr,"_temp_sel_pur_finfreqs.csv ",par[1]," ",par[2]," ../Data/data_central_red.txt selection",sep=""),intern=TRUE)
    system(paste("rm ",rrr,"_temp_sel_pur_finfreqs.csv ",rrr,"_temp_sel_pur_simres.csv",sep=""))
    sp <- sapply(strsplit(sout,"\t"),"[",2)
    sp1 <- sapply(strsplit(sout,"\t"),"[",1)
    x <- as.numeric(sp)
    ma <- max(x,na.rm=TRUE)
    w <- which(x==ma)
    cat("Maximum",ma,"at km", sp1[w],"\n")
    write(c(par,ma),file="opp_sel_p.txt",append=TRUE)  
    -ma
}

if (FALSE) { 
    library(optimParallel)
    cl <- makeCluster(12)
    setDefaultCluster(cl=cl)
    opp <- optimParallel(c(0.5,0.5,0.9),funk,lower=c(0.001,0.001,0.001),upper=c(1.0,1.0,1.0),hessian=TRUE,
                         control=list(parallel=list(loginfo=TRUE)))
    
    write(opp$par,file="opp_sel_pp.txt",append=TRUE)


    cl <- makeCluster(12)
    setDefaultCluster(cl=cl)
    opp <- optimParallel(c(0.7,0.2,0.99),funk,lower=c(0.001,0.001,0.001),upper=c(1.0,1.0,1.0),hessian=TRUE,
                         control=list(parallel=list(loginfo=TRUE)))
    
    write(opp$par,file="opp_sel_pp.txt",append=TRUE)
    
} else {
    opp <- optim(c(0.5,0.5,0.9),funk,lower=c(0.001,0.001,0.001),upper=c(1.0,1.0,1.0),method="L-BFGS-B")
    write(opp$par,file="opp_sel_pur.txt",append=TRUE)
    
    opp <- optim(c(0.7,0.2,0.99),funk,lower=c(0.001,0.001,0.001),upper=c(1.0,1.0,1.0),method="L-BFGS-B")
    write(opp$par,file="opp_sel_pur.txt",append=TRUE)
}
