if(FALSE) {
    ## allele 2 can be polymorphic in the west, which means that not all crows are black there
    funk <- function(par) {
        print(par)
        if(par[2]<0 | par[2]>1 |par[3]<0 | par[3]>1) {
            return(1e10)
        }
        system(paste("../c++/hybridzone_simulator ",par[1]," ",par[2]," ",par[3]," ",par[1],par[2],par[3],"temp_sel_add_ selection additive > /dev/null 2> /dev/null",sep=""))
        sout <- system(paste("../c++/genotype_llh ",par[1],par[2],par[3],"temp_sel_add_finfreqs.csv ",par[1]," ",par[2]," ../Data/data_central_red.txt selection",sep=""),intern=TRUE)
        system(paste("rm ",par[1],par[2],par[3],"temp_sel_add_finfreqs.csv ",par[1],par[2],par[3],"temp_sel_add_simres.csv",sep=""))
        sp <- sapply(strsplit(sout,"\t"),"[",2)
        sp1 <- sapply(strsplit(sout,"\t"),"[",1)
        x <- as.numeric(sp)
        ma <- max(x,na.rm=TRUE)
        w <- which(x==ma)
        cat("Maximum",ma,"at km", sp1[w],"\n")
        write(c(par,ma),file="op_sel_add.txt",append=TRUE)  
        min(-ma,1e10)
    }
    
    if (FALSE) { 
        library(optimParallel)
        cl <- makeCluster(8)
        setDefaultCluster(cl=cl)
        opp <- optimParallel(c(2.7,0.25,0.5),funk,lower=c(0.0001,0.0000001,0.0001),upper=c(1e6,0.9999999,0.9999999),hessian=TRUE,
                             control=list(parscale=c(100,1,1)),parallel=list(loginfo=TRUE))
        
        write(opp$par,file="opp_sel_add.txt",append=TRUE)
        
    } else {
        opp <- optim(c(2.7,0.25,0.5),funk)
        write(opp$par,file="op_sel_add.txt",append=TRUE)
    }
    ## best parameters found:
} else {
    funk <- function(par) {
        ## restriction that initially all crows in the West must be black
        print(par)
        if(par[2]<0 | par[2]>1) {
            return(1e10)
        }
        system(paste("../c++/hybridzone_simulator ",par[1]," ",par[2]," 1.0 ",par[1],par[2],"temp_sel_add_bl_ selection additive > /dev/null 2> /dev/null",sep=""))
        sout <- system(paste("../c++/genotype_llh ",par[1],par[2],"temp_sel_add_bl_finfreqs.csv ",par[1]," ",par[2]," ../Data/data_central_red.txt selection",sep=""),intern=TRUE)
        system(paste("rm ",par[1],par[2],"temp_sel_add_bl_finfreqs.csv ",par[1],par[2],"temp_sel_add_bl_simres.csv",sep=""))
        sp <- sapply(strsplit(sout,"\t"),"[",2)
        sp1 <- sapply(strsplit(sout,"\t"),"[",1)
        x <- as.numeric(sp)
        ma <- max(x,na.rm=TRUE)
        w <- which(x==ma)
        cat("Maximum",ma,"at km", sp1[w],"\n")
        write(c(par,ma),file="op_sel_add_bl.txt",append=TRUE)  
        min(-ma,1e10)
    }
    
    opp <- optim(c(2.7,0.25),funk)
    write(opp$par,file="op_sel_add_bl.txt",append=TRUE)
}
