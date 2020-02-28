par(mfcol=c(1,1))

## for 
w <- matrix(c(1,0.01,0.01,0.0229039,0.01,0.01,0.01,0.01,0.01,
0.01,1,0.01,0.609487,0.0286452,0.01,0.01,0.01,0.01,
0.01,0.01,1,0.01,0.648728,0.2143,0.0290243,0.0290243,0.0290243,
0.0229039,0.609487,0.01,1,0.01,0.01,0.01,0.01,0.01,
0.01,0.0286452,0.648728,0.01,1,0.0271598,0.01,0.01,0.01,
0.01,0.01,0.2143,0.01,0.0271598,1,0.663688,0.663688,0.663688,
0.01,0.01,0.0290243,0.01,0.01,0.663688,1,1,1,
0.01,0.01,0.0290243,0.01,0.01,0.663688,1,1,1,
0.01,0.01,0.0290243,0.01,0.01,0.663688,1,1,1),ncol=9,
            dimnames=
            list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
                 c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))

## for simulresults_sigma10_mmp10perc.csv:
w <- matrix(c(1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
0.1,1,0.1,0.609487,0.1,0.1,0.1,0.1,0.1,
0.1,0.1,1,0.1,0.648728,0.2143,0.1,0.1,0.1,
0.1,0.609487,0.1,1,0.1,0.1,0.1,0.1,0.1,
0.1,0.1,0.648728,0.1,1,0.1,0.1,0.1,0.1,
0.1,0.1,0.2143,0.1,0.1,1,0.663688,0.663688,0.663688,
0.1,0.1,0.1,0.1,0.1,0.663688,1,1,1,
0.1,0.1,0.1,0.1,0.1,0.663688,1,1,1,
0.1,0.1,0.1,0.1,0.1,0.663688,1,1,1),ncol=9,
            dimnames=
            list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
                 c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))

## for simulresults_sigma10_200mmp5perc.csv
w <- matrix(c(1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
0.05,1,0.05,0.137993,0.05,0.05,0.05,0.05,0.05,
0.05,0.05,1,0.05,0.177113,0.05,0.05,0.05,0.05,
0.05,0.137993,0.05,1,0.05,0.05,0.05,0.05,0.05,
0.05,0.05,0.177113,0.05,1,0.05,0.05,0.05,0.05,
0.05,0.05,0.05,0.05,0.05,1,0.194024,0.194024,0.194024,
0.05,0.05,0.05,0.05,0.05,0.194024,1,1,1,
0.05,0.05,0.05,0.05,0.05,0.194024,1,1,1,
0.05,0.05,0.05,0.05,0.05,0.194024,1,1,1),ncol=9,
            dimnames=
            list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
                 c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))

## 100 0.01
w <- matrix(c(1,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,
0.01,1,0.01,0.371474,0.01,0.01,0.01,0.01,0.01,
0.01,0.01,1,0.01,0.420848,0.0459244,0.01,0.01,0.01,
0.01,0.371474,0.01,1,0.01,0.01,0.01,0.01,0.01,
0.01,0.01,0.420848,0.01,1,0.01,0.01,0.01,0.01,
0.01,0.01,0.0459244,0.01,0.01,1,0.440482,0.440482,0.440482,
0.01,0.01,0.01,0.01,0.01,0.440482,1,1,1,
0.01,0.01,0.01,0.01,0.01,0.440482,1,1,1,
0.01,0.01,0.01,0.01,0.01,0.440482,1,1,1 ),ncol=9,
dimnames=
    list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
         c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))


image(w,xaxt="n",yaxt="n",col=grey.colors(100,start=0,end=1),
      breaks=0:100/100,
      main="Mating probabilities (white: yes; absolute black: no; grey=maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[1+23*(w[i,j]>0.5)])}}


image(w,xaxt="n",yaxt="n",col=rgb(99:0/99,99:0/99,1),
      breaks=0:100/100,
      main="Mating probabilities (white: no; blue: yes; light blue: maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[24-23*(w[i,j]>0.5)])}}

## 300 0.05
w <- matrix(c(1,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
0.05,1,0.05,0.0512609,0.05,0.05,0.05,0.05,0.05,
0.05,0.05,1,0.05,0.0745375,0.05,0.05,0.05,0.05,
0.05,0.0512609,0.05,1,0.05,0.05,0.05,0.05,0.05,
0.05,0.05,0.0745375,0.05,1,0.05,0.05,0.05,0.05,
0.05,0.05,0.05,0.05,0.05,1,0.085464,0.085464,0.085464,
0.05,0.05,0.05,0.05,0.05,0.085464,1,1,1,
0.05,0.05,0.05,0.05,0.05,0.085464,1,1,1,
0.05,0.05,0.05,0.05,0.05,0.085464,1,1,1),ncol=9,
dimnames=
    list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
         c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))


image(w,xaxt="n",yaxt="n",col=grey.colors(100,start=0,end=1),
      breaks=0:100/100,
      main="Mating probabilities (white: yes; absolute black: no; grey=maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[1+23*(w[i,j]>0.5)])}}

## 300 0.001
w <- matrix(c(1,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
0.001,1,0.001,0.0512609,0.001,0.001,0.001,0.001,0.001,
0.001,0.001,1,0.001,0.0745375,0.001,0.001,0.001,0.001,
0.001,0.0512609,0.001,1,0.001,0.001,0.001,0.001,0.001,
0.001,0.001,0.0745375,0.001,1,0.001,0.001,0.001,0.001,
0.001,0.001,0.001,0.001,0.001,1,0.085464,0.085464,0.085464,
0.001,0.001,0.001,0.001,0.001,0.085464,1,1,1,
0.001,0.001,0.001,0.001,0.001,0.085464,1,1,1,
0.001,0.001,0.001,0.001,0.001,0.085464,1,1,1),ncol=9,
dimnames=
    list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
         c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))


image(w,xaxt="n",yaxt="n",col=grey.colors(100,start=0,end=1),
      breaks=0:100/100,
      main="Mating probabilities (white: yes; absolute black: no; grey=maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],3),
                                   col=colors()[1+23*(w[i,j]>0.5)])}}

image(w,xaxt="n",yaxt="n",col=rgb(99:0/99,99:0/99,1),
      breaks=0:100/100,
      main="Mating probabilities (white: no; blue: yes; light blue: maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[24-23*(w[i,j]>0.5)])}}

## 5.439683 0.02339205 0.6187571 best_nosel_ (no selection)
w <- matrix(c(1,0.466612,0.0533909,0.663084,0.10706,0.0233921,0.0233921,0.0233921,0.0233921,
0.466612,1,0.494904,0.947557,0.679418,0.210633,0.118911,0.118911,0.118911,
0.0533909,0.494904,1,0.317725,0.954011,0.845706,0.68039,0.68039,0.68039,
0.663084,0.947557,0.317725,1,0.482403,0.111823,0.0572338,0.0572338,0.0572338,
0.10706,0.679418,0.954011,0.482403,1,0.675493,0.49586,0.49586,0.49586,
0.0233921,0.210633,0.845706,0.111823,0.675493,1,0.956381,0.956381,0.956381,
0.0233921,0.118911,0.68039,0.0572338,0.49586,0.956381,1,1,1,
0.0233921,0.118911,0.68039,0.0572338,0.49586,0.956381,1,1,1,
0.0233921,0.118911,0.68039,0.0572338,0.49586,0.956381,1,1,1),ncol=9,
dimnames=
    list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
         c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))


image(w,xaxt="n",yaxt="n",col=grey.colors(100,start=0,end=1),
      breaks=0:100/100,
      main="Mating probabilities (white: yes; absolute black: no; grey=maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[1+23*(w[i,j]>0.5)])}}


## 5.527374 0.02368335 best2_nosel_ (no selection)
w <- matrix(c(
1,0.460913,0.0509276,0.658707,0.103272,0.0236833,0.0236833,0.0236833,0.0236833,
0.460913,1,0.489324,0.946735,0.675198,0.20541,0.114898,0.114898,0.114898,
0.0509276,0.489324,1,0.311907,0.953288,0.843424,0.67618,0.67618,0.67618,
0.658707,0.946735,0.311907,1,0.476767,0.107942,0.0546544,0.0546544,0.0546544,
0.103272,0.675198,0.953288,0.476767,1,0.671235,0.490284,0.490284,0.490284,
0.0236833,0.20541,0.843424,0.107942,0.671235,1,0.955693,0.955693,0.955693,
0.0236833,0.114898,0.67618,0.0546544,0.490284,0.955693,1,1,1,
0.0236833,0.114898,0.67618,0.0546544,0.490284,0.955693,1,1,1,
0.0236833,0.114898,0.67618,0.0546544,0.490284,0.955693,1,1,1),ncol=9,
dimnames=
    list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
         c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))


image(w,xaxt="n",yaxt="n",col=grey.colors(100,start=0,end=1),
      breaks=0:100/100,
      main="Mating probabilities (white: yes; absolute black: no; grey=maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[1+23*(w[i,j]>0.5)])}}

image(w,xaxt="n",yaxt="n",col=rgb(99:0/99,99:0/99,1),
      breaks=0:100/100,
      main="Mating probabilities (white: no; blue: yes; light blue: maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[24-23*(w[i,j]>0.5)])}}

## best with additive genomic architecture (and selection):
w <- matrix(c(
1,0.987,0.987,0.987,0.987,0.987,0.987,0.987,0.987,
0.987,1,0.987,0.987,0.987,0.987,0.987,0.987,0.987,
0.987,0.987,1,0.987,0.987,0.987,0.987,0.987,0.987,
0.987,0.987,0.987,1,0.987,0.987,0.987,0.987,0.987,
0.987,0.987,0.987,0.987,1,0.987,0.987,0.987,0.987,
0.987,0.987,0.987,0.987,0.987,1,0.987,0.987,0.987,
0.987,0.987,0.987,0.987,0.987,0.987,1,0.987,0.987,
0.987,0.987,0.987,0.987,0.987,0.987,0.987,1,0.987,
0.987,0.987,0.987,0.987,0.987,0.987,0.987,0.987,1 
),ncol=9,
            dimnames=
            list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
                 c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))


image(w,xaxt="n",yaxt="n",col=grey.colors(100,start=0,end=1),
      breaks=0:100/100,
      main="Mating probabilities (white: yes; absolute black: no; grey=maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[1+23*(w[i,j]>0.5)])}}

## 2.181321 0.612858 0.7717959 -175.619 (selection model)
w <- matrix(c(1,0.736632,0.612858,0.848103,0.612858,0.612858,0.612858,0.612858,0.612858,
              0.736632,1,0.754228,0.978631,0.85642,0.612858,0.612858,0.612858,0.612858,
              0.612858,0.754228,1,0.631425,0.981298,0.935007,0.856911,0.856911,0.856911,
              0.848103,0.978631,0.631425,1,0.746529,0.612858,0.612858,0.612858,0.612858,
              0.612858,0.85642,0.981298,0.746529,1,0.854433,0.754811,0.754811,0.754811,
              0.612858,0.612858,0.935007,0.612858,0.854433,1,0.982275,0.982275,0.982275,
              0.612858,0.612858,0.856911,0.612858,0.754811,0.982275,1,1,1,
              0.612858,0.612858,0.856911,0.612858,0.754811,0.982275,1,1,1,
              0.612858,0.612858,0.856911,0.612858,0.754811,0.982275,1,1,1),ncol=9,
            dimnames=
            list(c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"),
                 c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd")))


image(w,xaxt="n",yaxt="n",col=grey.colors(100,start=0,end=1),
      breaks=0:100/100,
      main="Mating probabilities (white: yes; absolute black: no; grey=maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[1+23*(w[i,j]>0.5)])}}

image(w,xaxt="n",yaxt="n",col=rgb(99:0/99,99:0/99,1),
      breaks=0:100/100,
      main="Mating probabilities (white: no; blue: yes; light blue: maybe)")
axis(1,0:8/8,dimnames(w)[[1]])
axis(2,0:8/8,dimnames(w)[[2]])
for(i in 1:9) {for(j in 1:9) {text((i-1)/8,(j-1)/8,round(w[i,j],2),
                                   col=colors()[24-23*(w[i,j]>0.5)])}}

## Dispersal probabilities:
## disp <- c(0.00441961,0.00800287,0.0136133,0.0217541,0.0326567,0.0460534,0.0610108,0.0759292,0.0887702,0.097495,0.10059,0.097495,0.0887702,0.0759292,0.0610108,0.0460534,0.0326567,0.0217541,0.0136133,0.00800287,0.00441961)

## for simulresults_sigma10_10000.csv and simulresults_sigma10_mmp10perc.csv:
disp <- c(7.4336e-07,7.99187e-06,6.69151e-05,0.000436341,0.00221592,0.00876415,0.0269955,0.0647588,0.120985,0.176033,0.199471,0.176033,0.120985,0.0647588,0.0269955,0.00876415,0.00221592,0.000436341,6.69151e-05,7.99187e-06,7.4336e-07)

## for *_sigma5_*
disp <- c(7.6946e-23,1.02798e-18,5.05227e-15,9.13472e-12,6.07588e-09,1.48672e-06,0.00013383,0.00443185,0.053991,0.241971,0.398942,0.241971,0.053991,0.00443185,0.00013383,1.48672e-06,6.07588e-09,9.13472e-12,5.05227e-15,1.02798e-18,7.6946e-23)

disp <- disp/sum(disp)
barplot(names=-10:10*5,height=disp,
     main="East-West dispersal probability",
     xlab="Kilometers [only East-West direction counts]",ylab="Probability")


d <- read.csv("../c++/simulresults.csv")

dim(d)
x <- (1:200-99.5)*5
par(mfcol=c(1,1))
## par(mfcol=c(4,1))
## for(n in c(seq(1,dim(d)[1],1),dim(d)[1])) {
## for(n in 1:dim(d)[1]) {
for(n in which(d[1]==2000)){
## for(n in c(2,20,100,dim(d)[1])) {
    plot(x,unlist(d[n,2:201]),xlab="West     <---  Distance from initial contact line [km]  --->     East",
         ##main=paste("After",d[n,1],"generations"),
         ylab="Frequency of dark allele",t="l",col="red",ylim=c(0,1),
         ## xlim=c(-650,650)
         )
    abline(h=c(0,1))
    abline(v=0,lty=2)
    legend("topright",col=c("red","blue"),lty=c(1,2),legend=c("Locus 1","Locus 2"),bg="white")
    lines(x,unlist(d[n,202:401]),col="blue",lty=2)
}

hz <- 642
axis(3,-4:8*100-242,0:12*100)

hz <- 430
axis(3,-4:4*100-30,0:8*100)

hz <- 302
axis(3,-4:4*100+98,0:8*100)

hz <- 500
axis(3,-4:6*100-100,0:10*100)

km <- sapply(split(dat$DistanceAlongTransect,dat$NestID),mean)-hz
##km[km < -500] <- -500
##km[km >  500] <-  500
points(jitter(km,10),jitter(sapply(split(dat$NewHybridChr18,dat$NestID),mean)/2),col=rgb(1,0,0,0.5),pch=1)
points(jitter(km,10),jitter(sapply(split(dat$NDP,dat$NestID),mean)/2),col=rgb(0,0,1,0.5),pch=3)


## calculating shift of hybrid zones; that is infection point

inflex <- function(y) {
    x <- which.min(y[3:200]-y[1:198]) + 1
    z <- lm(y[(x-2):(x+2)]~I((x-2):(x+2))+I(((x-2):(x+2))^2)+I(((x-2):(x+2))^3))$coef
    (-z[[3]]/3/z[[4]]-100.5)*5
}

d <- read.csv("../c++/best_sel_simres.csv",h=F)
(which.min(d[200,4:201]-d[200,2:199])-100.5)*5
inflex(unlist(d[200,2:201]))
## 189.5637
(which.min(d[200,203:401]-d[200,201:399])-100.5)*5
inflex(unlist(d[200,202:401]))
## 188.2985

d <- read.csv("../c++/best2_nosel_simres.csv",h=F)
(which.min(d[200,4:201]-d[200,2:199])-100.5)*5
inflex(unlist(d[200,2:201]))
## 10.9898
(which.min(d[200,203:401]-d[200,201:399])-100.5)*5
inflex(unlist(d[200,202:401]))
## 7.189441

d <- read.csv("../c++/sel_add_bl_simres.csv",h=F)
(which.min(d[200,4:201]-d[200,2:199])-100.5)*5
inflex(unlist(d[200,2:201]))
## -0.1190476
(which.min(d[200,203:401]-d[200,201:399])-100.5)*5
inflex(unlist(d[200,202:401]))
## -0.1190476

d <- read.csv("../SimResults/pure/0.645_0.318_1.0_pure_2000simres.csv",h=F)
(which.min(d[200,4:201]-d[200,2:199])-100.5)*5
inflex(unlist(d[200,2:201]))
## -142.0535
(which.min(d[200,203:401]-d[200,201:399])-100.5)*5
inflex(unlist(d[200,202:401]))
## -141.695

## (end of calculation of inflection points)

x <- (1:200-99.5)*5
simres <- system("ls ../SimResults/*simres.csv ../SimResults/*/*simres.csv",intern=TRUE)
for(sr in simres) {
    d <- read.csv(sr)
    n <- dim(d)[1]
    plot(x,unlist(d[n,2:201]),xlab="West     <---  Distance from initial contact line [km]  --->     East",
         main=paste(sr," after",d[n,1],"generations"),ylab="Frequency of dark allele",t="l",col="red",ylim=c(0,1))
    abline(h=c(0,1))
    abline(v=0,lty=2)
    legend("bottomleft",col=c("red","blue"),lty=c(1,2),legend=c("Locus 1","Locus 2"),bg="white")
    lines(x,unlist(d[n,202:401]),col="blue",lty=2)
}


n <- 199
plot(x,unlist(d[n,2:201]),xlab="West     <---  Distance from initial contact line [km]  --->     East",
     main=paste("After",d[n,1],"generations\n"),ylab="Frequency of dark allele",t="l",col="red",ylim=c(0,1))
abline(h=c(0,1))
abline(v=0,lty=2)
legend("topright",col=c("red","blue"),lty=c(1,2),legend=c("Locus 1","Locus 2"),bg="white")
lines(x,unlist(d[n,202:401]),col="blue",lty=2)

dat <- read.table("../Data/data_central_red.txt",h=TRUE)
s <- sapply(split(dat[,5:6],round(dat$DistanceAlongTransect/50)),
            function(x){
                s <- rep(0,9)
                for(i in 1:nrow(x)) {
                    s[x[i,1]*3+x[i,2]+1] = s[x[i,1]*3+x[i,2]+1]+1
                }
                s
            }
            )
s <- sweep(s,2,apply(s,2,sum),"/")
barplot(s,width=c(6,1,1,1,1,1,6),space=0,col=grey(0.8-0.8*c(0.0, 0.3743380, 0.7339315, 0.2748254, 0.6409001, 0.9094524,
          1.0000000, 1.0000000, 1.0000000)),xaxt="n",border=rgb(0.4,0.4,0.6),
        legend.text=c("LLll","LLdl","LLdd","DLll","DLdl","DLdd","DDll","DDdl","DDdd"),args.legend=c(bg="white"))
s


## calculate frequencies under HWE and LD=0 assumptions:
af <- cbind(unlist(d[n,2:201]),unlist(d[n,202:401]))
gfreq <- cbind((1-af[,1])^2*(1-af[,2])^2,(1-af[,1])^2*2*(1-af[,2])*af[,2],(1-af[,1])^2*af[,2]^2,
               2*af[,1]*(1-af[,1])*(1-af[,2])^2,2*af[,1]*(1-af[,1])*2*(1-af[,2])*af[,2],2*af[,1]*(1-af[,1])*af[,2]^2,
               af[,1]^2*(1-af[,2])^2,af[,1]^2*2*(1-af[,2])*af[,2],af[,1]^2*af[,2]^2)
barplot(t(gfreq),col=grey(0.8-0.8*c(0.0, 0.3743380, 0.7339315, 0.2748254, 0.6409001, 0.9094524,
                 1.0000000, 1.0000000, 1.0000000)),xaxt="n",border=rgb(0.4,0.4,0.6))

gf <- read.csv("../c++/best_sel_finfreqs.csv",h=F)
barplot(t(gf[,10:2]),
        col=grey(0.8-0.8*c(1.0000000, 1.0000000, 1.0000000, 0.9094524,
        0.6409001, 0.2748254, 0.7339315, 0.3743380, 0.0)),
        xaxt="n",border=rgb(0.4,0.4,0.6))

barplot(sweep(as.array(t(gf[,2:10])),2,apply(as.array(t(gf[,2:10])),2,mean),"/"),col=grey(0.8-0.8*c(0.0, 0.3743380, 0.7339315, 0.2748254, 0.6409001, 0.9094524,1.0000000, 1.0000000, 1.0000000)),xaxt="n",border=rgb(0.4,0.4,0.6))

finfr <- system("ls ../SimResults/*finfreqs.csv ../SimResults/*/*finfreqs.csv",intern=TRUE)
barplot(t(t(rep(1,9))),col=grey(0.8-0.8*c(0.0, 0.3743380, 0.7339315, 0.2748254, 0.6409001, 0.9094524,
                           1.0000000, 1.0000000, 1.0000000)),yaxt="n",border=rgb(0.4,0.4,0.6))
text(0.7,1:9-0.5,col="white",c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"))

for(ffr in finfr) {
    gf <- read.csv(ffr,h=F)
    barplot(t(gf[,2:10]),col=grey(0.8-0.8*c(0.0, 0.3743380, 0.7339315, 0.2748254, 0.6409001, 0.9094524,
                             1.0000000, 1.0000000, 1.0000000)),xaxt="n",border=rgb(0.4,0.4,0.6),main=ffr)
}

yellowblue=c(rgb(1,1,1),rgb(1,1,0.5),rgb(0.8,1,0),rgb(0.7,0.1,1),rgb(0.7,0.1,0.5),rgb(0.7,0.1,0),rgb(0,0,1),rgb(0,0,0.5),rgb(0,0,0))
barplot(t(t(rep(1,9))),col=yellowblue,yaxt="n")
text(0.7,1:9-0.5,col=c(rep("black",4),rep("white",5)),c("LLll","LLld","LLdd","LDll","LDld","LDdd","DDll","DDld","DDdd"))

for(ffr in finfr) {
    gf <- read.csv(ffr,h=F)
    barplot(t(gf[,2:10]),col=yellowblue,xaxt="n",main=ffr)
}




par(mfrow=c(2,2))
plot(c(),xlim=c(300,600),ylim=c(-350,-200),ylab="log likelihood",xlab="initial contact km from western point",
     main="initially both alles fixed in their subpop")
## plot(c(),xlim=c(300,800),ylim=c(-2500,-200))
cn <- 1
for(str in c("sbd_5_finfreqs.csv sbd 0.05",
            "100_0_finfreqs.csv 100 0.0",
             "100_01_finfreqs.csv 100 0.001",
             "100_1_finfreqs.csv 100 0.01",
             "100_3_finfreqs.csv 100 0.03",
             "200_0_finfreqs.csv 200 0.0",
             "200_01_finfreqs.csv 200 0.001",
             "200_1_finfreqs.csv 200 0.01",
             "200_3_finfreqs.csv 200 0.03",
             "300_0_finfreqs.csv 300 0.0",
             "300_01_finfreqs.csv 300 0.001",
             "300_1_finfreqs.csv 300 0.01",
             "300_3_finfreqs.csv 300 0.03",
             "400_0_finfreqs.csv 400 0.0",
             "400_01_finfreqs.csv 400 0.001",
             "400_1_finfreqs.csv 400 0.01",
             "400_3_finfreqs.csv 400 0.03",
             "300_5_finfreqs.csv 300 0.05",
             "300_10_finfreqs.csv 300 0.1",
             "300_20_finfreqs.csv 300 0.2") ) {
    sout <- system(paste("../c++/genotype_llh ../SimResults/",str," ../Data/data_central_red.txt",sep=""),intern=TRUE)
    print(sout[50])
    v <- strsplit(sout[1],'\t')[[1]]
    for(i in 2:length(sout)) {
        w <- strsplit(sout[i],'\t')[[1]]
        if(length(v)>1) {
            vn <- as.numeric(v)
            wn <- as.numeric(w)
            lines(c(vn[1],wn[1]),c(vn[2],wn[2]),col=cn,lty=round(cn/7)+1)
        }
        v <- w
    }
    cn <- cn+1
}
legend("topright",lty=round(1:20/7)+1,col=1:20, legend=c("sbd 5", "100 0.0","100 0.001","100 0.01","100 0.03",
                                  "200 0.0","200 0.001","200 0.01","200 0.03",
                                  "300 0.0","300 0.001","300 0.01","300 0.03",
                                  "400 0.0","400 0.001","400 0.01","400 0.03",
                                  "300 0.05","300 0.1","300 0.2"),bg="white")


plot(c(),xlim=c(300,600),ylim=c(-350,-200),ylab="log likelihood",xlab="initial contact km from western point",
     main="initially loc2 95% dark in western subpop (loc1 fixed)")
## plot(c(),xlim=c(300,800),ylim=c(-2500,-200))
cn <- 1
for(str in  c("100_0_finfreqs.csv 100 0.0",
             "100_01_finfreqs.csv 100 0.001",
             "100_1_finfreqs.csv 100 0.01",
             "100_3_finfreqs.csv 100 0.03",
             "200_0_finfreqs.csv 200 0.0",
             "200_01_finfreqs.csv 200 0.001",
             "200_1_finfreqs.csv 200 0.01",
             "200_3_finfreqs.csv 200 0.03",
             "300_0_finfreqs.csv 300 0.0",
             "300_01_finfreqs.csv 300 0.001",
             "300_1_finfreqs.csv 300 0.01",
             "300_3_finfreqs.csv 300 0.03",
             "400_0_finfreqs.csv 400 0.0",
             "400_01_finfreqs.csv 400 0.001",
             "400_1_finfreqs.csv 400 0.01",
             "400_3_finfreqs.csv 400 0.03",
             "300_5_finfreqs.csv 300 0.05",
             "300_10_finfreqs.csv 300 0.1",
             "300_20_finfreqs.csv 300 0.2") ) {
    sout <- system(paste("../c++/genotype_llh ../SimResults/init_2nd_loc_het/",str," ../Data/data_central_red.txt",sep=""),intern=TRUE)
    print(sout[50])
    v <- strsplit(sout[1],'\t')[[1]]
    for(i in 2:length(sout)) {
        w <- strsplit(sout[i],'\t')[[1]]
        if(length(v)>1) {
            vn <- as.numeric(v)
            wn <- as.numeric(w)
            lines(c(vn[1],wn[1]),c(vn[2],wn[2]),col=cn,lty=round(cn/7)+1)
        }
        v <- w
    }
    cn <- cn+1
}
legend("topright",lty=round(1:19/7)+1,col=1:19, legend=c("100 0.0","100 0.001","100 0.01","100 0.03",
                                                    "200 0.0","200 0.001","200 0.01","200 0.03",
                                                    "300 0.0","300 0.001","300 0.01","300 0.03",
                                                    "400 0.0","400 0.001","400 0.01","400 0.03",
                                                    "300 0.05","300 0.1","300 0.2"),bg="white")

plot(c(),xlim=c(300,600),ylim=c(-350,-200),ylab="log likelihood",xlab="initial contact km from western point",
     main="initially loc2 80% dark in western subpop (loc1 fixed)")
## plot(c(),xlim=c(300,800),ylim=c(-2500,-200))
cn <- 1
for(str in  c(    "10000_01_finfreqs.csv 100 0.001",
             "10000_1_finfreqs.csv 100 0.01",
             "10000_3_finfreqs.csv 100 0.03",
             "10000_5_finfreqs.csv 100 0.01",
             "10000_10_finfreqs.csv 100 0.03",
              "10000_20_finfreqs.csv 100 0.03",
              "300_01_finfreqs.csv 300 0.001",
             "300_1_finfreqs.csv 300 0.01",
             "300_3_finfreqs.csv 300 0.03",
             "300_5_finfreqs.csv 300 0.05",
             "300_10_finfreqs.csv 300 0.1",
              "300_20_finfreqs.csv 300 0.2") ) {
    sout <- system(paste("../c++/genotype_llh ../SimResults/init_2nd_loc_het80/",str," ../Data/data_central_red.txt",sep=""),intern=TRUE)
    print(sout[50])
    v <- strsplit(sout[1],'\t')[[1]]
    for(i in 2:length(sout)) {
        w <- strsplit(sout[i],'\t')[[1]]
        if(length(v)>1) {
            vn <- as.numeric(v)
            wn <- as.numeric(w)
            lines(c(vn[1],wn[1]),c(vn[2],wn[2]),col=cn,lty=round(cn/7)+1)
        }
        v <- w
    }
    cn <- cn+1
}
legend("topright",lty=round(1:12/7)+1,col=1:12, legend=c("10000 0.001","10000 0.01","10000 0.03",
                                                    "10000 0.05","10000 0.1","10000 0.2",
                                                    "300 0.001","300 0.01","300 0.03",
                                                         "300 0.05","300 0.1","300 0.2"),bg="white")


plot(c(),xlim=c(300,600),ylim=c(-350,-200),ylab="log likelihood",xlab="initial contact km from western point",
     main="initially loc2 60% dark in western subpop (loc1 fixed)")
## plot(c(),xlim=c(300,800),ylim=c(-2500,-200))
cn <- 1
for(str in  c("sbd_5_finfreqs.csv sbd 0.05",
                "10000_01_finfreqs.csv 100 0.001",
             "10000_1_finfreqs.csv 100 0.01",
             "10000_3_finfreqs.csv 100 0.03",
             "10000_5_finfreqs.csv 100 0.01",
             "10000_10_finfreqs.csv 100 0.03",
              "10000_20_finfreqs.csv 100 0.03",
              "300_01_finfreqs.csv 300 0.001",
             "300_1_finfreqs.csv 300 0.01",
             "300_3_finfreqs.csv 300 0.03",
             "300_5_finfreqs.csv 300 0.05",
             "300_10_finfreqs.csv 300 0.1",
              "300_20_finfreqs.csv 300 0.2") ) {
    sout <- system(paste("../c++/genotype_llh ../SimResults/init_2nd_loc_het60/",str," ../Data/data_central_red.txt",sep=""),intern=TRUE)
    print(sout[50])
    v <- strsplit(sout[1],'\t')[[1]]
    for(i in 2:length(sout)) {
        w <- strsplit(sout[i],'\t')[[1]]
        if(length(v)>1) {
            vn <- as.numeric(v)
            wn <- as.numeric(w)
            lines(c(vn[1],wn[1]),c(vn[2],wn[2]),col=cn,lty=round(cn/7)+1)
        }
        v <- w
    }
    cn <- cn+1
}
legend("topright",lty=round(1:13/7)+1,col=1:13, legend=c("sbd 5","10000 0.001","10000 0.01","10000 0.03",
                                                    "10000 0.05","10000 0.1","10000 0.2",
                                                    "300 0.001","300 0.01","300 0.03",
                                                         "300 0.05","300 0.1","300 0.2"),bg="white")

## dev.copy2pdf(file="likelihood_curves.pdf")

### three locus simulations

d <- read.csv("../c++/threeloci/best_sel_simres.csv")
d0 <- read.csv("../c++/threeloci/allmate_simres.csv")
x <- (1:200-99.5)*5
par(mfcol=c(2,1))
for(n in 1:dim(d0)[1]) {
    plot(x,unlist(d[n,2:201]),xlab="West     <---  Distance from initial contact line [km]  --->     East",
         main=paste("After",d[n,1],"generations; locus 3 neutral; all unlinked"),
         ylab="Frequency of dark allele",t="l",col="red",ylim=c(0,1))
    abline(h=c(0,1))
    abline(v=0,lty=2)
    legend("topright",col=c("red","blue","black"),lty=c(1,2,3),lwd=c(1,1,4),legend=c("Locus 1","Locus 2","Locus 3"),bg="white")
    lines(x,unlist(d[n,202:401]),col="blue",lty=2)
    lines(x,unlist(d[n,402:601]),col="black",lwd=4,lty=2)
    
    plot(x,unlist(d0[n,2:201]),xlab="West     <---  Distance from initial contact line [km]  --->     East",
         main=paste("After",d0[n,1],"generations; all three loci neutral and unlinked"),ylab="Frequency of dark allele",t="l",col="red",ylim=c(0,1))
    abline(h=c(0,1))
    abline(v=0,lty=2)
    legend("topright",col=c("red","blue","black"),lty=c(1,2,3),lwd=c(1,1,4),legend=c("Locus 1","Locus 2","Locus 3"),bg="white")
    lines(x,unlist(d0[n,202:401]),col="blue",lty=2)
    lines(x,unlist(d0[n,402:601]),col="black",lwd=4,lty=2)
}

## dev.copy2pdf(file="threelocicomparison.pdf")

par(mfrow=c(4,2))
for(n in c(1,49,99,199)) {
       plot(x,unlist(d[n,2:201]),xlab="West     <---  Distance from initial contact line [km]  --->     East",
         main=paste("After",d[n,1],"generations"),
         ylab="Frequency of dark allele",t="l",col="red",ylim=c(0,1))
    abline(h=c(0,1))
    abline(v=0,lty=2)
    # legend("topright",col=c("red","blue","black"),lty=c(1,2,3),lwd=c(1,1,4),legend=c("Locus 1","Locus 2","Locus 3"),bg="white")
    lines(x,unlist(d[n,202:401]),col="blue",lty=2)
    lines(x,unlist(d[n,402:601]),col="black",lwd=4,lty=2)
    
    plot(x,unlist(d0[n,2:201]),xlab="West     <---  Distance from initial contact line [km]  --->     East",
         main=paste("After",d0[n,1],"generations"),ylab="",t="l",col="red",ylim=c(0,1))
    abline(h=c(0,1))
    abline(v=0,lty=2)
    # legend("topright",col=c("red","blue","black"),lty=c(1,2,3),lwd=c(1,1,4),legend=c("Locus 1","Locus 2","Locus 3"),bg="white")
    lines(x,unlist(d0[n,202:401]),col="blue",lty=2)
    lines(x,unlist(d0[n,402:601]),col="black",lwd=4,lty=2)
}
## dev.copy2pdf(file="threelocithroughtime.pdf")

