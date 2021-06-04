## ran ../../c++/finite_sim/nest_grid_sim > fs.txt and repeated for fs2.txt, fs3.txt fs4.txt
## and ../../c++/finite_sim/nest_grid_imprint_sim > fsi.txt and repeated for fsi2.txt, fsi3.txt fsi4.txt
## and ../../c++/finite_sim/nest_grid_imprint_f_sim > fsif.txt and repeated for fsif2.txt
## and ../../c++/finite_sim/nest_grid_imprint_m_sim > fsim.txt and repeated for fsim2.txt

d <- read.table("fs.txt.gz",sep=",")
par(mfcol=c(3,1))
image(t(as.matrix(d[1:500,])),col=grey(10:0/11),axes=FALSE,main="Crow nest mean phenotype in finite population model 1 generation after initial contact",
      xlab="Distance to initial contact zone [km]")
axis(1,c(0,0.250,0.500,0.750,1),c(-500,-250,0,250,500))
image(t(as.matrix(d[2001:2500,])),col=grey(10:0/11),axes=FALSE,main="Crow nest mean phenotype in finite population model 2000 generation after initial contact",
      xlab="Distance to initial contact zone [km]")
axis(1,c(0,0.250,0.500,0.750,1),c(-500,-250,0,250,500))
image(t(as.matrix(d[3001:3500,])),col=grey(10:0/11),axes=FALSE,main="Crow nest mean phenotype in finite population model 5000 generation after initial contact",
      xlab="Distance to initial contact zone [km]")
axis(1,c(0,0.250,0.500,0.750,1),c(-500,-250,0,250,500))
## dev.copy(file="finiteimage.png",dev=png,width=15,height=30,units="cm",res=500); dev.off()

image(t((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10)))

par(mfcol=c(4,1))
for(f in c("fs.txt.gz", "fs2.txt.gz","fs3.txt.gz","fs4.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Self-reference model after 2000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4,col="blue")
}

for(f in c("fs.txt.gz", "fs2.txt.gz","fs3.txt.gz","fs4.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Reference matching model after 5000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    lines(-500:499+0.5,apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4,col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4,col="blue")
}

par(mfcol=c(3,1))
d <- read.table("fs.txt.gz",sep=",")
plot(c(), main="Finite-populations self-reference model 1 generation after initial contact",
     panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
     xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1))
lines(-500:499+0.5,apply((as.matrix(d[501:1000,]) %/% 1000) + ((as.matrix(d[501:1000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[501:1000,]) %/% 100) %% 10) + (as.matrix(d[501:1000,]) %% 10),2,mean)/4,col="blue")
 plot(c(), main="After 2000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
      xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1))
lines(-500:499+0.5,apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4,col="blue")
plot(c(), main="After 5000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
     xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1))
lines(-500:499+0.5,apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4,col="red")
lines(-500:499+0.5,apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4,col="blue")
## dev.copy2pdf(file="finitereference.pdf")

clines2000 <- data.frame(pos=-500:499+0.5)
clines2000$sr1 <- apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4
clines2000$sr2 <- apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4

d <- read.table("fsi.txt.gz",sep=",")
image(t(as.matrix(d[1:500,])),col=grey(10:0/10),breaks=0:11/10-0.05)
image(t(as.matrix(d[1001:1500,])),col=grey(10:0/10),breaks=0:11/10-0.05)
image(t(as.matrix(d[1501:2000,])),col=grey(10:0/10),breaks=0:11/10-0.05)
image(t(as.matrix(d[2001:2500,])),col=grey(10:0/10),breaks=0:11/10-0.05)
image(t(as.matrix(d[3001:3500,])),col=grey(10:0/10),breaks=0:11/10-0.05)

par(mfcol=c(3,1))
d <- read.table("fsi.txt.gz",sep=",")
plot(c(), main="Mid-parental imprinting model 1 generation after initial contact",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
     xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
lines(-500:499+0.5,apply((as.matrix(d[501:1000,]) %/% 1000) + ((as.matrix(d[501:1000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[501:1000,]) %/% 100) %% 10) + (as.matrix(d[501:1000,]) %% 10),2,mean)/4,col="blue")
plot(c(), main="After 1000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
lines(-500:499+0.5,apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4,col="blue")
plot(c(), main="After 2000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
     xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
lines(-500:499+0.5,apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4,col="blue")
## dev.copy2pdf(file="midparentclines.pdf")

clines2000$mp1 <- apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4
clines2000$mp2 <- apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4



par(mfcol=c(5,1))
for(f in c("fsi.txt.gz", "fsi2.txt.gz","fsi3.txt.gz","fsi4.txt.gz","fsi5.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Imprinting model 1 generation after initial contact",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[501:1000,]) %/% 1000) + ((as.matrix(d[501:1000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[501:1000,]) %/% 100) %% 10) + (as.matrix(d[501:1000,]) %% 10),2,mean)/4,col="blue")
}

par(mfcol=c(5,1))
for(f in c("fsi.txt.gz", "fsi2.txt.gz","fsi3.txt.gz","fsi4.txt.gz","fsi5.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Imprinting model after 1000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4,col="blue")
}

par(mfcol=c(5,1))
for(f in c("fsi.txt.gz", "fsi2.txt.gz","fsi3.txt.gz","fsi4.txt.gz","fsi5.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Imprinting model after 2000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4,col="blue")
}

par(mfcol=c(3,1))
d <- read.table("fsif.txt.gz",sep=",")
plot(c(), main="Father-imprinting model 1 generation after initial contact",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
     xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
lines(-500:499+0.5,apply((as.matrix(d[501:1000,]) %/% 1000) + ((as.matrix(d[501:1000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[501:1000,]) %/% 100) %% 10) + (as.matrix(d[501:1000,]) %% 10),2,mean)/4,col="blue")
plot(c(), main="After 2000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
     xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
lines(-500:499+0.5,apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4,col="blue")
plot(c(), main="After 5000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
     xlab="Distance from initial contact line [km]",ylab="Allele frequency")
abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
lines(-500:499+0.5,apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
lines(-500:499+0.5,apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4,col="blue")
## dev.copy2pdf(file="fatherimprintclines.pdf")


clines2000$fi1 <- apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4
clines2000$fi2 <- apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4
## save(clines2000,file="clines2000.RData")

par(mfcol=c(3,1))
plot(clines2000$pos,clines2000$sr1,col="red",t="l",main="Self-referenced",ylim=c(0,1))
lines(clines2000$pos,clines2000$sr2,col="blue")
plot(clines2000$pos,clines2000$fi1,col="red",t="l",main="Imprinting on father",ylim=c(0,1))
lines(clines2000$pos,clines2000$fi2,col="blue")
plot(clines2000$pos,clines2000$mp1,col="red",t="l",main="Imprinting on mid-parental phenotype",ylim=c(0,1))
lines(clines2000$pos,clines2000$mp2,col="blue")


par(mfcol=c(4,1))
for(f in c("fsif.txt.gz", "fsif2.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Father-imprinting model after 2000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4,col="blue")
}
for(f in c("fsim.txt.gz","fsim2.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Mother-imprinting model after 2000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[2501:3000,]) %/% 1000) + ((as.matrix(d[2501:3000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[2501:3000,]) %/% 100) %% 10) + (as.matrix(d[2501:3000,]) %% 10),2,mean)/4,col="blue")
}
## dev.copy2pdf(file="mofaimprint2000.pdf")

par(mfcol=c(4,1))
for(f in c("fsif.txt.gz", "fsif2.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Father-imprinting model after 5000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4,col="blue")
}
for(f in c("fsim.txt.gz","fsim2.txt.gz")) {
    d <- read.table(f,sep=",")
    plot(c(), main="Mother-imprinting model after 5000 generations",panel.first=grid(2,2),ylim=c(0,1),xlim=c(-500,500),
         xlab="Distance from initial contact line [km]",ylab="Allele frequency")
    abline(h=c(0,1),ylab="Allele frequency",xlab="Distance to initial contact line")
    lines(-500:499+0.5,apply((as.matrix(d[3501:4000,]) %/% 1000) + ((as.matrix(d[3501:4000,]) %/% 10) %% 10),2,mean)/4,t="l",col="red")
    lines(-500:499+0.5,apply(((as.matrix(d[3501:4000,]) %/% 100) %% 10) + (as.matrix(d[3501:4000,]) %% 10),2,mean)/4,col="blue")
}
## dev.copy2pdf(file="mofaimprint5000.pdf")
