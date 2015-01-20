# Read data files
b1=read.csv("base.song",sep=" ",header=F)
b2=read.csv("bases.1kp",sep=" ",header=F)
b3=read.csv("bases.avian",sep=" ",header=F)
r1=read.csv("rates.song",sep=" ",header=F)
r2=read.csv("rates.1kp",sep=" ",header=F)
r3=read.csv("rates.avian",sep=" ",header=F)

# Sample 400 rows from each dataset
s1=sample(nrow(b1),400)
s2=sample(nrow(b2),400)
s3e=sample(which(b3$V7=="exon" & 5000>b3$V9& b3$V9>1000),400)
s3i=sample(which(b3$V7=="intron" & 5000>b3$V9 & b3$V9>1000),400)
s3u=sample(which(b3$V7=="uce" & 5000>b3$V9 & b3$V9>1000),400)
b1=b1[s1,3:6]
b2=b2[s2,3:6]
b3=rbind(b3[s3e,3:6], b3[s3i,3:6], b3[s3u,3:6])
cols=c(2,14,12,15,10,13,11)
r1=r1[s1,cols]
r2=r2[s2,cols]
r3=rbind(r3[s3e,cols], r3[s3i,cols], r3[s3u,cols])

require(sirt)

print ("########## estimating base frequencies")
round(dirichlet.mle(b1)$alpha)
round(dirichlet.mle(b2)$alpha)
round(dirichlet.mle(b3)$alpha)
round(dirichlet.mle(rbind(b1,b2,b3))$alpha)

print ("########## estimating gtr matrices")
round(dirichlet.mle(r1[,2:7])$alpha)
round(dirichlet.mle(r2[,2:7])$alpha)
round(dirichlet.mle(r3[,2:7])$alpha)
round(dirichlet.mle(rbind(r1[,2:7],r2[,2:7],r3[,2:7]))$alpha)

print ("########## estimating alhpa shape parameter")
1/mean(r1$V2)
1/mean(r2$V2)
1/mean(r3$V2)
1/mean(rbind(r1$V2,r2$V2,r3$V2))
