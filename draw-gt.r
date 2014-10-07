require(ggplot2)
gt=read.csv("gt-err.stat",sep=" ",header=F);
gt$V3<-as.factor(gt$V3);gt$V2<-as.factor(gt$V2);gt$V1<-as.factor(gt$V1);
qplot(V3,V6,data=gt,geom="violin",ylab="Gene tree estimation error",xlab="Replicates")+facet_grid(V1~V2)+theme_bw()
