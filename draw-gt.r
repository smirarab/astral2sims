require(ggplot2)
gt=read.csv("gt-err.stat",sep=" ",header=F);
gt$V3<-as.factor(gt$V3);gt$V2<-as.factor(gt$V2);gt$V1<-as.factor(gt$V1);
pdf("figs/gt-err.pdf")
qplot(V3,V6,data=gt,geom="violin",ylab="Gene tree estimation error",xlab="Replicates")+facet_grid(V1~V2)+theme_bw()
qplot(V3,V6,data=gt,geom="boxplot",ylab="Gene tree estimation error",xlab="Replicates",outlier.shape=1)+facet_grid(V1~V2)+theme_bw()
qplot(V6,data=gt,geom="density",xlab="Gene tree estimation error",color=V1,linetype=V2)+theme_bw()+scale_color_discrete(name="tree length")+scale_linetype_discrete(name="speciation rate")
qplot(V6,data=gt,geom="density",xlab="Gene tree estimation error",linetype=V2)+theme_bw()+scale_linetype_discrete(name="speciation rate")+facet_wrap(~V1,ncol=1)+theme(legend.position="bottom")
dev.off()
