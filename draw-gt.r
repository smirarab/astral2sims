require(ggplot2)
gt=read.csv("gt-err.stat",sep=" ",header=F);
gt=gt[gt$V4<26,]
gt$V4<-as.factor(gt$V4);gt$V3<-as.factor(gt$V3);gt$V2<-as.factor(gt$V2);
pdf("figs/gt-err.pdf",width=9)
qplot(V4,V8,data=gt,geom="violin",ylab="Gene tree estimation error (FN)",xlab="Replicates")+facet_grid(V2~V3)+theme_bw()
qplot(V4,V8,data=gt,geom="boxplot",ylab="Gene tree estimation error (FN)",xlab="Replicates",outlier.shape=1)+facet_grid(V2~V3)+theme_bw()
qplot(V8,data=gt,geom="density",xlab="Gene tree estimation error (FN)",color=V2,linetype=V3)+theme_bw()+scale_color_discrete(name="tree length")+scale_linetype_discrete(name="speciation rate")
qplot(V8,data=gt,geom="density",xlab="Gene tree estimation error (FN)",linetype=V3)+theme_bw()+scale_linetype_discrete(name="speciation rate")+facet_wrap(~V2,ncol=1)+theme(legend.position="bottom")

qplot(V4,V11,data=gt,geom="violin",ylab="Gene tree estimation error (FP)",xlab="Replicates")+facet_grid(V2~V3)+theme_bw()
qplot(V4,V11,data=gt,geom="boxplot",ylab="Gene tree estimation error (FP)",xlab="Replicates",outlier.shape=1)+facet_grid(V2~V3)+theme_bw()
qplot(V11,data=gt,geom="density",xlab="Gene tree estimation error (FP)",color=V2,linetype=V3)+theme_bw()+scale_color_discrete(name="tree length")+scale_linetype_discrete(name="speciation rate")
qplot(V11,data=gt,geom="density",xlab="Gene tree estimation error (FP)",linetype=V3)+theme_bw()+scale_linetype_discrete(name="speciation rate")+facet_wrap(~V2,ncol=1)+theme(legend.position="bottom")
dev.off()
