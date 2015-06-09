source("./setup.r")


if (TRUE){
gt=read.csv("gt-err.stat",sep=" ",header=F);
gt$V4<-as.factor(gt$V4);gt$V3<-factor(gt$V3,levels=c(1e-06, 1e-07));gt$V2<-as.factor(gt$V2);
gt$rf=(gt$V8+gt$V11)/2
gt.sum = summarySE(gt,measurevar="rf",groupvars=c("V1","V2","V3","V4"),na.rm=T,conf.interval=.95)
}

ils=read.csv("gt-truedist.stat",sep=" ",header=F)
ils$V3=factor(ils$V3,levels=c(1e-06, 1e-07))
pdf("figs/simulated-ils.pdf",width=5,height=3)
qplot(V7,data=ils[ils$V1=="200",],geom="density",linetype=V3,color=V2,xlab="RF distance (true species tree vs true gene trees)")+scale_linetype_discrete(name="rate")+scale_color_brewer(name="tree height",palette="Dark2")+theme_bw()+theme(legend.position=c(.6,.87),legend.direction=1,legend.box="horizontal")+scale_x_continuous(labels=percent)+ggtitle("(a) True gene tree discordance")
qplot(rf,data=gt[gt$V1=="200",],geom="density",color=V2,xlab="RF distance (true vs estimated)",linetype=V3)+theme_bw()+scale_linetype_discrete(name="rate")+theme(legend.position=c(.75,.9))+scale_x_continuous(labels=percent)+theme(legend.position="none")+scale_color_brewer(name="tree height",palette="Dark2")+
ggtitle("(b) Gene tree estimation error")
dev.off()

pdf("figs/ils-size.pdf",width=3,height=6)
qplot(V7,data=ils[ils$V1=="200",],geom="density",adjust=2,linetype=V3,xlab="RF distance\n(true species tree vs true gene trees)")+scale_linetype_discrete(name="rate")+theme_bw()+theme(legend.position=c(.65,.92),legend.direction=1,legend.box="horizontal")+scale_x_continuous(labels=percent)+ggtitle("(a) Varying tree shape")+facet_wrap(~V2,ncol=1,scales="free_y")
qplot(V7,data=ils[ils$V2=="2M" & ils$V3==1e-06,],geom="density",adjust=2,xlab="RF distance\n(true species tree vs true gene trees)")+theme_bw()+scale_x_continuous(labels=percent)+ggtitle("(b) Varying number of taxa")+facet_wrap(~V1,ncol=1,scales="free_y")
dev.off()

pdf("figs/gt-err-sum.pdf",width=3,height=6)
qplot(rf,data=gt[gt$V1=="200",],geom="density",main="Varying tree shape",xlab="RF distance (true vs estimated)",linetype=V3)+theme_bw()+scale_linetype_discrete(name="rate")+facet_wrap(~V2,ncol=1)+theme(legend.position=c(.75,.9))+scale_x_continuous(labels=percent)+theme(plot.margin=unit(c(2,0,2,0),"mm"))
qplot(rf,data=gt[gt$V2=="2M" & gt$V3 == 1e-06,],geom="density",main="Varying number of taxa", xlab="RF distance (true vs estimated)")+theme_bw()+facet_wrap(~V1,ncol=1)+theme(legend.position="bottom")+scale_x_continuous(labels=percent)+theme(plot.margin=unit(c(2,0,2,0),"mm"))
#dc<-summarySE(scm,measurevar="fn",groupvars=c("genes","rate","height","method"),na.rm=F,conf.interval=.95);dc$sum=paste(format(dc$fn),format(dc$se),sep="$\\pm$");latex(recast(rate+height+genes~method,measure.var="sum",data=dc),rowname=NULL,file="manuscript/scm.tex",caption="Species tree error for 200 taxa and varying levels of ILS and number of genes. Average and standard errors are given for FN rates over 50 replicates, with the exception of 500K, 1e-6, where 3 replicates had to be removed and 47 replicates are used.",label="scm")
dev.off()

pdf("figs/gt-err-rep.pdf",width=10,height=3)
gt.sum$rep=interaction(gt.sum$V1,gt.sum$V2,gt.sum$V3,gt.sum$V4)
qplot(reorder(rep,rf),rf,data=with(gt.sum,gt.sum[V1==200,]),ylab="RF (true vs estimated gene trees)",xlab="replicates ordered by gene tree error",main="gene tree error for replicates of various model conditions with 200 taxa")+facet_wrap(~V1+V2+V3,scales="free_x",ncol=6)+geom_errorbar(aes(ymin=rf-sd, ymax=rf+sd),color="blue")+theme_bw()+th4+theme(plot.margin=unit(c(0,0,0,0),"mm"))
qplot(reorder(rep,rf),rf,data=with(gt.sum,gt.sum[V2=="2M" & V3 == 1e-06,]),ylab="RF (true vs estimated gene trees)",xlab="replicates ordered by gene tree error",main="gene tree error for replicates of various model conditions with 2M/1e-06 tree shape")+facet_wrap(~V1+V2+V3,scales="free_x",ncol=6)+geom_errorbar(aes(ymin=rf-sd, ymax=rf+sd),color="blue")+theme_bw()+th4+theme(plot.margin=unit(c(0,0,0,0),"mm"))
dev.off()

q()

pdf("figs/gt-err.pdf",width=9)
#gt=gt[gt$V4<26,]
qplot(V4,V8,data=gt,geom="violin",ylab="Gene tree estimation error (FN)",xlab="Replicates")+facet_grid(V2~V3)+theme_bw()
qplot(V4,V8,data=gt,geom="boxplot",ylab="Gene tree estimation error (FN)",xlab="Replicates",outlier.shape=1)+facet_grid(V2~V3)+theme_bw()
qplot(V8,data=gt,geom="density",xlab="Gene tree estimation error (FN)",color=V2,linetype=V3)+theme_bw()+scale_color_discrete(name="tree length")+scale_linetype_discrete(name="speciation rate")+facet_wrap(~V1)
qplot(V8,data=gt,geom="density",xlab="Gene tree estimation error (FN)",linetype=V3)+theme_bw()+scale_linetype_discrete(name="speciation rate")+facet_wrap(~V2,ncol=1)+theme(legend.position="bottom")+facet_wrap(~V1)

qplot(V4,V11,data=gt,geom="violin",ylab="Gene tree estimation error (FP)",xlab="Replicates")+facet_grid(V2~V3)+theme_bw()
qplot(V4,V11,data=gt,geom="boxplot",ylab="Gene tree estimation error (FP)",xlab="Replicates",outlier.shape=1)+facet_grid(V2~V3)+theme_bw()
qplot(V11,data=gt,geom="density",xlab="Gene tree estimation error (FP)",color=V2,linetype=V3)+theme_bw()+scale_color_discrete(name="tree length")+scale_linetype_discrete(name="speciation rate")
qplot(V11,data=gt,geom="density",xlab="Gene tree estimation error (FP)",linetype=V3)+theme_bw()+scale_linetype_discrete(name="speciation rate")+facet_wrap(~V2,ncol=1)+theme(legend.position="bottom")
dev.off()
