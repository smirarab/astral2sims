source("./setup.r")

######### SETUP ############
colr=c("#ffb380"  ,"#217821" ,"#5f8dd3", "#e76bf3","#782121")
rename<-list( 
        "ASTRAL-472-fast" ="astral-v47-p1-new" ,"ASTRAL-472-slow" = "astral-v47-p2-new" , 
        "ASTRAL-473-b-fast" ="astral-v47b-p1" ,"ASTRAL-473-b-slow" = "astral-v47b-p2",
        "ASTRAL-473-c-fast" ="astral-v47c-p1" ,"ASTRAL-473-c-slow" = "astral-v47c-p2",
        "ASTRAL-473-fast" ="astral-v473-p1" ,"ASTRAL-473-slow" = "astral-v473-p2",
	"astral-1"="astral1","ASTRAL-I"="astral311", 
        "ASTRAL-p0+true"="astral-v47-p0-true","ASTRAL-p0+true"="astral-v474-p0-true",
        "ASTRAL-p0"="astral-v47-p0-new","ASTRAL-p0"="astral-v474-p0",
        "ASTRAL-II" ="astral-v476-p1","ASTRAL-II" ="astral-v474-p1" ,"ASTRAL-474-slow" = "astral-v474-p2",
        "ASTRAL-II + true st"="astral-v474-p1-true","ASTRAL-II + true st"="astral-v476-p1-true",
        "ASTRAL-II (true gt)"="astral-v474-p1-truegenetrees",
        "ASTRAL-II + random resolution of polytomies"="astral-v474-p1-fullresolved",
        "NJst"="njst","MRL"="mrl","Recursive-greedy" = "rgreedy", "Greedy"= "greedy","MP-EST"="mpest",
	"CA-ML"="concatenatedtree")
contl <- list("No-contraction"="-1","10%"="10","33%"="33","50%"="50")
astralvars <- c("ASTRAL-II + true st","ASTRAL-II","ASTRAL-I")
methods <- c("CA-ML","NJst","ASTRAL-II","MP-EST")
methodcons <- c("NJst","ASTRAL-II")

#############################
draw <-function (sc,name,xfactor,xlab,ylab,ffactor=sc$genes,v="fn",tm=theme(legend.position="bottom")){
    dt <- sc; dt$x <- xfactor; dt$f <- ffactor
    a=dcast(dt,x+rep+f~method,value.var=v,drop=TRUE);
    dt=melt(data=a,id.vars=c("x","rep","f"));
    dc <- summarySE(dt,measurevar="value",groupvars=c("x","f","variable"),na.rm=F,conf.interval=.95)
    #print(dc)
    print(head(dc))
    sp=0.7
    ba1 <- ggplot(aes(x=x,y=value,fill=variable),data=dc)+xlab(xlab)+ylab(ylab)+geom_bar(position=position_dodge(sp),width=sp,color=1,stat="identity")+theme_bw()+scale_fill_discrete(name="")+geom_errorbar(position=position_dodge(sp), size=.5, width=.25, aes(ymin=value-se, ymax=value+se))+facet_wrap(~f,ncol=3)+tm
    ba2 <- ggplot(aes(x=x,y=value,fill=variable),data=dc)+xlab(xlab)+ylab(ylab)+geom_bar(position=position_dodge(sp),width=sp,color=1,stat="identity")+theme_bw()+scale_fill_discrete(name="")+coord_cartesian(ylim=c(-0.02,.15))+geom_errorbar(position=position_dodge(sp), size=.5, width=.25, aes(ymin=value-se, ymax=value+se))+facet_wrap(~f,ncol=3)+tm
    bl <- qplot(x,value,geom=c("point","line"),data=dc,colour=variable,shape=variable,group=variable,ylab=(ylab),xlab=xlab)+theme_bw()+scale_shape_manual(name="",values=c(19,23,15,2,1,0,4:14))+scale_colour_discrete(name="")+geom_errorbar( width=.1,size=.2, aes(ymin=value-se, ymax=value+se))+facet_wrap(~f,ncol=3)+tm
    bx <- ggplot(aes(x=x,y=value,fill=variable),data=dt)+xlab(xlab)+ylab(ylab)+geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),position=position_dodge(.8))+theme_bw()+scale_fill_discrete(name="")+facet_wrap(~f,ncol=3)+tm
    ec <- qplot(value,data=dt,stat="ecdf",geom="step",color=variable)+xlab(ylab)+ylab(xlab)+theme_bw()+facet_grid(f~x)+tm
    v <- qplot(x,value,fill=variable,data=dt,geom="violin")+xlab("Tree length (controls ILS)")+ylab(ylab)+theme_bw()+scale_fill_discrete(name="")+facet_wrap(~f,ncol=3)+tm
    d <- qplot(value,geom="density",data=dt,colour=variable)+xlab(ylab)+theme_bw()+scale_fill_discrete(name="")+facet_grid(f~x)+tm
    pdf(name,width=9,height=6)
    print(bl);
    print(ba1); #print(ba2);
    print(ec);
     print(bx); print(d); print(v);
    dev.off()
}

draw2 <-function (sc,name,xlab,ylab,v="fn",thm=theme(legend.position="bottom")){
}

############### load data ###################
if (TRUE) {
sc=read.csv("scores.stat",sep=" ",header=F)
names(sc) <- c("taxa","height","rate","rep","method","genes","contract","all","fnc","fn","res","fpc","fp")
levels(sc$method)<-rename
sc$contract <- as.factor(sc$contract)
levels(sc$contract) <- contl
sc$inds = "1";sc[grep("-2",sc$taxa),"inds"] = "2";sc[grep("-5",sc$taxa),"inds"] = "5"
levels(sc$taxa )<-list("50"="50","100"="100", "200"="200-2", "200"="200-5", "200"="200","500"="500", "1000"="1000")
sc$rate<-factor(sc$rate,levels=c(1e-06,1e-07))
sc$rf<-(sc$fn+sc$fp)/2
#sc$genes<-as.factor(paste(sc$genes,"genes"))
sc$genes<-as.factor(sc$genes)

gt=read.csv("gt-err.stat",sep=" ",header=F);
gt$V4<-as.factor(gt$V4);gt$V3<-factor(gt$V3,levels=c(1e-06, 1e-07));gt$V2<-as.factor(gt$V2);
gt$rf=(gt$V8+gt$V11)/2
gt.sum = summarySE(gt,measurevar="rf",groupvars=c("V1","V2","V3","V4"),na.rm=T,conf.interval=.95)
names(gt.sum)<-c(names(sc)[1:4],"N","gtrf",names(gt.sum)[7:9])
sc=merge(sc,gt.sum[,c(1:4,6)])
sc$gterr="medium";sc[sc$gtrf>0.4,"gterr"]="high"; sc[sc$gtrf<=0.25,"gterr"]="low";
sc$gterr <- factor(sc$gterr, levels=c("low","medium","high"))

scd1 = sc[sc$taxa == 200 & sc$inds ==1 & sc$contract == "No-contraction",]
scd2 = sc[sc$height == "2M" & sc$inds ==1 & sc$contract == "No-contraction" & sc$rate == 1e-06 ,]
scm = scd1[scd1$method %in% methods,]
scs = scd2[scd2$method %in% methods,]
scv = scd1[scd1$method %in% astralvars ,]
scvs = scd2[scd2$method %in% astralvars,]
scvc = sc[sc$taxa == 200 & sc$inds ==1  & sc$method == "ASTRAL-II" ,]
sciv = sc[sc$taxa == 200 & sc$contract == "No-contraction" & sc$height == "500K" & sc$rate == 1e-06 & sc$method %in% astralvars ,]

rn=read.csv("runtime.stat",sep=" ",header=F);
names(rn)<-c(names(sc)[1:7],"time","junk")
rn$contract <- as.factor(rn$contract)
levels(rn$contract) <- contl
levels(rn$method)<-rename
rn$inds = "1";rn[grep("-2",rn$taxa),"inds"] = "2";rn[grep("-5",rn$taxa),"inds"] = "5"
levels(rn$taxa )<-list("50"="50","100"="100", "200"="200-2", "200"="200-5", "200"="200","500"="500", "1000"="1000")
#rn$genes<-as.factor(paste(rn$genes,"genes"))
rn$genes<-as.factor(rn$genes)
rn$rate<-factor(rn$rate,levels=c(1e-06,1e-07))
rn$time <- rn$time/3600
rnv <- rn[rn$taxa == 200 & rn$inds ==1 & rn$method  %in% astralvars & rn$contract == "No-contraction",c(1:6,8)]
rnvc <- rn[rn$taxa == 200 & rn$inds ==1 & rn$method  == "ASTRAL-II",c(1:7,8)]
rnm <- rn[rn$taxa == 200 & rn$inds ==1 & rn$method %in% methods & rn$contract == "No-contraction",c(1:6,8)]
rns <- rn[rn$height == "2M" & rn$inds ==1 & rn$rate == 1e-06 &  rn$method %in% methods & rn$contract == "No-contraction",c(1:6,8)];
rnvs = rn[rn$height == "2M" & rn$inds ==1 & rn$rate == 1e-06 & rn$method %in% astralvars & rn$contract == "No-contraction",c(1:6,8)]

}

################## Draw "full" figures
#draw(scv,"figs/st-astral_variant-ILS.pdf", with(scv,interaction(rate,height,sep="/")),"model condition (length/rate)","Species tree error (FN)",tm=th2)
#draw(scvs,"figs/st-astral_variant-size.pdf", scvs$taxa,"number of taxa","Species tree error (FN)",tm=th2)
#draw(scm,"figs/st-all_methods-ILS.pdf", with(scm,interaction(rate,height,sep="/")),"model condition (length/rate)","Species tree error (FN)",tm=th2)
##draw(scm,"figs/st-all_methods-ILS-fp.pdf", with(scm,interaction(rate,height,sep="/")),"model condition (length/rate)","Species tree error (FP)","fp",tm=th2)
##draw(scm,"figs/st-all_methods-ILS-RF.pdf", with(scm,interaction(rate,height,sep="/")),"model condition (length/rate)","Species tree error (RF)","rf",tm=th2)
#draw(scs,"figs/st-all_methods-size.pdf",scs$taxa,"number of taxa","Species tree error (FN)")
#draw(sciv,"figs/st-astral_variant-individuals.pdf",sciv$inds,"number of individuals","Species tree error (FN)")

############### Draw "full" running time ###################
#draw(rnv,"figs/rt-astral_variant-ILS.pdf", with(rnv,interaction(rate,height,sep="/")),xlab="model condition (length/rate)",ylab="Running time (hours)",v="time",tm=th3)
#draw(rnvs,"figs/rt-astral_variant-size.pdf", rnvs$taxa ,"number of taxa","Running time (hours)",v="time",tm=th2)
#draw(rnm,"figs/rt-all_methods-ILS.pdf", with(rnm,interaction(rate,height,sep="/")),"model condition (length/rate)","Running time (hours)",v="time",tm=th2)
#draw(rns,"figs/rt-all_methods-size.pdf",rns$taxa,"number of taxa","Running time (hours)",v="time",tm=th1)


################### Run statistical tests
dtb=scv;an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II","ASTRAL-I") ,c(1:6,10)]));nrow(an); summary(aov(fn ~ method*(height+rate+genes),data=an))
dtb=scvc;an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II") ,c(1:7,10)]));nrow(an); summary(aov(fn ~ contract*(height+rate+genes),data=an))
dtb=scvc;an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II") & contract %in% c("No-contraction","33%") ,c(1:7,10)]));nrow(an); summary(aov(fn ~ contract*(height+rate+genes),data=an))
dtb=rbind(scv,scvs);an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II","ASTRAL-II + true st") ,c(1:6,10)]));nrow(an); summary(aov(fn ~ method*(taxa+rate+height+genes),data=an))
dtb=scs;an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II","NJst") ,c(1:6,10)]));nrow(an); summary(aov(fn ~ method*(taxa+genes),data=an))
dtb=scs;an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II","CA-ML") ,c(1:6,10)]));nrow(an); summary(aov(fn ~ method*(taxa+genes),data=an))
dtb=scm;an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II","CA-ML") ,c(1:6,10)]));nrow(an); summary(aov(fn ~ method*(height+rate+genes),data=an))
dtb=scm;an=droplevels(with(dtb,dtb[method %in% c("ASTRAL-II","NJst") ,c(1:6,10)]));nrow(an); summary(aov(fn ~ method*(height+rate+genes),data=an))

##################### draw pretty figures ####################################
thm=theme(legend.position="bottom")
thm2=theme(legend.position=c(.85,.2))
thm3=theme(legend.position="none")
xlab="genes"

pdf("figs/methods-fn.pdf",width=8,height=4.8)
 dt=melt(data=dcast(scm,genes+rep+height+rate~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuPu")+facet_grid(rate~height)+thm3
 dt=melt(data=dcast(scs,genes+rep+taxa~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","taxa"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuPu")+facet_wrap(~taxa,nrow=2)+thm2
dev.off()

pdf("figs/methods-rt.pdf",width=8,height=4.8)
 dt=melt(data=dcast(rnm,genes+rep+height+rate~method,value.var="time",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Running time (hours)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuPu")+facet_grid(rate~height)+thm
 dt=melt(data=dcast(rns,genes+rep+taxa~method,value.var="time",drop=TRUE),id.vars=c("genes","rep","taxa"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Running time (hours)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuPu")+facet_wrap(~taxa,scales="free",nrow=2)+thm2
dev.off()
pdf("figs/methods-size-rt.pdf",width=4,height=4)
   qplot(taxa,value,stat="summary",fun.y=mean,data=dt,color=variable,linetype=variable,geom=c("point","line"),group=interaction(genes,variable),shape=genes)+
   xlab("number of taxa")+ylab("Running time (hours)")+theme_bw()+scale_color_brewer(name="method",palette="Dark2")+scale_linetype_discrete(name="method")+
   theme(legend.position=c(.4,.8),legend.box.just=0,legend.direction=1)
dev.off()



pdf("figs/methods-gtcategories.pdf",width=7,height=9)
 dt=melt(data=dcast(scm,genes+rep+height+rate+gterr~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","rate","height","gterr"));
    ggplot(aes(x=gterr,y=value,fill=variable),data=dt)+xlab("gene tree error")+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuPu")+facet_grid(height+rate~genes)+thm
 dt=melt(data=dcast(scs,genes+rep+taxa+gterr~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","taxa","gterr"));
    ggplot(aes(x=gterr,y=value,fill=variable),data=dt)+xlab("gene tree error")+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuPu")+facet_grid(taxa~genes)+thm
dev.off()


pdf("figs/astral-contract.pdf",width=8,height=4.5)
 dt=melt(data=dcast(scvc,genes+rep+height+rate~contract,value.var="fn",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="PuRd")+facet_grid(rate~height)+thm3
 a=dcast(scvc,genes+rep+height+rate~contract,value.var="fn",drop=TRUE);a$`10%`=a$`10%`-a$`No-contraction`;a$`33%`=a$`33%`-a$`No-contraction`;
 a$`50%`=a$`50%`-a$`No-contraction`;a$`No-contraction`=0;
 dt=melt(a,id.vars=c("genes","rep","rate","height"))
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Delta FN due to branch contraction")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="PuRd")+facet_grid(rate~height)+thm3
 dt=melt(data=dcast(rnvc,genes+rep+height+rate~contract,value.var="time",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Running time (hours)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="PuRd")+facet_grid(rate~height)+thm
dev.off()

pdf("figs/astral-variants.pdf",width=8,height=4)
 dt=melt(data=dcast(scv,genes+rep+height+rate~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuGn")+facet_grid(rate~height)+thm
 dt=melt(data=dcast(rnv,genes+rep+height+rate~method,value.var="time",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Running time (hours)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuGn")+facet_grid(rate~height)+thm
 dt=melt(data=dcast(scvs,genes+rep+taxa~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","taxa"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuGn")+facet_wrap(~taxa,scales="free",nrow=2)+thm2
 dt=melt(data=dcast(rnvs,genes+rep+taxa~method,value.var="time",drop=TRUE),id.vars=c("genes","rep","taxa"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Running time (hours)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuGn")+facet_wrap(~taxa,scales="free",nrow=2)+thm2
dev.off()
pdf("figs/variants-size-rt.pdf",width=4,height=4)
   qplot(taxa,value,stat="summary",fun.y=mean,data=dt,color=variable,linetype=variable,geom=c("point","line"),group=interaction(genes,variable),shape=genes)+
   xlab("number of taxa")+ylab("Running time (hours)")+theme_bw()+scale_color_brewer(name="method",palette="Dark2")+scale_linetype_discrete(name="method")+
   theme(legend.position=c(.2,.7),legend.box.just=0)
dev.off()

pdf("figs/astral-rand-res.pdf",width=8,height=4)
 dt=melt(data=dcast(scd1[scd1$method %in% c("ASTRAL-II + random resolution of polytomies","ASTRAL-II"),],genes+rep+height+rate~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="BuGn")+facet_grid(rate~height)+thm
dev.off()

pdf("figs/astral-truegt.pdf",width=8,height=4)
 dt=melt(data=dcast(scd1[scd1$method %in% c("ASTRAL-II","ASTRAL-II (true gt)","CA-ML"),],genes+rep+height+rate~method,value.var="fn",drop=TRUE),id.vars=c("genes","rep","rate","height"));
    ggplot(aes(x=genes,y=value,fill=variable),data=dt)+xlab(xlab)+ylab("Species tree topological error (FN)")+
    geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),colour="black",position=position_dodge(.8))+
    theme_bw()+scale_fill_brewer(name="",palette="OrRd")+facet_grid(rate~height)+thm
dev.off()


lm_eqn = function(df){
    m = lm(rf ~ gtrf, df);
    eq <- substitute(italic(y) == a + b ~ italic(x)*","~~italic(r)^2~"="~r2,
         list(a = format(coef(m)[1], digits = 1),
              b = format(coef(m)[2], digits = 2),
             r2 = format(summary(m)$r.squared, digits = 2)))
    as.character(as.expression(eq));
}

pdf("figs/gt-vs-st.pdf",width=10,height=6)
 eq <- ddply(scm[scm$method=="ASTRAL-II",],genes~height+rate,lm_eqn)
    ggplot(aes(x=gtrf,y=rf),data=scm[scm$method=="ASTRAL-II",])+theme_bw()+th1+geom_point(size=2,alpha=0.5)+geom_smooth(method=lm)+
    facet_grid(genes~height+rate)+geom_text(data=eq,aes(x = .4, y = .38,label=V1),size=2.5, parse = TRUE, inherit.aes=FALSE)+
    xlab("true versus estimated gene trees (normalized RF distance)")+ylab("Species tree error (normalized FN)")+ggtitle("ASTRAL-II")
 eq <- ddply(scm[scm$method=="NJst",],genes~height+rate,lm_eqn)
    ggplot(aes(x=gtrf,y=rf),data=scm[scm$method=="NJst",])+theme_bw()+th1+geom_point(size=2,alpha=0.5)+geom_smooth(method=lm)+
    facet_grid(genes~height+rate)+geom_text(data=eq,aes(x = .4, y = .38,label=V1),size=2.5, parse = TRUE, inherit.aes=FALSE)+
    xlab("true versus estimated gene trees (normalized RF distance)")+ylab("Species tree error (normalized FN)")+ggtitle("NJst")
 eq <- ddply(scm[scm$method=="CA-ML",],genes~height+rate,lm_eqn)
    ggplot(aes(x=gtrf,y=rf),data=scm[scm$method=="CA-ML",])+theme_bw()+th1+geom_point(size=2,alpha=0.5)+geom_smooth(method=lm)+
    facet_grid(genes~height+rate)+geom_text(data=eq,aes(x = .4, y = .38,label=V1),size=2.5, parse = TRUE, inherit.aes=FALSE)+
    xlab("true versus estimated gene trees (normalized RF distance)")+ylab("Species tree error (normalized FN)")+ggtitle("CA-ML")
dev.off()
pdf("figs/gt-vs-st-mpest.pdf",width=8,height=5)
 eq <- ddply(scs[scs$method=="MP-EST",],genes~taxa,lm_eqn)
    ggplot(aes(x=gtrf,y=rf),data=scs[scs$method=="MP-EST",])+theme_bw()+th1+geom_point(size=2,alpha=0.5)+geom_smooth(method=lm)+
    facet_grid(taxa~genes)+geom_text(data=eq,aes(x = .3, y = .42,label=V1),size=3.5, parse = TRUE, inherit.aes=FALSE)+
    xlab("true versus estimated gene trees (normalized RF distance)")+ylab("Species tree error (normalized FN)")+ggtitle("MP-EST")
dev.off()

scm$fn=scm$fn*100; dc <- summarySE(scm,measurevar="fn",groupvars=c("genes","rate","height","method"),na.rm=F,conf.interval=.95);dc$sum=paste(format(dc$fn),format(dc$se),sep="$\\pm$");l<-latex(recast(rate+height+genes~method,measure.var="sum",data=dc),rowname=NULL,file="manuscript/scm.tex",caption="Species tree error for 200 taxa and varying levels of ILS and number of genes. Average and standard errors are given for FN rates over 50 replicates, with the exception of 500K, 1e-6, where 3 replicates had to be removed and 47 replicates are used.",label="scm")

scs$fn=scs$fn*100; dc <- summarySE(scs,measurevar="fn",groupvars=c("taxa","genes","method"),na.rm=F,conf.interval=.95);dc$sum=paste(format(dc$fn),format(dc$se),sep="$\\pm$");l<-latex(recast(taxa+genes~method,measure.var="sum",data=dc),rowname=NULL,file="manuscript/scs.tex",caption="Species tree error for varying number of taxa with 200M/1e-6 tree shape and varying number of genes. Average and standard errors are given for FN rates over 50 replicates, with the exceptions of 50 and 100 taxa, where 3 and 2 replicates had to be removed and 47 and 48 replicates are used.",label="scs")

tmp=recast(scv,genes+rate+height+rep~method, measure.var="fn"); tmp$delta=100*(tmp$`ASTRAL-I`-tmp$`ASTRAL-II`); dc <- summarySE(tmp,measurevar="delta",groupvars=c("genes","rate","height"),na.rm=F,conf.interval=.95); dc$sum=paste(format(dc$delta,digits=1),format(dc$se,digits=1),sep="$\\pm$");l<-latex(recast(genes~interaction(height,rate,sep="/"),measure.var="sum",data=dc),rowname=NULL,file="manuscript/astral-vars.tex",caption="Delta species tree error between ASTRAL-I and ASTRAL-II for dataset I. Average and standard errors are given for the difference between FN rates of ASTRAL-II and ASTRAL-I.",label="tab:varils")

tmp=recast(scvs,genes+taxa+rep~method, measure.var="fn"); tmp$delta=100*(tmp$`ASTRAL-I`-tmp$`ASTRAL-II`); dc <- summarySE(tmp,measurevar="delta",groupvars=c("genes","taxa"),na.rm=F,conf.interval=.95); dc$sum=paste(format(dc$delta,digits=1),format(dc$se,digits=1),sep="$\\pm$");l<-latex(recast(genes~taxa,measure.var="sum",data=dc),rowname=NULL,file="manuscript/astral-vars-size.tex",caption="Delta species tree error between ASTRAL-I and ASTRAL-II for dataset II. Average and standard errors are given for the difference between FN rates of ASTRAL-II and ASTRAL-I.",label="tab:varsize")
