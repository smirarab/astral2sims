library(scales)
require(reshape2)
require(ggplot2)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

######### SETUP ############
colr=c("#ffb380"  ,"#217821" ,"#5f8dd3", "#e76bf3","#782121")
rename<-list("ASTRAL-472-no-add"="astral-v47-p0-new", 
        "ASTRAL-472-fast" ="astral-v47-p1-new" ,"ASTRAL-472-slow" = "astral-v47-p2-new" , 
        "ASTRAL-473-b-fast" ="astral-v47b-p1" ,"ASTRAL-473-b-slow" = "astral-v47b-p2",
        "ASTRAL-473-c-fast" ="astral-v47c-p1" ,"ASTRAL-473-c-slow" = "astral-v47c-p2",
        "ASTRAL-473-fast" ="astral-v473-p1" ,"ASTRAL-473-slow" = "astral-v473-p2",
        "ASTRAL-474-default" ="astral-v474-p1" ,"ASTRAL-474-slow" = "astral-v474-p2",
        "ASTRAL-with-true"="astral-v47-p0-true","ASTRAL-with-true"="astral-v474-p0-true",
        "NJst"="njst","MRL"="mrl","Recursive-greedy" = "rgreedy", "Greedy"= "greedy")
methods <- c("Greedy","MRL","NJst","ASTRAL-474-default")
methodcons <- c("NJst","ASTRAL-474-default")
th1=theme(legend.position="bottom")
th2=theme(legend.position="bottom",axis.text.x=element_text(angle=45, hjust = 1))
th3=theme(legend.position="right",axis.text.x=element_text(angle=45, hjust = 1))

#############################
drawST <-function (sc,name,xfactor,xlab,ylab,v="V9",tm=theme(legend.position="bottom")){
    sc2 <- sc; sc2$x <- xfactor 
    a=dcast(sc2,x+V4+V6~V5,value.var=v,drop=TRUE);
    sc2=melt(data=a,id.vars=c("x","V4","V6"));
    dc <- summarySE(sc2,measurevar="value",groupvars=c("x","V6","variable"),na.rm=F,conf.interval=.95)
    print(dc)
    print(dc[,c("x","V6","variable","N")])
    print(head(dc))
    sp=0.7
    ba1 <- ggplot(aes(x=x,y=value,fill=variable),data=dc)+xlab(xlab)+ylab(ylab)+geom_bar(position=position_dodge(sp),width=sp,color=1,stat="identity")+theme_bw()+scale_fill_discrete(name="")+geom_errorbar(position=position_dodge(sp), size=.5, width=.25, aes(ymin=value-se, ymax=value+se))+facet_wrap(~V6,ncol=3)+tm
    ba2 <- ggplot(aes(x=x,y=value,fill=variable),data=dc)+xlab(xlab)+ylab(ylab)+geom_bar(position=position_dodge(sp),width=sp,color=1,stat="identity")+theme_bw()+scale_fill_discrete(name="")+coord_cartesian(ylim=c(-0.02,.15))+geom_errorbar(position=position_dodge(sp), size=.5, width=.25, aes(ymin=value-se, ymax=value+se))+facet_wrap(~V6,ncol=3)+tm
    bl <- qplot(x,value,geom=c("point","line"),data=dc,colour=variable,shape=variable,group=variable,ylab=(ylab),xlab=xlab)+theme_bw()+scale_shape_manual(name="",values=c(19,23,15,2,1,0,2,4:14))+scale_colour_discrete(name="")+geom_errorbar( width=.1,size=.2, aes(ymin=value-se, ymax=value+se))+facet_wrap(~V6,ncol=3)+tm
    bx <- ggplot(aes(x=x,y=value,fill=variable),data=sc2)+xlab(xlab)+ylab(ylab)+geom_boxplot(outlier.size=1.2,outlier.colour=rgb(0.01,.01,.01,.5),position=position_dodge(.8))+theme_bw()+scale_fill_discrete(name="")+facet_wrap(~V6,ncol=3)+tm
    v <- qplot(x,value,fill=variable,data=sc2,geom="violin")+xlab("Tree length (controls ILS)")+ylab(ylab)+theme_bw()+scale_fill_discrete(name="")+facet_wrap(~V6,ncol=3)+tm
    d <- qplot(value,geom="density",data=sc2,colour=variable)+xlab(ylab)+theme_bw()+scale_fill_discrete(name="")+facet_grid(V6~x)+tm
    pdf(name,width=9,height=6)
    print(ba1); #print(ba2);
    print(bl); print(bx); print(d); print(v);
    dev.off()
}


############### Draw species tree FN  ###################
sc=read.csv("scores.stat",sep=" ",header=F)
levels(sc$V5)<-rename
sc$V1<-factor(sc$V1)
sc$V3<-factor(sc$V3,levels=c(1e-06,1e-07))
sc$V13<-(sc$V9+sc$V12)/2
sc$V6<-as.factor(paste(sc$V6,"genes"))
scv = sc[sc$V1 == 200 & grepl("ASTRAL",sc$V5),]
scvs = sc[sc$V2 == "2M" & sc$V3 == 1e-06 & grepl("ASTRAL",sc$V5),]
scm = sc[sc$V1 == 200 & sc$V5 %in% methods,]
scs = sc[sc$V2 == "2M" & sc$V3 == 1e-06 &  sc$V5 %in% methods,]
scm2 = sc[sc$V1 == 200 & sc$V5 %in% methodcons,]
scs2 = sc[sc$V2 == "2M" & sc$V3 == 1e-06 &  sc$V5 %in% methodcons,]
drawST(scv,"figs/st-astral-variant-200.pdf", with(scv,interaction(V3,V2,sep="/")),"model condition (length/rate)","Species tree error (FN)",tm=th3)
drawST(scvs,"figs/st-astral-variant-size.pdf", scvs$V1,"number of taxa","Species tree error (FN)",tm=th3)
#drawST(scm,"figs/st-methods-200.pdf", with(scm,interaction(V3,V2,sep="/")),"model condition (length/rate)","Species tree error (FN)",tm=th2)
#drawST(scm,"figs/st-methods-200-fp.pdf", with(scm,interaction(V3,V2,sep="/")),"model condition (length/rate)","Species tree error (FP)","V12",tm=th2)
#drawST(scm,"figs/st-methods-200-RF.pdf", with(scm,interaction(V3,V2,sep="/")),"model condition (length/rate)","Species tree error (RF)","V13",tm=th2)
drawST(scs,"figs/st-size.pdf",scs$V1,"number of taxa","Species tree error (FN)")
drawST(scm2,"figs/st-methods-200-sc.pdf", with(scm2,interaction(V3,V2,sep="/")),"model condition (length/rate)","Species tree error (FN)",tm=th2)
drawST(scs2,"figs/st-size-sc.pdf",scs2$V1,"number of taxa","Species tree error (FN)")


############### Draw running time ###################
rn=read.csv("runtime.stat",sep=" ",header=F);
levels(rn$V5)<-rename
rn$V6<-as.factor(paste(rn$V6,"genes"))
rn$V3<-factor(rn$V3,levels=c(1e-06,1e-07))
rn$V1<-factor(rn$V1)
rn$V7 <- rn$V7/3600
rnv <- rn[rn$V1 == 200 & grepl("ASTRAL",rn$V5),]
rnm <- rn[rn$V1 == 200 & rn$V5 %in% methods,]
rns <- rn[rn$V2 == "2M" & rn$V3 == 1e-06 &  rn$V5 %in% methods,];

#drawRT(rnv,"figs/rt-astral-variant-200.pdf", with(rnv,interaction(V3,V2,sep="/")),"model condition (length/rate)")
#drawRT(rnm,"figs/rt-methods-200.pdf", with(rnm,interaction(V3,V2,sep="/")),"model condition (length/rate)")
#drawRT(rns,"figs/rt-size.pdf",rns$V1,"number of taxa")
drawST(rnv,"figs/rt-astral-variant-200.pdf", with(rnv,interaction(V3,V2,sep="/")),"model condition (length/rate)","Running time (hours)","V7",tm=th3)
drawST(rnm,"figs/rt-methods-200.pdf", with(rnm,interaction(V3,V2,sep="/")),"model condition (length/rate)","Running time (hours)","V7",tm=th2)
drawST(rns,"figs/rt-size.pdf",rns$V1,"number of taxa","Running time (hours)","V7",tm=th1)

