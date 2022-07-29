forestplot_v1=function(dat,show_95CI=T,zero = 1,boxsize = 0.4,lineheight =5,colgap =2,lwd.zero=2,lwd.ci=2
                       ,box_col='#458B00',summary_col="#8B008B",lines_col='black',zero_col='#7AC5CD'
                       ,xlab='HR',lwd.xaxis=2,lty.ci = "solid",graph.pos = 2,xlim=NULL,xlog=F){
  nc=ncol(dat)
  nr=nrow(dat)
  library(forestplot)
  col=fpColors(box=box_col,summary=summary_col,lines = lines_col,zero = zero_col)
  
  if(nc>3){
    hr=as.numeric(dat[,nc-2])
    lower=as.numeric(dat[,nc-1])
    upper=as.numeric(dat[,nc])
    
    if(is.null(xlim)){
      xlim=c(min(lower,na.rm = T),max(upper,na.rm = T))
    }
    if(is.infinite(max(xlim))){
      xlim=c(xlim[1],5)
    }
    if(is.infinite(min(xlim))){
      xlim=c(0,xlim[2])
    }
    
    if(min(xlim)<=0){
      xlog=F
    }
    smary=rep(F,length(hr))
    nind=which(is.na(lower)|is.na(upper)|is.na(hr))
    smary[nind]=T
    labeltext=as.matrix(dat[,1:(nc-3)])
    if(show_95CI){
      adt=paste0(round(hr,5),'(',round(lower,5),',',round(upper,5),')')
      adt[nind]=''
      labeltext=cbind(labeltext,adt)
      colnames(labeltext)=c(colnames(labeltext)[1:(ncol(labeltext)-1)],'Hazard Ratio(95% CI)')
    }
    if(graph.pos>ncol(labeltext)+1){
      labeltext=ncol(labeltext)+1
    }else if(graph.pos<2){
      graph.pos=2
    }
    hz_list=list('2'=gpar(lty=1,col=summary_col),
                 '3'=gpar(lty=1,col=summary_col)
    )
    names(hz_list)=c(2,nrow(labeltext)+2)
    p=forestplot(labeltext = rbind(colnames(labeltext),labeltext),
                 hrzl_lines = hz_list,
                 mean = c(NA,hr),
                 lower =c(NA,lower),
                 upper = c(NA,upper),
                 is.summary=c(T,smary),
                 zero = zero,
                 fn.ci_norm="fpDrawDiamondCI",
                 boxsize = boxsize, 
                 lineheight = unit(lineheight,'mm'),
                 colgap = unit(colgap,'mm'),
                 lwd.zero = lwd.zero,
                 lwd.ci = lwd.ci,
                 col=col,
                 xlab=xlab,
                 lwd.xaxis=lwd.xaxis,
                 lty.ci = lty.ci,
                 clip = xlim,
                 #xlog=xlog,
                 mar=unit(rep(1.25, times = 4), "cm"),
                 txt_gp = fpTxtGp(ticks = gpar(cex = 0.8), xlab = gpar(cex = 1), cex = 0.8),
                 graph.pos = graph.pos,
                 new_page = F
    )
    return(p)
  }else{
    return(mg_getplot_bank('data must be greater than 3 column'))
  }
}

load('xxx.RData')
load('Clinical.RData')
comsamples=intersect(colnames(dt),c(paste0(rownames(clinical),'-01'),paste0(rownames(clinical),'-03'),paste0(rownames(clinical),'-06')))

genes=c('PUS1','PUS3','PUS7','PUS8')
data=log2(dt[,comsamples]+1)
cli.cut=clinical[substr(comsamples,1,12),]
data.final=data[which(apply(data,1,function(x){return(sum(x>0))})>0.1*ncol(data)),]

times=cli.cut$OS.time/365
status=cli.cut$OS
dat111=data.frame(time=times,status=status,t(data.final[genes,]))

Age=as.numeric(cli.cut$A01_Age)
Stage=gsub('[ABC12]','',cli.cut$A07_Stage,perl = T)
Stage[which(Stage=='0')]=NA
Stage[which(Stage=='IS')]=NA
Stage[which(Stage=='X')]=NA

fors1=coxFun(data.frame(dat111$time,dat111$status,ifelse(as.numeric(dat111[,3])>median(as.numeric(dat111[,3])),1,0)))
fors2=coxFun(data.frame(dat111$time,dat111$status,ifelse(as.numeric(dat111[,4])>median(as.numeric(dat111[,4])),1,0)))
fors3=coxFun(data.frame(dat111$time,dat111$status,ifelse(as.numeric(dat111[,5])>median(as.numeric(dat111[,5])),1,0)))
fors4=coxFun(data.frame(dat111$time,dat111$status,ifelse(as.numeric(dat111[,6])>median(as.numeric(dat111[,6])),1,0)))
fors5=coxFun(data.frame(dat111$time,dat111$status,ifelse(Age <= 60,0,1)))
Stage1=Stage
Stage1[Stage1=='I']=0
Stage1[Stage1=='II']=0
Stage1[Stage1=='III']=1
Stage1[Stage1=='IV']=1
fors6=coxFun(data.frame(dat111$time,dat111$status,Stage1)[which(Stage1==0|Stage1==1),])

dat.duo=data.frame(as.numeric(dat111$time),as.numeric(dat111$status),ifelse(as.numeric(dat111[,3])>median(as.numeric(dat111[,3])),1,0)
                   ,ifelse(as.numeric(dat111[,4])>median(as.numeric(dat111[,4])),1,0),ifelse(as.numeric(dat111[,5])>median(as.numeric(dat111[,5])),1,0)
                   ,ifelse(as.numeric(dat111[,6])>median(as.numeric(dat111[,6])),1,0),as.numeric(ifelse(Age <= 60,0,1)),as.numeric(Stage1))
dat.duo=dat.duo[which(dat.duo[,1]!='NA'&dat.duo[,2]!='NA'&dat.duo[,3]!='NA'&dat.duo[,4]!='NA'&dat.duo[,5]!='NA'&dat.duo[,6]!='NA'&dat.duo[,7]!='NA'&dat.duo[,8]!='NA'),]
colnames(dat.duo)=c('time','status','PUS1','PUS3','PUS7','PUS8','Age','Stage')

fmla1 <- as.formula(paste0("Surv(time, status) ~"
                           ,paste0(c('PUS1','PUS3','PUS7','PUS8','Age','Stage'),collapse = '+')))
cox1 <- coxph(fmla1, data =as.data.frame(dat.duo))

HR=c(round(fors1[2],3),round(fors2[2],3),round(fors4[2],3),round(fors4[2],3),round(fors5[2],3),round(fors6[2],3))
lower=c(round(fors1[3],3),round(fors2[3],3),round(fors3[3],3),round(fors4[3],3),round(fors5[3],3),round(fors6[3],3))
upper=c(round(fors1[4],3),round(fors2[4],3),round(fors3[4],3),round(fors4[4],3),round(fors5[4],3),round(fors6[4],3))
p.value=c(round(fors1[1],3),round(fors2[1],3),round(fors3[1],3),round(fors4[1],3),round(fors5[1],3),round(fors6[1],3))
data.test=data.frame(Names=c('PUS1','PUS3','PUS7','PUS8','Age','Stage'),p.value,HR,lower,upper)

forestplot_v1(data.test,xlog = T,colgap = 6,lineheight = 12,xlab = 'Hazard Ratio',box_col=mypal[2],summary_col='black',graph.pos=4)
data.test1=data.frame(Names=names(cox1$coefficients),'p.value'=round(summary(cox1)[[7]][,5],3)
                      ,round(summary(cox1)[[7]][,2],3),round(summary(cox1)[[8]][,3],3),round(summary(cox1)[[8]][,4],3))
forestplot_v1(data.test1,xlog = T,colgap = 6,lineheight = 12,xlab = 'Hazard Ratio',box_col=mypal[2],summary_col='black',graph.pos=4)

library(rms)
library(Hmisc) 

dat=data.frame(dat111$time,dat111$status,dat111[,genes],Age,Stage)
colnames(dat)=c('time','status','PUS1','PUS3','PUS7','PUS8','Age','Stage')

fmla <- as.formula(paste0("Surv(time, status) ~",paste0(c('PUS1','PUS3','PUS7','PUS8','Age','Stage'),collapse = '+')))
cox3 <- cph(fmla, data = dat,surv = T,x = T,y = T)
pbccox <- coxph(fmla, data = dat)

library(regplot)

regplot(pbccox#对观测2的六个指标在列线图上进行计分展示
        #,observation=pbc[1,] #也可以不展示
        #预测3年和5年的死亡风险，此处单位是day
        ,plots = c("bean", #可选"no plot" "density" "boxes" "ecdf" "bars" "boxplot" "violin" "bean" "spikes"
                   "bars")
        ,failtime = c(1,3,5)
        ,subticks = TRUE
        #,clickable=TRUE
        ,prfail = TRUE #cox回归中需要TRUE
        ,showP = T #是否展示统计学差异
        ,droplines = F#观测2示例计分是否画线
        ,colors = mypal #用前面自己定义的颜色
        #,rank="range" #根据统计学差异的显著性进行变量的排序
        ,interval="confidence"
        ,points=TRUE) #展示观测的可信区间

######
f1<-cph(formula = fmla,data=dat,x=T,y=T,surv = T,na.action=na.delete,time.inc = 1) 
cal1<-calibrate(f1, cmethod="KM", method="boot",u=1,m=200,B=200) 
f2<-cph(formula = fmla,data=dat,x=T,y=T,surv = T,na.action=na.delete,time.inc = 3) 
cal2<-calibrate(f2, cmethod="KM", method="boot",u=3,m=200,B=200) 
f3<-cph(formula = fmla,data=dat,x=T,y=T,surv = T,na.action=na.delete,time.inc = 5) 
cal3<-calibrate(f3, cmethod="KM", method="boot",u=5,m=200,B=200) 

cal1 <- calibrate(cox3, u=1, cmethod='KM', m=200, B=200)
cal3 <- calibrate(cox3, u=3, cmethod='KM', m=200, B=200)
cal5 <- calibrate(cox3, u=5, cmethod='KM', m=200, B=200)


plot(cal1,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#A6CEE3"),
     xlim = c(0,1),ylim= c(0,1),col = c("#A6CEE3"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#A6CEE3"), pch = 16)

abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year",'5-year'), #图例文字
       col =c("#2166AC","#B2182B",'#A6CEE3'), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框


