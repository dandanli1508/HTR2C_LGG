load('Immu_score.RData')
im.score=Immu_score[colnames(data.ana),grep('EPI',colnames(Immu_score))]

cor_point=function(x,y,method='Pearson',top_col='#D55E00',right_col='#009E73'
                   ,ylab='y expression',xlab='x expression',title=NULL
                   ,marginal.type=c("histogram", "boxplot", "density", "violin", "densigram")[1]){
  library(ggstatsplot)
  dat=data.frame(X=x,Y=y)
  tp='nonparametric'
  if(method=='Pearson'|method=='pearson'){
    tp='parametric'
  }
  g1=ggscatterstats(data = dat, 
                    x = X, 
                    y = Y
                    ,type = tp
                    ,xfill = top_col
                    ,yfill = right_col
                    ,xlab = xlab
                    ,ylab=ylab
                    ,marginal.type = marginal.type
                    ,title = title)
  return(g1)  
}
#################################
plot.rs=list()
for (aaa in colnames(im.score)) {
  plot.rs[[aaa]]=cor_point(x=as.numeric(risk),y=log2(as.numeric(im.score[,aaa])+1),top_col='red',right_col='blue'
                           ,xlab=paste0('Riskscore')
                           ,ylab=paste0('Log2 (',aaa,' expression)')
                           ,marginal.type='density',method = 'spearman')
}

plot.rs


GS=as.data.frame(data.map$Groups)####分组文件
rownames(GS)=rownames(im.score)
colnames(GS)='Groups'

plotdata <- t(scale(im.score,center = T))
plotdata[plotdata > 2] <- 2
plotdata[plotdata < -2] <- -2

blank <- "   "
p.value <- pvalues
sigcode <- cut(as.numeric(p.value), c(0, 0.001, 0.01, 0.05, 0.1, 1),labels=c('***', '**', '*', '', ''))
sig.label <- as.character(sigcode)
p.label <- formatC(p.value,format = "e",digits = 2)
add.label <- str_pad(paste0(rownames(plotdata),sig.label), # 固定行名宽度并再右侧补齐" "
                     max(nchar(paste0(rownames(plotdata),sig.label))), 
                     side = "right")
plotdata=plotdata[,rownames(GS)]
#annCol=data.gs
#colnames(annCol) <- 'Group'
# colnames(annCol)[1] <- paste(str_pad(colnames(annCol)[1], # 注释列名补上"P-value"，宽度和刚才一致
#                                      max(nchar(paste0(rownames(plotdata),sig.label))), 
#                                      side = "right"),"P-value",sep = blank)
#names(annColors) <- colnames(annCol)[1]
aa1=pheatmap(mat = plotdata, # 输入数据
                   scale = "none", # 不标准化因为数据已经被标准化
                   annotation_col = GS, # 列注释信息
                   color = colorRampPalette(c(mypal[2], "white", mypal[1]))(100),
                   #annotation_colors = annColors, # 列注释对应的颜色
                   cluster_cols = F, # 列不聚类
                   cluster_rows = T, # 行不聚类
                   show_colnames = F, # 不显示列名
                   show_rownames = T, # 显示基因名
                   #annotation_legend = F, # 不显示图例
                   labels_row = paste(add.label, p.label, sep=blank), # 自定义样本名义blank作间隔
                   fontfamily = "mono")

#####################################################
load('XXX.RData')
cli.eac=read.csv('PMC6066282-TCGA-CDR-clinical.txt',stringsAsFactors = F,row.names = 1,check.names = F,sep = '\t')
comsample=intersect(paste0(rownames(cli.eac),'-01'),colnames(dt))

data.ana=log2(dt[,comsample]+1)
data.cli=cli.eac[substr(comsample,1,12),]

genes=c('CD274','CTLA4','HAVCR2','LAG3','PDCD1','PDCD1LG2','TIGIT','SIGLEC15')
dt.che=data.ana[genes,]

ttt0= dt.che %>%
      rownames_to_column('Samples') %>%
      pivot_longer(cols = 2:(ncol(dt.che)+1),names_to='Celltype',values_to='Values')

datas.final=as.data.frame(rbind(datas,datas.nor))
pvalues <- sapply(unique(datas.final$Type), function(x) {
      if (length(unique(datas.final$Group))==2) {
        res <- wilcox.test(Values ~ Group, data = datas.final[which(datas.final$Type=='CD274'),])
        res$p.value
      }else if (length(unique(datas.final$Group))>=3){
        res <- kruskal.test(Values ~ Group, data = datas.final[which(datas.final$Type==x),])
        res$p.value
      }
    })
pv <- data.frame(Type = unique(datas.final$Type), pvalue = pvalues)
pv$sigcode <- cut(as.numeric(pv$pvalue), c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                 labels=c('***', '**', '*', '', ''))
					  
p1=ggboxplot(datas.final,x='Type',y='Values',color='Group',palette = "nejm",size = 0.5,shape=16,
                   add = "jitter",add.params=list(size=0.1),bxp.errorbar =T,width=0.55,outlier.shape=NA,ggtheme=theme_bw())
p1=p1+geom_text(aes(Type, y=max(datas.final$Values)*1.1,label=pv$sigcode),data=pv, inherit.aes=F,size=6) + xlab(NULL)+ylab('Immune checkpoint')
p1=p1+theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="serif",size = 10)
                ,axis.text.y=element_text(family="serif",face="plain",size = 10)
                ,axis.title.y=element_text(family="serif",face="plain",size = 10)
                ,panel.border = element_blank(),axis.line = element_line(colour = "black")
                ,legend.text=element_text(face="plain", family="serif", colour="black",size = 10)
                ,legend.title=element_text(face="plain", family="serif", colour="black",size = 10)
                #,legend.justification=c(1,1), legend.position=c(1,1)
                ,legend.background = element_rect(fill = NA, colour = NA)
                ,panel.grid.major = element_blank()
                ,panel.grid.minor = element_blank())
if (length(unique(datas.final$Group))==2) {
    fin.box=p1+labs(fill =paste0("  wilcox.test","\n\n"," * p < 0.05","\n\n","** p < 0.01","\n\n","Groups"))+theme(text=element_text(size=8,family="serif"))
}else if (length(unique(datas.final$Group))>=3){
    fin.box=p1+labs(fill =paste0("  kruskal.test","\n\n"," * p < 0.05","\n\n","** p < 0.01","\n\n","Groups"))+theme(text=element_text(size=8,family="serif"))
}


