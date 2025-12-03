library(factoextra)
library(FactoMineR)
library(magrittr)
library(LorMe)
library(dplyr)
library(tidyr)
library(vegan)

kegg=read.csv("./data/keggL3.csv")
kegg_cato=read.csv("./data/KEGG_catogory.csv") 
kegg_cato=unique(kegg_cato[,c(3:5)])
kegg_pct=sweep(kegg[,-1],colSums(kegg[,-1]),"/",MARGIN=2)
metafile=read.table("./data/metafile.txt",header = T,sep="\t")
plant_col=color_scheme("Plan7",names=unique(metafile$Plant))

res.pca <- PCA(t(kegg_pct), 
               scale.unit = F,# 默认对数据进行scale标准化，最后响应变量间的均值为0，sd=1。
               ##相同尺度的标准化避免了一些变量仅仅因为其较大的度量单位而成为主导变量。使变量具有可比性。
               ncp = 11, # 结果中默认保留5，nrow(env[-ind.sup,])-1(样本数)和ncol(env[2:14][-c(quanti.sup,quali.sup)])(响应变量数)中的最小值个主成分
               #ind.sup = 1, # 附加样本，1表示数据表的第一行，不是必须的数据，默认为NULL。附加样本的坐标将根据PCA分析结果进行预测。
               quanti.sup = NULL, # 附加的定量响应变量数据，不是必须的，默认为NULL。基于PCA分析结果预测坐标。
               quali.sup = NULL, # 附加定性变量，分类信息最好不要用数字表示，在行数区间的数字，运行会报错。
               ##1:2表示定性变量在env[2:14]中的列数。不是必须的，默认为NULL。基于PCA分析结果预测坐标，可用于对样本进行分组着色。
               row.w = NULL, col.w = NULL, # 分别给样本和响应变量设定权重，默认所有样本或所有响应变量的权重都是一样的。默认相当于对变量数据进行了中心化，变量均值接近0。可以设置为1个数字，也可以设置为分别与行数或列数相等的矢量。
               graph = FALSE)

sam=data.frame(res.pca$ind$coord) %>% .[1:2] %>% data.frame(metafile)
cent = aggregate(sam[, 1:2], by = list(sam$Plant),FUN = mean)
colnames(cent) = c("Plant", "cent1", "cent2")
PCAframe = left_join(data.frame(sam, joint = sam$Plant), 
                     data.frame(cent[, 2:3], joint = cent$Plant)) %>% suppressMessages()


ggplot(PCAframe,aes(x=Dim.1,y=Dim.2,fill=factor(Plant)))+
  geom_hline(yintercept=0,linetype=4,color="grey") +
  geom_vline(xintercept=0,linetype=4,color="grey") +
  geom_segment(aes(xend = cent1, yend = cent2,color=Plant),show.legend = FALSE, size = 0.3, alpha = 0.6)+
  geom_point(size=1.5,alpha=.8,aes(shape=sex))+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = plant_col)+
  scale_fill_manual(values = plant_col)+
  labs(x=paste0("PC1:",get_eigenvalue(res.pca) %>% .[1,2]%>% round(.,2),"%"),
       y=paste0("PC1:",get_eigenvalue(res.pca) %>% .[2,2]%>% round(.,2),"%"),
       fill="")+
  # stat_ellipse(level=0.95)+
  theme_zg()+guides(color="none")#+
  #theme(legend.position = c(.85,.7),legend.background = element_rect(color="black",linetype = 2,linewidth = .5))
adonis2(t(kegg_pct)~metafile$Plant*metafile$sex,method="euclidean",by = "terms")
ggsave("./results/KEGG_function/function_pcoa.pdf",width = 4.5,height = 2.8)

indic_data=as.data.frame(t(kegg_pct))
colnames(indic_data)=paste0("Kegg",1:ncol(indic_data))
                         
kegg_indicator=indicspecies::multipatt(indic_data, 
                        metafile$Plant, func = "r.g", control = permute::how(nperm = 1000)) %$% 
  as.data.frame(sign)
kegg_indicator=data.frame(keggID=rownames(kegg_indicator),kegg_indicator)
taxonomy=data.frame(keggID=paste0("Kegg",1:ncol(indic_data)),Category3=kegg$Category) %>%
  left_join(kegg_cato)
kegg_indicator_anno = left_join(kegg_indicator, taxonomy)

kegg_indicator_anno$tag="None"
for (i in unique(metafile$Plant)) {
  kegg_indicator_anno$tag[kegg_indicator_anno[,paste0("s.", i)] == 1 & 
                               rowSums(kegg_indicator_anno[, 1:length(unique(metafile$Plant))] == 1) & 
                            kegg_indicator_anno$p.value <  0.05] = i
}


cate2_mean=aggregate(kegg_pct,by=list(Category2=taxonomy$Category2),FUN=sum)

cate2_stat=data.frame()
for(i in (unique(cate2_mean$Category2)%>% as.character())){
  subdata=cate2_mean[cate2_mean$Category2==i,-1] %>% as.numeric()
  stat=aov(subdata~Plant*sex,data=metafile) %>% summary() %>% .[[1]] %>% as.data.frame()
  stat$Item=rownames(stat)
  stat$signif=case_when(
    stat$`Pr(>F)`>0.05 ~"ns",
    stat$`Pr(>F)`>0.01&stat$`Pr(>F)`<0.05 ~"*",
    stat$`Pr(>F)`>0.001&stat$`Pr(>F)`<0.01 ~"**",
    stat$`Pr(>F)`<0.001 ~"***",
    is.na( stat$`Pr(>F)`)==T~NA
  )
  stat$Category2=i
  cate2_stat=rbind(cate2_stat,stat)
}
kegg_anno=data.frame(Category3=kegg$Category) %>% left_join(kegg_cato)
cate2_mean_anno=left_join(cate2_mean,unique(kegg_anno[,2:3]))
cate2_mean_anno=cate2_mean_anno[cate2_mean_anno$Category1=="Metabolism",]
cate2_mean_anno=cate2_mean_anno[cate2_mean_anno$Category2 %in% (cate2_stat$Category2[cate2_stat$Item=="Plant      "&cate2_stat$`Pr(>F)`<0.05]),]
cate2_mean_anno_long=gather(cate2_mean_anno[,-50],"Sample","Value",-c(Category2)) %>% left_join(metafile[,1:3])

inputdata<- cate2_mean_anno_long 
facet_compare<-data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]
for(i in unique(inputdata[,1])){
  sub_facet<-inputdata[which((inputdata[,1])==i),] 
  results<-auto_signif_test(data =sub_facet,treatment_col =4,value_col =3,prior = T)
  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=max(1.3*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std'])),facetlabel=i))
}

colnames(facet_compare)[7]=colnames(inputdata)[1]
mean_frame<- aggregate(inputdata[,3],by=list(inputdata[,4],inputdata[,1]),FUN=mean)
Sd<- aggregate(inputdata[,3],by=list(inputdata[,4],inputdata[,1]),FUN=sd)%>% .[,'x']
Treatment_Name<- mean_frame$Group.1
N<- table(inputdata[,4]) %>% as.numeric()
Mean<- mean_frame[,'x']
SEM<- Sd/(N^0.5)
input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])
colnames(input_mean_frame)[6]=colnames(inputdata)[1]
ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
  theme_zg()+
  geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name)),color='#000000',alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
  scale_y_continuous(expand = c(0.01,0.01))+
  scale_color_manual(values=plant_col)+
  scale_fill_manual(values=plant_col)+
  labs(x='',fill='',y=colnames(inputdata)[3])+
  facet_wrap(~Category2,scales ='free_y',strip.position ='top',as.table =T)+
  geom_text(data=facet_compare,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
  theme(legend.position = 'right',
  )

ggsave("./results/KEGG_function/key_function_bar.pdf",width = 15,height = 4)

cate3=data.frame(Category=kegg$Category,kegg_pct)
cate3_stat=data.frame()
for(i in (unique(cate3$Category)%>% as.character())){
  subdata=cate3[cate3$Category==i,-1] %>% as.numeric()
  stat=aov(subdata~Plant*sex,data=metafile) %>% summary() %>% .[[1]] %>% as.data.frame()
  stat$Item=rownames(stat)
  stat$signif=case_when(
    stat$`Pr(>F)`>0.05 ~"ns",
    stat$`Pr(>F)`>0.01&stat$`Pr(>F)`<0.05 ~"*",
    stat$`Pr(>F)`>0.001&stat$`Pr(>F)`<0.01 ~"**",
    stat$`Pr(>F)`<0.001 ~"***",
    is.na( stat$`Pr(>F)`)==T~NA
  )
  stat$Category3=i
  cate3_stat=rbind(cate3_stat,stat)
}


cate3_stat=left_join(cate3_stat,taxonomy)
write.table(cate3_stat,"./results/KEGG_function/kegg_anova.txt",row.names = F,sep="\t")


cate3_select=cate3_stat[cate3_stat$`Pr(>F)`<0.05&cate3_stat$Item=="Plant      "&cate3_stat$Category1=="Metabolism",]
cate3_select=cate3_select[order(cate3_select$`Pr(>F)`,decreasing = F),]
select_cate=cate3_select$Category3[1:10]

rownames(kegg_pct)=taxonomy$keggID
select_cate_long=data.frame(metafile[,1:2],t(kegg_pct[kegg$Category %in% select_cate,])) %>% gather("keggID","Value",-c(Sample,Plant))
select_cate_long=left_join(select_cate_long,taxonomy)


inputdata<- select_cate_long 
facet_compare<-data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]
for(i in unique(inputdata[,5])){
  sub_facet<-inputdata[which((inputdata[,5])==i),] 
  results<-auto_signif_test(data =sub_facet,treatment_col =2,value_col =4,prior = T)
  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=(1.3*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std']))%>%max(),facetlabel=i))
}

colnames(facet_compare)[7]=colnames(inputdata)[5]
mean_frame<- aggregate(inputdata[,4],by=list(inputdata[,2],inputdata[,5]),FUN=mean)
Sd<- aggregate(inputdata[,4],by=list(inputdata[,2],inputdata[,5]),FUN=sd)%>% .[,'x']
Treatment_Name<- mean_frame$Group.1
N<- table(inputdata[,2]) %>% as.numeric()
Mean<- mean_frame[,'x']
SEM<- Sd/(N^0.5)
input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])
colnames(input_mean_frame)[6]=colnames(inputdata)[5]
ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
  theme_zg()+
  geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name)),color='#000000',alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values=plant_col)+
  scale_fill_manual(values=plant_col)+
  scale_y_continuous(labels = scales::label_percent())+
  labs(x='',fill='',y="Relative abundance")+
  facet_wrap(~get(colnames(inputdata)[5]),scales ='free_y',strip.position ='top',as.table =T,nrow=2)+
  geom_text(data=facet_compare,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
  theme(legend.position = 'right',
  )
ggsave("./results/KEGG_function/key_function_bar.pdf",width = 15,height = 4)


##core kegg####
kegg_trans=t(kegg_pct)
kegg_trans=data.frame(metafile,kegg_trans)
speices_mean_plant=aggregate(kegg_trans[,-c(1:6)],by=list(Plant=kegg_trans$Plant),FUN=mean)

###exist function####
exist=function(x){
  which(x>0) %>% length() %>% return()
}
exist_kegg=apply(speices_mean_plant[,-1], 1, exist)
###high abundant kegg####
threshold=1e-4
high_abundant=function(x){
  which(x>threshold) %>% length() %>% return()
}
high_abundant_kegg=apply(speices_mean_plant[,-1], 1, high_abundant)

###core high abundant kegg####
core_high_kegg=intersect(
  colnames(speices_mean_plant)[which(speices_mean_plant[1,-1]>threshold)],
  colnames(speices_mean_plant)[which(speices_mean_plant[2,-1]>threshold)]
)
i=3
while(i<=8){
  core_high_kegg=intersect(
    core_high_kegg,
    colnames(speices_mean_plant)[which(speices_mean_plant[i,-1]>threshold)]
  )
  i=i+1
}
###flower visualizaiton####
core_stat=data.frame(Plant=speices_mean_plant$Plant,
                     All=exist_kegg,
                     High=high_abundant_kegg,
                     core=length(core_high_kegg))
core_stat$rare=core_stat$All-core_stat$High
core_stat$high_left=core_stat$High-core_stat$core

library(tidyverse)
x<-1:(8*30)
y<-sin(8*x*pi/240)
df2<-data.frame(x1=x,
                y1=abs(y),
                var=gl(8,30,labels = as.character(core_stat$Plant)%>% sort)) %>% merge(.,core_stat,by.x = 'var',by.y = 'Plant')
df2$varname=paste0(df2$var,df2$all)
ggplot()+
  geom_area(data = data.frame(x=1:240,y=0.5),
            aes(x=x,y=y),
            fill="gray90",
            alpha=1)+
  ggnewscale::new_scale_fill()+
  geom_area(data=df2,aes(x=x1,y=y1,fill=var),
            #fill="blue",
            alpha=0.8)+
  scale_fill_manual(values = plant_col)+
  coord_polar()+
  geom_area(data = data.frame(x=1:240,y=0.2),
            aes(x=x,y=y),
            fill="gray90",
            alpha=1)+
  theme_void()+
  geom_text(data=df2 %>% filter(y1==1),
            aes(x=x1,y=y1+0.2,label=varname),
            fontface="italic",size=4)+
  geom_text(data=df2 %>% filter(y1==1),
            aes(x=x1,y=y1-0.2,label=rare),size=3)+
  geom_text(data=df2 %>% filter(y1==1),
            aes(x=x1,y=y1-0.6,label=high_left),size=3)+
  annotate(geom = "text",x=1,y=0,label=paste0("Core ",unique(df2$core)),size=3)+
  #labs(title="High ")+
  theme(legend.position = "none",plot.title=element_text(hjust=0.5))
ggsave("./results/KEGG_function/kegg_core_flower.pdf",width = 5,height = 4)



##Top 20 Kegg####
kegg_pct_order=kegg_pct[order(rowMeans(kegg_pct),decreasing = T),] %>% .[1:20,]
kegg_anno=data.frame(Category3=kegg$Category) %>% left_join(kegg_cato)
kegg_anno_order=kegg_anno[order(rowMeans(kegg_pct),decreasing = T),] %>% .[1:20,]
kegg_anno_order$mean=rowMeans(kegg_pct)[order(rowMeans(kegg_pct),decreasing = T)]%>% .[1:20]
  
kegg_cate3_mean=aggregate(t(kegg_pct_order),by=list(Plant=metafile$Plant),FUN=mean)

kegg_cate3_mean_heat=t(kegg_cate3_mean[,-1]) %>% as.data.frame()
colnames(kegg_cate3_mean_heat)=kegg_cate3_mean$Plant
rownames(kegg_cate3_mean_heat)=kegg_anno_order$Category3

kegg_anno_order_heat=data.frame(Category1=kegg_anno_order$Category1,row.names =kegg_anno_order$Category3 )
ann_colors=list(Category1=color_scheme("Plan2",4,names = unique(kegg_anno_order_heat$Category1)))

#heatamp
library(pheatmap)
pheatmap(kegg_cate3_mean_heat,scale="row",cluster_rows = F,
         annotation_row =kegg_anno_order_heat,annotation_colors = ann_colors,border_color = NA,
         cellwidth  = 15,cellheight = 15)

#signif stat
top20_stat=cate3_stat[cate3_stat$Category3 %in% kegg_anno_order$Category3,]
top20_stat=top20_stat[top20_stat$Item !="Residuals  ",]
top20_stat$Category3=factor(top20_stat$Category3,levels =rownames(kegg_anno_order_heat) %>% rev() )
top20_stat$Item=factor(top20_stat$Item,levels = unique(top20_stat$Item))
signif_col= color_scheme("Plan3") %>% .[1:3] %>% color_scheme(expand = 4)
names(signif_col)=c("***","**","*","ns")

p1=ggplot(top20_stat,
       aes(x=Item,y=Category3))+
  geom_tile(color=NA,size=0.3,aes(fill = signif),alpha=.6)+
  geom_text(aes(label=signif))+
  scale_fill_manual(values =signif_col)+
  theme_zg()+
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank())
#relative abundance
kegg_anno_order$Category3=factor(kegg_anno_order$Category3,levels =rownames(kegg_anno_order_heat) %>% rev() )
p2=ggplot(kegg_anno_order,aes(x=mean,y=Category3))+
  geom_vline(xintercept = c(0.05,0.1,0.15),size=.2,lty=2)+
  geom_col()+
  scale_x_continuous(expand = c(0,0),labels = scales::label_percent())+
  theme_zg()+
  theme(axis.ticks.y =element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = "black",size=0.2),
        axis.title=element_blank())
library(patchwork)
(p1|p2)+plot_layout(guides = "collect",widths = c(1,3))
ggsave("./top20kegg_stat.pdf",width = 4,height = 5)

write.table(PCAframe,"Fig4a.txt",row.names = F,sep="\t")
write.table(kegg_cate3_mean_heat,"Fig4c.txt",row.names = T,sep="\t")
write.table(top20_stat,"Fig4cmiddle.txt",row.names = F,sep="\t")
write.table(kegg_anno_order,"Fig4cright.txt",row.names = F,sep="\t")
