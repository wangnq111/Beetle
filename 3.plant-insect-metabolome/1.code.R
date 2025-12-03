library(ropls)
library(tidyverse)
library(LorMe)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(vegan)
library(tidyverse)
library(magrittr)

metabolite=read.csv("./input/metabolism.csv",header = T,encoding = "UTF-8",as.is = TRUE)
metafile=read.table("./input/metafile.txt",header = T,sep="\t")
plant_col=color_scheme("Plan7",names=unique(metafile$Plant))
plant_order=c("Locust","Rose","Apricot","Seabuckthorn","Willow","Apple","Elm","Filbert")
sex_col=c("F"="#FE5D5D",FM="#6376A0")

kegg_anno=data.frame(KeggID=paste0("kegg",1:nrow(metabolite)),metabolite[,1:34])

metabolite_pct=sweep(metabolite[,35:82],colSums(metabolite[,35:82]),"/",MARGIN=2)




# OPLSDA
plsda_input=metabolite[,35:82] %>% t()
colnames(plsda_input)=metabolite$Compound_ID

plsda<-opls(plsda_input, metafile$Plant,
            predI = 2, orthoI = 0,log10L = T,crossvalI =3,scaleC="pareto",permI=999)#,#orthoI = 1:opls-da; orthoI = 0:Pls-da
#log10L = F,
#crossvalI =10,
#scaleC="pareto",#pareto scaling
#fig.pdfC = c("none", "interactive", "outputfile/plsda.pdf")[3],
#permI=200)

# data and score data
## score data
plsdaScore <- data.frame(
  t1 =plsda@scoreMN[,1],
  to1 =plsda@scoreMN[,2] ) %>% 
  scale(center = T,scale = T) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'Sample') %>% 
  merge(metafile,by="Sample")

t1Weight=sprintf("%.1f%%", plsda@modelDF[1,1] * 100);t1Weight
to1Weight=sprintf("%.1f%%", plsda@modelDF[2,1] * 100);to1Weight

R2X=plsda@modelDF[1,1]+plsda@modelDF[2,1]
R2Y=plsda@modelDF[1,3]+plsda@modelDF[2,3]
Q2Y=plsda@modelDF[1,6]+plsda@modelDF[2,6]

subTitle <- paste0("R2X=",R2X,"  R2Y=",R2Y,"  Q2Y=",Q2Y)
cent = aggregate(plsdaScore[, 2:3], by = list(Plant=plsdaScore$Plant),FUN = mean)
colnames(cent)[2:3] = c("cent1", "cent2")
plsdaScore=merge(plsdaScore,cent,by=c("Plant"))
## opls-da 
oplsda=ggplot(plsdaScore,aes(x=t1,y=to1,fill=Plant))+
  
  geom_segment(aes(xend = cent1, yend = cent2,color=Plant),show.legend = FALSE, size = 0.3, alpha = 0.8)+
  geom_point(size=1.5,alpha=.8,aes(shape=factor(sex)))+
  geom_hline(yintercept = 0, linetype = 4,color = "grey") +
  geom_vline(xintercept = 0, linetype = 4,color = "grey")+
  theme_zg()+
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values=plant_col)+
  scale_color_manual(values=plant_col)+
  labs(#title = "OPLS-DA", 
    #subtitle = subTitle, 
    x = paste0("T score1 (",t1Weight,")"),
    y = paste0("T score2 (",to1Weight,")"))+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)), 
    plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)) 
  )
#ggsave("./results/structure/plsda.pdf",width = 5,height = 4)


##stat
set.seed(999)
adonis2(plsda_input~metafile$Plant*metafile$sex, by = "terms",method="euclidean")


##composition####
metabolite_summary=aggregate(metabolite_pct,by=list(metabolite$ClassI),FUN=sum) 
metabolite_summary$Group.1[rowMeans(metabolite_summary[,-1])<=(rowMeans(metabolite_summary[,-1]) %>% sort(decreasing = T) %>% .[10])]="Others"
metabolite_summary=aggregate(metabolite_summary[,-1],by=list(metabolite_summary$Group.1),FUN=sum)
metabolite_summary_combine=data.frame(metafile,t(metabolite_summary[,-1]))
colnames(metabolite_summary_combine)[-c(1:6)]=metabolite_summary$Group.1
metabolite_summary_combine_mean=aggregate(metabolite_summary_combine[,-c(1:6)],by=list(metabolite_summary_combine$Plant),FUN=mean)  #,metabolite_summary_combine$sex
metabolite_summary_long=gather(metabolite_summary_combine_mean,"classify","Rel",-c(Group.1))

metabolite_order=metabolite_summary$Group.1[rowMeans(metabolite_summary[,-1]) %>% order(.,decreasing=T)] %>%.[c(1:8,10,9)]
metabolite_summary_long$Group.1=factor(metabolite_summary_long$Group.1,levels =plant_order)
#metabolite_summary_long$Group.2=factor(metabolite_summary_long$Group.2,levels =c('Susceptible','Tolerant'))
metabolite_summary_long$classify=factor(metabolite_summary_long$classify,levels=metabolite_order)
library(ggalluvial)
classfy_col=color_scheme("Plan8",names = metabolite_order)


comp=ggplot(metabolite_summary_long, aes(x = Group.1, y = Rel, fill = classify)) + 
  #facet_grid(~Group.1,scale="free")+
  geom_col(position = "stack", width = 0.8)+
  #geom_flow(alpha = 0.5, width = 0.1)+
  scale_y_continuous(expand = c(0,0),labels = scales::percent)+
  scale_x_discrete(expand = c(0.1,0))+
  theme_zg()+
  labs(x="",y="Relative abundance")+
  scale_fill_manual(values=classfy_col)+
  theme(legend.text = element_text(size = 8, face = "bold"), 
        strip.text = element_text(size = 8, face = "bold"), 
        axis.ticks.x = element_blank(), panel.background = element_rect(color = NA), 
        axis.line.y = element_line(size = 0.4, color = "black"), 
        axis.ticks.length.y = unit(0.4, "lines"), axis.ticks.y = element_line(color = "black", 
                                                                              size = 0.4))
comp
ggsave("./results/community/composition.pdf",width = 6,height = 3)

#compunds stat####
phylum_stat=data.frame()
for(i in colnames(metabolite_summary_combine)[-c(1:6,15)]){
  y=metabolite_summary_combine[,i]
  stat=aov(y~Plant*sex,data=metabolite_summary_combine) %>% summary()%>% .[[1]] %>% as.data.frame()
  stat$Item=rownames(stat)
  stat$signif=case_when(
    stat$`Pr(>F)`>0.05 ~"ns",
    stat$`Pr(>F)`>0.01&stat$`Pr(>F)`<0.05 ~"*",
    stat$`Pr(>F)`>0.001&stat$`Pr(>F)`<0.01 ~"**",
    stat$`Pr(>F)`<0.001 ~"***",
    is.na( stat$`Pr(>F)`)==T~NA
  )
  stat$phylum=i
  phylum_stat=rbind(phylum_stat,stat)
}




phylum_stat$sign=case_when(phylum_stat$`Pr(>F)`>0.05~'ns',
                           phylum_stat$`Pr(>F)`<0.05&phylum_stat$`Pr(>F)`>0.01~'*',
                           phylum_stat$`Pr(>F)`<0.01&phylum_stat$`Pr(>F)`>0.001~'**',
                           phylum_stat$`Pr(>F)`<0.001~'***')

write.table(phylum_stat,"./results/community/classI_anova.txt",row.names = F,sep="\t")

#pro test####
library(vegan)
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1
host_matrix <- model.matrix(~ Plant + sex - 1, data = metafile)
host_hell <- decostand(host_matrix, "hellinger")

# 2
host_ord <- rda(host_hell)
host_coords <- scores(host_ord, choices = 1:2)

# 3
metab_hell <- decostand(t(metabolite_pct), "hellinger")
metab_ord <- rda(metab_hell)
metab_coords <- scores(metab_ord, choices = 1:2)

# 4
proc <- protest(host_ord, metab_ord, permutations = 999)
summary(proc)

# 5
df_host <- as.data.frame(host_coords$sites) %>%
  mutate(ID = rownames(.),
         Group = "Plant+sex",
         HostType = metafile$Plant,
         Sex = metafile$sex)%>% data.frame(metafile)

df_metab <- as.data.frame(metab_coords$sites) %>%
  mutate(ID = rownames(.),
         Group = "Metabolites") %>% as.data.frame()%>% left_join(data.frame(ID=metafile$Sample,metafile))

df_all <- bind_rows(df_host, df_metab)

# 
arrows <- data.frame(
  x = df_host$PC1, y = df_host$PC2,
  xend = df_metab$PC1, yend = df_metab$PC2,
  ID = df_host$ID
) %>% data.frame(metafile)

# 7
pro_plot=ggplot() +
  geom_point(data = df_metab, aes(x = PC1, y = PC2,fill=Plant),
             alpha = 0.8, size = 1.5,pch=24) +
  geom_point(data = df_host, aes(x = PC1, y = PC2,
                                 fill = HostType, shape = Sex),
             size = 1.5,alpha = 0.8) +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values=plant_col)+
  scale_color_manual(values=plant_col)+
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend,color=Plant),
               arrow = arrow(length = unit(0.01, "cm")),
               alpha = 0.5,linewidth=.5) +
  theme_zg() +
  labs(x = "Procrustes PC1", y = "Procrustes PC2",
       color = "HostType", shape = "Sex")
#core kegg####
kegg_mean_plant=aggregate(t(metabolite_pct),by=list(Plant=metafile$Plant),FUN=mean)
###exist species####
exist=function(x){
  which(x>0) %>% length() %>% return()
}
exist_kegg=apply(kegg_mean_plant[,-1], 1, exist)


###high abundant kegg####
threshold=1e-4
high_abundant=function(x){
  which(x>threshold) %>% length() %>% return()
}
high_abundant_kegg=apply(kegg_mean_plant[,-1], 1, high_abundant)

###core high abundant kegg####
core_high_kegg=intersect(
  colnames(kegg_mean_plant)[which(kegg_mean_plant[1,-1]>threshold)+1],
  colnames(kegg_mean_plant)[which(kegg_mean_plant[2,-1]>threshold)+1]
)
i=3
while(i<=8){
  core_high_kegg=intersect(
    core_high_kegg,
    colnames(kegg_mean_plant)[which(kegg_mean_plant[i,-1]>threshold)+1]
  )
  i=i+1
}

###flower visualizaiton####
core_stat=data.frame(Plant=kegg_mean_plant$Plant,
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
flower=ggplot()+
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
#ggsave("./results/kegg_function/core_flower_kegg.pdf",width = 5,height = 4)






###indicator####
indic_data=as.data.frame(t(metabolite_pct))
colnames(indic_data)=paste0("kegg",1:ncol(indic_data))

kegg_indicator=indicspecies::multipatt(indic_data, 
                                       metafile$Plant, func = "r.g", control = permute::how(nperm = 1000)) %$% 
  as.data.frame(sign)
kegg_indicator=data.frame(KeggID=rownames(kegg_indicator),kegg_indicator)
kegg_indicator_anno = left_join(kegg_indicator, kegg_anno)

pie_list=list()
for (i in unique(metafile$Plant)) {
  kegg_indicator_anno_select=
    kegg_indicator_anno[kegg_indicator_anno[,paste0("s.", i) ]== 1 & kegg_indicator_anno$p.value <  0.05&kegg_indicator_anno$stat>0.5,] #
  type_frame=table(kegg_indicator_anno_select$ClassI) %>% sort(decreasing = TRUE) %>% as.data.frame()
  type_frame$Var1=as.character(type_frame$Var1)
  type_frame$Var1[!type_frame$Var1 %in%levels(metabolite_summary_long$classify) ]="Others"
  type_frame=aggregate(type_frame$Freq,by=list(Var1=type_frame$Var1),FUN=sum)
  type_frame$Var1=factor(type_frame$Var1,levels = levels(metabolite_summary_long$classify))
  colnames(type_frame)[2]="Freq"
  
  p=ggplot(type_frame, aes_string(x = "0",y = "Freq", fill = "Var1")) + 
    geom_bar(width = 1,  stat = "identity", show.legend = TRUE) + 
    coord_polar("y", start = 0) + 
    scale_fill_manual(values = classfy_col) + 
    labs(x = "", y = "", title = i, fill = "Type") + 
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), 
          axis.text = element_blank(),  
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 8))
  pie_list=c(pie_list,list(p))
  names(pie_list)[length(pie_list)]=i
}

up_pie=(pie_list$Apple|pie_list$Locust|pie_list$Willow|pie_list$Rose)*guides(fill="none")
down_pie=(pie_list$Seabuckthorn|pie_list$Filbert|pie_list$Apricot|pie_list$Elm)*guides(fill="none")


up_pie/down_pie
#ggsave("./results/kegg_function/indicator_pie.pdf",width = 8,height = 4)

#top kegg####
##Top 20 metabolite####
metabolite_pct_order=metabolite_pct[order(rowMeans(metabolite_pct),decreasing = T),] %>% .[1:20,]
metabolite_anno=metabolite[,1:18]
metabolite_anno_order=metabolite_anno[order(rowMeans(metabolite_pct),decreasing = T),] %>% .[1:20,]
metabolite_anno_order$mean=rowMeans(metabolite_pct)[order(rowMeans(metabolite_pct),decreasing = T)]%>% .[1:20]
rownames(metabolite_pct_order)=metabolite_anno_order$Name

metabolite_cate3_mean=aggregate(t(metabolite_pct_order),by=list(Plant=metafile$Plant),FUN=mean)

metabolite_cate3_mean_heat=t(metabolite_cate3_mean[,-1]) %>% as.data.frame()
colnames(metabolite_cate3_mean_heat)=metabolite_cate3_mean$Plant
rownames(metabolite_cate3_mean_heat)=metabolite_anno_order$Name

metabolite_anno_order_heat=data.frame(Category1=metabolite_anno_order$ClassI,row.names =metabolite_anno_order$Name )
ann_colors=list(Category1=classfy_col)
library(pheatmap)
#heatamp
pheatmap(metabolite_cate3_mean_heat,scale="row",cluster_rows = F,
         annotation_row =metabolite_anno_order_heat,annotation_colors = ann_colors,border_color = NA,
         cellwidth  = 15,cellheight = 15)

#signif stat
cate3_stat=data.frame()
for(i in rownames(metabolite_pct_order)){
  subdata=metabolite_pct_order[i,] %>% as.numeric()
  stat=aov(subdata~Plant*sex,data=metafile) %>% summary() %>% .[[1]] %>% as.data.frame()
  stat$Item=rownames(stat)
  stat$signif=case_when(
    stat$`Pr(>F)`>0.05 ~"ns",
    stat$`Pr(>F)`>0.01&stat$`Pr(>F)`<0.05 ~"*",
    stat$`Pr(>F)`>0.001&stat$`Pr(>F)`<0.01 ~"**",
    stat$`Pr(>F)`<0.001 ~"***",
    is.na( stat$`Pr(>F)`)==T~NA
  )
  stat$Name=i
  cate3_stat=rbind(cate3_stat,stat)
}


cate3_stat=left_join(cate3_stat,metabolite_anno)
#write.table(cate3_stat,"./results/KEGG_function/kegg_anova.txt",row.names = F,sep="\t")
##
top20_stat=cate3_stat
top20_stat=top20_stat[top20_stat$Item !="Residuals  ",]
top20_stat$Name=factor(top20_stat$Name,levels =rownames(metabolite_anno_order_heat) %>% rev() )
top20_stat$Item=factor(top20_stat$Item,levels = unique(top20_stat$Item))
signif_col= color_scheme("Plan3") %>% .[1:3] %>% color_scheme(expand = 4)
names(signif_col)=c("***","**","*","ns")

p1=ggplot(top20_stat,
          aes(x=Item,y=Name))+
  geom_tile(color=NA,size=0.3,aes(fill = signif),alpha=.6)+
  geom_text(aes(label=signif))+
  scale_fill_manual(values =signif_col)+
  theme_zg()+
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank())
#relative abundance
metabolite_anno_order$Name=factor(metabolite_anno_order$Name,levels =rownames(metabolite_anno_order_heat) %>% rev() )
p2=ggplot(metabolite_anno_order,aes(x=mean,y=Name))+
  geom_vline(xintercept = c(0.05,0.1,0.15,0.2),size=.2,lty=2)+
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
#ggsave("./top20metabolite_stat.pdf",width = 5,height = 4.7)

(oplsda|pro_plot|comp|flower)*guides(fill="none",color="none",shape="none")+plot_layout(guides = "collect")
ggsave("./figure6_top.pdf",width = 9,height = 3)


save(classfy_col,metabolite_pct,metabolite,file="metabolite.rda")
