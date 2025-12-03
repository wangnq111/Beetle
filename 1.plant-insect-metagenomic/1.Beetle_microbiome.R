library(LorMe)
library(magrittr)
library(ggplot2)
library(dplyr)
library(picante)
library(ggbeeswarm)

Bac=read.csv("./data/bacnr.csv",header = T)
Bac=data.frame(TAXID=paste0("Taxon",1:nrow(Bac)),Bac[,-1],taxonomy=Bac$Tax)
metafile=read.table("./data/metafile.txt",header = T,sep="\t")

plant_col=color_scheme("Plan7",names=unique(metafile$Plant))
plant_order=c("Locust","Rose","Apricot","Seabuckthorn","Willow","Apple","Elm","Filbert")
sex_col=c("F"="#FE5D5D",FM="#6376A0")

#filtering####
Bac_clean=Filter_function(input = Bac,threshold = 5e-6,format = 4)
colSums(Bac_clean[,2:49]) %>% sort()

#packaging####
tax_obj=tax_summary(groupfile = metafile,inputtable = Bac_clean[,2:49],reads = T,taxonomytable = Bac_clean[,c(1,50)],outputtax = "standard")
tax_obj_plan1=object_config(taxobj =tax_obj,treat_location = 2,facet_location = 3,rep_location = 4 )
tax_obj_plan2=object_config(taxobj =tax_obj,treat_location = 3,facet_location =2,rep_location = 4 )
tax_obj_plan_plant=object_config(taxobj =tax_obj,treat_location = 2,rep_location = 4,
                                 treat_col = plant_col,treat_order = plant_order)
save(tax_obj_plan1,file="./tax_obj_plan1.rda")

#alpha diversity####
alpha_results=Alpha_diversity_calculator(taxobj =tax_obj_plan_plant,taxlevel = "Genus")
alphaframe=alpha_results$alphaframe
alphaframe[alphaframe$Indexname=="Shannon",] %>%
  aov(Indexvalue~Plant*sex,data=.) %>% summary()
alphaframe[alphaframe$Indexname=="Chao",] %>%
  aov(Indexvalue~Plant*sex,data=.) %>% summary()
alphaframe[alphaframe$Indexname=="Evenness",] %>%
  aov(Indexvalue~Plant*sex,data=.) %>% summary()

make_alpha_bar=function(inputdata,letterp,texty){
  signif_results<- auto_signif_test(data =inputdata,treatment_col =2,value_col =8,prior = T)
  letters<- data.frame(signif_results$comparison_letters,letterp=letterp)
  mean_frame<- aggregate(inputdata[,8],by=list(inputdata[,2]),FUN=mean)
  Sd<- aggregate(inputdata[,8],by=list(inputdata[,2]),FUN=sd) %>% .[,'x']
  Treatment_Name<- mean_frame$Group.1
  N<- table(inputdata[,2]) %>% as.numeric()
  Mean<- mean_frame[,'x']
  SEM<- Sd/(N^0.5)
  input_mean_frame<- data.frame(Treatment_Name,N,Mean,Sd,SEM)
  input_mean_frame$Treatment_Name=factor(input_mean_frame$Treatment_Name,levels = plant_order)
  ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
    theme_zg()+
    geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name),),color='#000000',alpha=0.8)+
    geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
    scale_y_continuous(expand = c(0.01,0.0))+
    scale_color_manual(values=plant_col)+
    scale_fill_manual(values=plant_col)+
    labs(x='',fill='',y=texty)+
    geom_quasirandom(data=inputdata,aes(x=Plant,y=Indexvalue,fill=Plant,shape=sex),alpha=0.8,size=1,show.legend = F,color="black")+
    scale_shape_manual(values = c(21,22))+
    geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
    theme(legend.position = 'right') %>% return()
}

alpha1=make_alpha_bar(alphaframe[alphaframe$Indexname=="Shannon",],5.5,"Shannon index")
alpha2=make_alpha_bar(alphaframe[alphaframe$Indexname=="Evenness",],0.7,"Evenness")

##outputs####
ggsave("./results/alpha/Shannon_bar.pdf",width = 5,height = 3)
alpha_results$plotlist$Plotobj_Chao$Barplot
ggsave("./results/alpha/Chao_bar.pdf",width = 5,height = 3)
alpha_results$plotlist$Plotobj_Evenness$Barplot
ggsave("./results/alpha/Evenness_bar.pdf",width = 5,height = 3)
alpha_results$plotlist$Plotobj_Shannon$Statistics$comparison_letters %>% write.table(file = "./results/alpha/Shannon_multiple_comp.txt",row.names = F,sep="\t",quote = T)
alpha_results$plotlist$Plotobj_Chao$Statistics$comparison_letters %>% write.table(file = "./results/alpha/Chao_multiple_comp.txt",row.names = F,sep="\t",quote = T)
alpha_results$plotlist$Plotobj_Evenness$Statistics$comparison_letters %>% write.table(file = "./results/alpha/Evenness_multiple_comp.txt",row.names = F,sep="\t",quote = T)


#make_sex_bar=function(Index){
#  inputdata<- subset(alphaframe,Indexname==Index)
#  for(i in unique(inputdata[,2])){
#    sub_facet<-inputdata[which((inputdata[,2])==i),] 
#    results<-auto_signif_test(data =sub_facet,treatment_col =3,value_col =8,prior = T)
#  }
#  mean_frame<- aggregate(inputdata[,8],by=list(inputdata[,3],inputdata[,2]),FUN=mean)
#  Sd<- aggregate(inputdata[,8],by=list(inputdata[,3],inputdata[,2]),FUN=sd)%>% .[,'x']
#  Treatment_Name<- mean_frame$Group.1
#  N<- table(inputdata[,3]) %>% as.numeric()
#  Mean<- mean_frame[,'x']
#  SEM<- Sd/(N^0.5)
#  input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])
#  colnames(input_mean_frame)[6]=colnames(inputdata)[2]
#  ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
#    theme_zg()+
#    geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name)),color='#000000',alpha=0.8)+
#    geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
#    scale_y_continuous(expand = c(0.00,0.00),limits = c(0,1.05*max(input_mean_frame$Mean+input_mean_frame$Sd)))+
#    scale_color_manual(values=sex_col)+
#    scale_fill_manual(values=sex_col)+
#    labs(x='',fill='',y=Index)+
#    facet_wrap(~get(colnames(inputdata)[2]),scales ='free_y',strip.position ='top',as.table =T,nrow=1)+
#    stat_compare_means(method = "t.test",data=inputdata,aes(x=as.factor(inputdata[,3]),y=inputdata[,8],label = paste0('p = ', after_stat(p.format))),size=4,label.x.npc = 'left',label.y.npc = 0.9)+
#    theme(legend.position = 'right',
#    ) %>% return()
#}
#
#shannon=make_sex_bar("Shannon")
#shannon
#ggsave("./results/alpha/F_FM_Shannon_bar.pdf",width = 10,height = 2)
#chao=make_sex_bar("Chao")
#chao
#ggsave("./results/alpha/F_FM_Chao_bar.pdf",width = 10,height = 2)
#
#Evenness=make_sex_bar("Evenness")
#Evenness
#ggsave("./results/alpha/F_FM_Evenness_bar.pdf",width = 10,height = 2)

#structure####
str_result=structure_plot(taxobj = tax_obj_plan2,taxlevel = "Genus")
str_result$PCoA_Plot
coord=str_result$PCoA_coordinates %>% .[,1:8]
str_result$PERMANOVA_statistics

library(pairwiseAdonis)
pairwise.adonis(t(tax_obj_plan1$Genus_percent[,-1]),metafile$Plant,sim.method = "bray")

cent = aggregate(coord[, 1:2], by = list(Plant=coord$Plant), FUN = mean)
colnames(cent)[2:3] = c( "cent1", "cent2")
PCoA_data = left_join(coord, cent) 


structure=ggplot(PCoA_data,aes(x=PCoA_data[,1],y=PCoA_data[,2],fill=Plant))+
  geom_point(color="black",aes(shape=sex))+
  geom_segment(aes(xend = cent1, yend = cent2, color = Plant), show.legend = FALSE, size = 0.3, alpha = 0.8)+
  #stat_ellipse(aes(color=Plant),level = 0.8)+
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(values = plant_col)+
  scale_color_manual(values = plant_col)+
  theme_zg()+labs(x="PCo1:34.47% ",y="Pco2:18.46%")
#ggsave("./results/structure/pcoa.pdf",width = 4.5,height = 3)
#community####  
com_result=community_plot(taxobj = tax_obj_plan2,taxlevel = "Phylum",n = 10,palette = "Paired",nrow = 1,rmprefix = "p__")
#com_result$barplot
#ggsave("./results/community/allsample_bar.pdf",width = 7,height = 5)
com_result=community_plot(taxobj = tax_obj_plan_plant,taxlevel = "Phylum",n = 10,palette = "Paired",nrow = 1,rmprefix = "p__")
community=com_result$mean_barplot
#ggsave("./results/community/mean_bar.pdf",width = 7,height = 3)


com_result_genus=community_plot(taxobj = tax_obj_plan2,taxlevel = "Genus",n = 15,nrow = 1,palette = "Set2",rmprefix = "g__")
com_result_genus$barplot
#ggsave("./results/community/allsample_bar_genus.pdf",width = 7,height = 5)
com_result_genus$mean_barplot
ggsave("./results/community/mean_bar_genus.pdf",width = 8,height = 3)

phylum_long=com_result$Grouped_Top10Phylum
phylum_stat=data.frame()
for(i in (unique(phylum_long$tax)%>% as.character())){
  subdata=phylum_long[phylum_long$tax==i,]
  stat=aov(rel_abundance~Plant*sex,data=subdata) %>% summary() %>% .[[1]] %>% as.data.frame()
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
write.table(phylum_stat,"./results/community/phylum_anova_stat.txt",row.names = F,sep="\t")

phylum_wide=com_result$Top10Phylum
colSums(phylum_wide[1:5,-1])
write.table(phylum_wide,"./results/community/phylum_rel.txt",row.names = F,sep="\t")
##Figure 2 top####
library(patchwork)
Fig2_top=(((alpha1|alpha2)|structure|community))*guides(fill="none",shape="none")+plot_layout(widths = c(1,1,1,2),guides = "collect")
ggsave("./Figures/Figure2/Figure2.pdf",width = 12,height = 2.5)
#core microbiome####

species=tax_obj_plan_plant$Species_percent
species_taxonomy=tax_obj_plan_plant$Species_taxonomy
species_trans=t(species[,-1])
colnames(species_trans)=species_taxonomy$SpeciesID
species_trans=data.frame(metafile,species_trans)
speices_mean_sex_plant=aggregate(species_trans[,-c(1:6)],by=list(Plant=species_trans$Plant,sex=species_trans$sex),FUN=mean)
speices_mean_plant=aggregate(species_trans[,-c(1:6)],by=list(Plant=species_trans$Plant),FUN=mean)

###exist species####
exist=function(x){
  which(x>0) %>% length() %>% return()
}
exist_species=apply(speices_mean_plant[,-1], 1, exist)
###high abundant species####
threshold=1e-4
high_abundant=function(x){
  which(x>threshold) %>% length() %>% return()
}
high_abundant_species=apply(speices_mean_plant[,-1], 1, high_abundant)

###core high abundant species####
core_high_species=intersect(
  colnames(speices_mean_plant)[which(speices_mean_plant[1,-1]>threshold)],
  colnames(speices_mean_plant)[which(speices_mean_plant[2,-1]>threshold)]
)
i=3
while(i<=8){
  core_high_species=intersect(
    core_high_species,
    colnames(speices_mean_plant)[which(speices_mean_plant[i,-1]>threshold)]
  )
  i=i+1
}

###flower visualizaiton####
core_stat=data.frame(Plant=speices_mean_plant$Plant,
           All=exist_species,
           High=high_abundant_species,
           core=length(core_high_species))
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
ggsave("./results/core/core_flower.pdf",width = 5,height = 4)

###core taxon composition####
library(stringr)
com_result_temp=community_plot(taxobj = tax_obj_plan_plant,taxlevel = "Phylum",n = 10,palette = "Paired",nrow = 1)
core_species_taxonomy=species_taxonomy[which(species_taxonomy$SpeciesID %in% core_high_species),]
core_phylum_stat=table(core_species_taxonomy$Phylum) %>% sort() %>% as.data.frame()
core_phylum_stat$Freq_pct=core_phylum_stat$Freq/sum(core_phylum_stat$Freq)
core_phylum_stat$Var1=as.character(core_phylum_stat$Var1)
core_phylum_stat$Var1[!core_phylum_stat$Var1 %in% names(com_result_temp$filled_color)]="Others"
core_phylum_stat_clean=aggregate(core_phylum_stat$Freq_pct,by=list(Phylum=core_phylum_stat$Var1),FUN=sum)
core_phylum_stat_clean$Phylum=factor(core_phylum_stat_clean$Phylum,levels = names(com_result_temp$filled_color))

ggplot(core_phylum_stat_clean, aes_string(x = "0", y = "x", fill = "Phylum")) +
  geom_bar(width = 1,  stat = "identity", show.legend = TRUE) + coord_polar("y",  start = 0) + 
  scale_fill_manual(values = com_result_temp$filled_color) + 
  labs(x = "", y = "", title = "Core taxon", fill = "Phylum") + 
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size=6),
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(hjust = 0.5,  size = 8))
ggsave("./results/core/core_pie.pdf",width = 4,height = 3)
write.table(core_phylum_stat,"./results/core/core_phylum_percent.txt",row.names = F,sep="\t")
##indicator####
set.seed(999)
indicator_results=indicator_analysis(taxobj =tax_obj_plan_plant,taxlevel = "Genus" )
library(dplyr)
library(tidyr)
library(stringr)
indicator_anno=subset(indicator_results,padj<0.05)
write.table(indicator_anno,"./results/indicator/indicator_annotation.txt",row.names = F,sep="\t")
#configuration

phylum_col=community_plot(taxobj = tax_obj_plan_plant,taxlevel = "Phylum",n = 10,palette = "Paired",nrow = 1) %$%
  filled_color

indicator_anno_long=gather(indicator_anno[,c(1:8,12)],"tag","value",-GenusID) %>% subset(.,value!=0)
indicator_anno_long=left_join(indicator_anno_long,tax_obj_plan_plant$Genus_taxonomy)
indicator_anno_long$tag=str_replace(indicator_anno_long$tag,"s.","")

###sort
indicator_anno_long1=data.frame()
for(i in unique(indicator_anno_long$tag)){
  for(j in names(phylum_col)){
    indicator_anno_long1=rbind(indicator_anno_long1,subset(indicator_anno_long,tag==i&Phylum==j))
  }
}

##add phylum color
for(i in unique(indicator_anno_long1$Phylum)){
  indicator_anno_long1$phylumcol[indicator_anno_long1$Phylum==i]=phylum_col[which(names(phylum_col)==i)]
}


#add relative abundance
i=1
while(i <=nrow(indicator_anno_long1)){
  ID=indicator_anno_long1$GenusID[i]
  plant=indicator_anno_long1$tag[i]
  Genus_percent=tax_obj_plan_plant$Genus_percent
  Genus_percent$Genus=tax_obj_plan_plant$Genus_taxonomy$GenusID
  indicator_anno_long1$rel[i]=
    Genus_percent[which(Genus_percent$Genus==ID),which(tax_obj_plan_plant$groupfile$Plant==plant)+1] %>%
    as.numeric() %>% mean()
  i=i+1
} #add relative abundance


#trans to numeric
for(i in unique(indicator_anno_long1$tag)){
  indicator_anno_long1$statvalue[indicator_anno_long1$tag==i]=1:length(which(indicator_anno_long1$tag==i))
}

##add specific tag
phylum_col=com_result$filled_color
indicator_anno_long1$sptag_col=(phylum_col[11] %>% as.character())
for(i in unique(indicator_anno_long1$tag)){
  indicator_anno_long1$sptag_col[indicator_anno_long1$tag==i]= as.character(plant_col[i])
}
##visualization
library(ggtreeExtra)
library(ggnewscale)
library(ggtree)
library(circlize)
pdf("./results/indicator/indicator.pdf",width = 6,height = 6)
circos.par(gap.degree=c(1,1,1,1,1,1,1,30))
circos.initialize(indicator_anno_long1$tag, x=indicator_anno_long1$statvalue) #initialize

circos.track(indicator_anno_long1$tag, ylim=c(0,1), bg.col =plant_col[unique(indicator_anno_long1$tag) %>% sort()],track.height=.12,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(-1), 
                           CELL_META$sector.index,niceFacing =F)
             }) #label layer
circos.track(indicator_anno_long1$tag, ylim=c(0,1),
             panel.fun = function(x, y) {
               #circos.text(CELL_META$xcenter, 
               #            CELL_META$cell.ylim[2] + mm_y(7), 
               #            CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             })  ##text layer

for(i in unique(indicator_anno_long1$tag)){
  range=indicator_anno_long1$statvalue[indicator_anno_long1$tag==i]
  breaks = seq(min(range), max(range), by = 1)
  n_breaks = length(breaks)
  
  circos.rect(sector.index = i,breaks[-n_breaks], rep(0, n_breaks - 1),
              breaks[-1], rep(1, n_breaks - 1),
              col = indicator_anno_long1$phylumcol[indicator_anno_long1$tag==i], border = NA)    
} 
circos.track(indicator_anno_long1$tag, ylim=c(0,6),
             track.height=.1,bg.border = "transparent",bg.col="gray95")  #"transparent"

for (i in unique(indicator_anno_long1$tag)){
  circos.barplot(sector.index = i,
                 value = log10(indicator_anno_long1[indicator_anno_long1$tag==i,]$rel)+6,
                 pos = indicator_anno_long1[indicator_anno_long1$tag==i,]$statvalue,
                 col = plant_col[i],
                 bar_width = 0.8,
                 border="transparent")
} #relative abundance bar
circos.yaxis(side =  "right",labels.cex = 0.4,labels = c(0,0.01,0.1))
circos.clear()
dev.off()





##network####
PG_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Apple")
HS_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Locust")
LS_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Willow")
YJ_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Rose")
SJ_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Seabuckthorn")
SHT_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Filbert")
XS_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Apricot")
YS_obj=sub_tax_summary(tax_obj_plan_plant,Plant=="Elm")

threshold=0.95
rel_threshold=5e-5


PG_obj_network=network_analysis(taxobj = PG_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)
HS_obj_network=network_analysis(taxobj = HS_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)
LS_obj_network=network_analysis(taxobj = LS_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)
YJ_obj_network=network_analysis(taxobj = YJ_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)
SJ_obj_network=network_analysis(taxobj = SJ_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)
XS_obj_network=network_analysis(taxobj = XS_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)
YS_obj_network=network_analysis(taxobj = YS_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)
SHT_obj_network=network_analysis(taxobj = SHT_obj,taxlevel = "Genus",reads = F,n = 6,threshold = threshold,rel_threshold = rel_threshold)

pdf("./results/network/network.pdf",width = 8,height =4 )
par(mfrow = c(2, 4), mar = c(0, 1, 1, 0), pty = "m")
HS_visual=network_visual(network_obj = HS_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])
LS_visual=network_visual(network_obj = LS_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])
PG_visual=network_visual(network_obj = PG_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])
SHT_visual=network_visual(network_obj = SHT_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])
SJ_visual=network_visual(network_obj = SJ_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])
XS_visual=network_visual(network_obj = XS_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])
YJ_visual=network_visual(network_obj = YJ_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])
YS_visual=network_visual(network_obj = YS_obj_network,mode = "major_tax",taxlevel = "Phylum",select_tax = names(com_result_temp$filled_color)[-11],palette = com_result_temp$filled_color[-11])

dev.off()

cat("####PG","\n")
cat(network_stat(PG_obj_network$Igraph_object))
cat("####HS","\n")
cat(network_stat(HS_obj_network$Igraph_object))
cat("####LS","\n")
cat(network_stat(LS_obj_network$Igraph_object))
cat("####YJ","\n")
cat(network_stat(YJ_obj_network$Igraph_object))
cat("####SJ","\n")
cat(network_stat(SJ_obj_network$Igraph_object))
cat("####SHT","\n")
cat(network_stat(SHT_obj_network$Igraph_object))
cat("####XS","\n")
cat(network_stat(XS_obj_network$Igraph_object))
cat("####YS","\n")
cat(network_stat(YS_obj_network$Igraph_object))

###
write.table(alphaframe,"Fig2a.txt",row.names = F,sep="\t")
write.table(PCoA_data,"Fig2b.txt",row.names = F,sep="\t")
write.table(com_result$Top10Phylum,"Fig2c.txt",row.names = F,sep="\t")
write.table(indicator_anno_long1,"Fig2d.txt",row.names = F,sep="\t")
