library(LorMe)
library(magrittr)
library(dplyr)
library(ggbeeswarm)
library(stringr)

load("leave_obj_plan.rda")
load("tax_obj_plan1.rda")
metafile=read.table("./input/metafile.txt",header = T,sep="\t")
plant_col=color_scheme("Plan7",names=unique(metafile$Plant))

leave_genus=leave_obj_plan$Genus
leave_genus$Genus=str_replace(leave_genus$Genus," ","")
meta_genus=tax_obj_plan1$Genus

leave_genus_new=data.frame()
meta_genus_new=data.frame()
for(i in unique(leave_genus$Genus)){
  if(i!="Unclassified"){
    cat(i,"\n")
    leave_genus_sub=leave_genus[leave_genus$Genus %in% i,]
    max_asv=max(rowSums(leave_genus_sub[,2:49]))
    Genus_ID=leave_genus_sub$Genus[which(rowSums(leave_genus_sub[,2:49])==max_asv)]%>% .[1]
    leave_genus_left=leave_genus_sub[1,]
    leave_genus_left[1,1]=Genus_ID
    leave_genus_left[,2:49]=colSums(leave_genus_sub[,2:49]) 
    leave_genus_left[1,50]=genus=leave_genus_sub$Genus[leave_genus_sub$Genus==Genus_ID]
    meta_genus_sub=meta_genus[grep(genus,meta_genus$Genus),]
    if(nrow(meta_genus_sub)==0){
    }else{
      meta_genus_sub_left=meta_genus_sub[1,]
      meta_genus_sub_left[1,]=Genus_ID
      meta_genus_sub_left[,2:49]=colSums(meta_genus_sub[,2:49])
      meta_genus_sub_left[1,50]=leave_genus_left[1,50]
      leave_genus_new=rbind(leave_genus_new,leave_genus_left)
      meta_genus_new=rbind(meta_genus_new,meta_genus_sub_left)
    }
  }
  
}

library(FEAST)

# default 1000
EM_iterations = 1000

# if you use different sources for each sink, different_sources_flag = 1, otherwise = 0
different_sources_flag = 0

data_prop=data.frame()
for(i in unique(metafile$Plant)){
  sub_metafile=data.frame(#SampleID=c(metafile$Sample[metafile$Plant==i],colnames(leave_genus_new)[which(metafile$Plant==i)+1]),
    Env=c(metafile$sex[metafile$Plant==i],rep("leaves",6)),
    SourceSink=c(rep("Sink",6),rep("Source",6)),
    id=c(metafile$Sample[metafile$Plant==i],rep(NA,6)),row.names = c(metafile$Sample[metafile$Plant==i],colnames(leave_genus_new)[which(metafile$Plant==i)+1]))
  asv_sub=data.frame(otuid=leave_genus_new$Genus,meta_genus_new[,which(metafile$Plant==i)+1]/30,leave_genus_new[,which(metafile$Plant==i)+1])
  #asv_sub[,-1]=apply(asv_sub[,-1],1,round)
  otus=t(asv_sub[,-1])
  otus=apply(otus,1,round)
  rownames(otus)=asv_sub$otuid
  for(sample in 1:6){
    FFEAST_output <- FEAST( t(otus)[c(sample,7:12),], metadata = sub_metafile[c(sample,7:12),],EM_iterations = EM_iterations,
                            different_sources_flag = different_sources_flag,dir_path = "./",outfile = "source")
    data_prop_temp=FFEAST_output$data_prop %>% .[1:6,]
    data_prop_temp$Sample=rownames(sub_metafile)[sample] 
    data_prop_temp= left_join(data_prop_temp,metafile)
    data_prop=rbind(data_prop,data_prop_temp)
  }
}
temp=auto_signif_test(data = data_prop,treatment_col =4,value_col = 1,prior = T)

aov(pred_emnoise_all~Plant*sex,data=data_prop) %>% summary()


data_prop$Plant=factor(data_prop$Plant,levels = plant_order)
signif_results<- auto_signif_test(data =data_prop,treatment_col =4,value_col =1,prior = T)
letters<- data.frame(signif_results$comparison_letters,letterp=0.73)
ggplot(data_prop,aes(x=as.factor(data_prop[,4]),y=data_prop[,1]))+
  scale_y_continuous(expand = c(0.05,0.01),labels = scales::label_percent())+
  scale_color_manual(values=plant_col)+
  scale_fill_manual(values=plant_col)+
  theme_zg()+
  geom_violin(aes(color=factor(data_prop[,4])),alpha=0.8,width=0.8,linewidth=0.5,trim = F)+
  labs(x='',color='',y="Proportions from leaves")+
  geom_quasirandom(aes(color=as.factor(data_prop[,4])),alpha=0.8,size=1,pch=16,show.legend = F)+
  geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
  theme(legend.position = 'right')+guides(fill="none",color="none") 
#ggsave("./source_proportion.pdf",width = 4,height = 3)
write.table(data_prop,"Fig3de_sourcedata.txt",row.names = F,sep="\t")

