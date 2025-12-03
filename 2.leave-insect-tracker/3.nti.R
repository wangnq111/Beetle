leave_asv=read.table("./input/leave_asv.txt",header = T,sep="\t")
load("./input/tax_obj.rda")
library(phyloseq)
phylogeny=read_tree("./input/phylogeny.tre")
library(stringr)
library(tidyr)
library(magrittr)
leave_asv$taxonomy= gsub(" ","",leave_asv$taxonomy)


meta_base=data.frame(tax_obj$Base,tax_obj$Base_taxonomy)
meta_base$taxonomy=gsub(" ","",meta_base$taxonomy)
metafile=tax_obj$Groupfile

leave_asv_new=data.frame()
meta_base_new=data.frame()
for(i in unique(leave_asv$taxonomy)){
  if(i!="Unclassified"){
    cat(i,"\n")
    leave_asv_sub=leave_asv[leave_asv$taxonomy %in% i,]
    max_asv=max(rowSums(leave_asv_sub[,2:49]))
    ASV_ID=leave_asv_sub$OTU.ID[which(rowSums(leave_asv_sub[,2:49])==max_asv)]%>% .[1]
    leave_asv_left=leave_asv_sub[1,]
    leave_asv_left[1,1]=ASV_ID
    leave_asv_left[,2:49]=colSums(leave_asv_sub[,2:49]) 
    leave_asv_left[1,50]=leave_asv_sub$taxonomy[leave_asv_sub$OTU.ID==ASV_ID]
    species=separate(leave_asv_left,col = "taxonomy",into=c(c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")),sep = ";") %$% Species
    meta_base_sub=meta_base[grep(species,meta_base$taxonomy),]
    if(nrow(meta_base_sub)==0){
      
    }else{
      meta_base_sub_left=meta_base_sub[1,]
      meta_base_sub_left[1,]=ASV_ID
      meta_base_sub_left[,2:49]=colSums(meta_base_sub[,2:49])
      meta_base_sub_left[1,50]=leave_asv_left[1,50]
      leave_asv_new=rbind(leave_asv_new,leave_asv_left)
      meta_base_new=rbind(meta_base_new,meta_base_sub_left)
    }
  }
  
}
library(NST)
library(picante)
for(i in unique(metafile$Plant)){ 
  for(j in unique(metafile$sex)){
    otutable=data.frame(meta_base_new[,which(metafile$Plant==i&metafile$sex==j)+1],
                        leave_asv_new[,which(metafile$Plant==i&metafile$sex==j)+1])
    otut=data.frame(t(otutable))
    colnames(otut)=meta_base_new$ID
    group=data.frame(row.names = colnames(otutable),group=c(rep("Beetle",3),rep("Leave",3)))
    tree = prune.sample(otut, phylogeny)
    cat("Start",i,"_",j,"\n")
    set.seed(123)
    pd.wd_path=paste0("./nti/",i)
    pnst <- pNST(comm = otut, tree = tree, group = group, phylo.shuffle = TRUE, taxo.null.model = NULL, 
                 pd.wd = pd.wd_path, abundance.weighted = TRUE, rand = 99, nworker = 2, SES = T, RC = TRUE)
    temp=pnst$index.pair
    tempID=paste0(i,"_",j)
    temp$Plant=i
    temp$sex=j
    #rmnum=intersect(grep("CK",temp$name1),grep("CK",temp$name2))
    #rmnum=c(rmnum,intersect(grep(tempID,temp$name1),grep(tempID,temp$name2)))
    #temp=temp[-rmnum,]
    savetxt=paste0("./allntitxt/pnst_",i,"_",j,".txt")
    write.table(temp,savetxt,row.names = F,sep="\t")
    savename=paste0("./allpnst/pnst_",i,"_",j,".rda")
    save(pnst,file=savename)
  }
}
####
allnti=read.table("./allntitxt/all.txt",header = T,sep="\t")


allnti$Letters=case_when(
  (allnti$bNTI.wt)>2 ~ "Homogeneous",
  (allnti$bNTI.wt) <(-2) ~ "Heterogeneous",
  (allnti$bNTI.wt) <=(2)&allnti$bNTI.wt >=(-2)&allnti$RC.bMNTD.wt>0.95~ "Dispersal",
  (allnti$bNTI.wt) <=(2)&allnti$bNTI.wt >=(-2)&allnti$RC.bMNTD.wt<0.95~ "Drift"
  #.default = as.character(dunnet_stat$pvalue)
)

nti_stat=table(allnti[,c(15,13)]) %>% as.data.frame()
nti_stat$Proportion=nti_stat$Freq/18
#nti_stat$combine=factor(paste0(nti_stat$Plant,nti_stat$sex),levels = )

nti_stat$Plant=factor(nti_stat$Plant,levels = plant_order)

p3=ggplot(nti_stat,aes(x=Plant,y=Proportion))+
  geom_col(position = "stack", width = 0.8,aes(fill=Letters))+
  #facet_grid(~Plant)+
  scale_fill_manual(values=color_scheme("Plan2") %>% .[c(3,1,2)])+
  scale_y_continuous(labels = scales::percent, expand = c(0, 0))+
  labs(fill="",y="Relative Proportion",x="")+
  theme_zg()
ggsave("./nti.pdf",width = 5,height = 2.5)


aov(bNTI.wt~Plant*sex,data=allnti) %>% summary()

allnti$Plant=factor(allnti$Plant,levels = plant_order)
signif_results<- auto_signif_test(data =allnti,treatment_col =13,value_col =11,prior = T)
letters<- data.frame(signif_results$comparison_letters,letterp=8)
p1=ggplot(allnti,aes(x=as.factor(allnti[,13]),y=allnti[,11]))+
  scale_y_continuous(expand = c(0.05,0.01),limits = c(-5,6))+
  scale_color_manual(values=plant_col)+
  scale_fill_manual(values=plant_col)+
  geom_hline(yintercept = 2,lty=2)+
  geom_hline(yintercept = -2,lty=2)+
  theme_zg()+
  geom_violin(aes(color=factor(allnti[,13])),alpha=0.8,width=0.8,linewidth=0.5,trim = F)+
  labs(x='',color='',y="βNTI")+
  geom_quasirandom(aes(color=as.factor(allnti[,13])),alpha=0.8,size=1,pch=16,show.legend = F)+
  #geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
  theme(legend.position = 'right')+guides(fill="none",color="none") 
ggsave("./ntidetail.pdf",width = 4,height = 3)

allnti$sex=factor(allnti$sex,levels = c("F","FM"))
p2=ggplot(allnti,aes(x=as.factor(allnti[,14]),y=allnti[,11]))+
  scale_y_continuous(expand = c(0.05,0.01),limits = c(-5,6))+
  scale_color_manual(values=sex_col)+
  scale_fill_manual(values=sex_col)+
  geom_hline(yintercept = 2,lty=2)+
  geom_hline(yintercept = -2,lty=2)+
  theme_zg()+
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "wilcox.test")+
  geom_violin(aes(color=factor(allnti[,14])),alpha=0.8,width=0.8,linewidth=0.5,trim = F)+
  labs(x='',color='',y="βNTI")+
  geom_quasirandom(aes(color=as.factor(allnti[,14])),alpha=0.8,size=1,pch=16,show.legend = F)+
  theme(legend.position = 'right')+guides(fill="none",color="none") 


(p1|p2|p3)*guides(fill="none")+plot_layout(widths = c(2,1,2))
ggsave("./figure4.pdf",width = 8,height = 3)


write.table(allnti,"Fig3ab_sourcedata.txt",row.names = F,sep="\t")
write.table(nti_stat,"Fig3c_sourcedata.txt",row.names = F,sep="\t")
