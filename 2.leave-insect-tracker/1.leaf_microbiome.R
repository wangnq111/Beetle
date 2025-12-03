library(LorMe)
library(patchwork)
library(ggplot2)
library(magrittr)
library(phyloseq)
leave_asv=read.table("./input/leave_asv.txt",header = T,sep="\t")
load("./input/tax_obj.rda")
phylogeny=read_tree("./input/phylogeny.tre")
metafile=read.table("./input/metafile.txt",header = T,sep="\t")
plant_col=color_scheme("Plan7",names=unique(metafile$Plant))
plant_order=c("Locust","Rose","Apricot","Seabuckthorn","Willow","Apple","Elm","Filbert")
sex_col=c("F"="#FE5D5D",FM="#6376A0")


leave_obj=tax_summary(groupfile = metafile,inputtable = leave_asv[,2:49],reads = T,taxonomytable = leave_asv[,c(1,50)],outputtax = "standard")
leave_obj_plan=object_config(taxobj = leave_obj,treat_location =2,rep_location =6 ,treat_col = plant_col,treat_order = plant_order)
save(leave_obj_plan,file="leave_obj_plan.rda")
leave_alpha=Alpha_diversity_calculator(taxobj = leave_obj_plan,taxlevel = "Genus")
p1=leave_alpha$plotlist$Plotobj_Shannon$Barplot+guides(fill="none")
p2=leave_alpha$plotlist$Plotobj_Evenness$Barplot+guides(fill="none")
alpha_frame=leave_alpha$alphaframe
aov(Indexvalue~Plant,data=alpha_frame[alpha_frame$Indexname=="Shannon",]) %>% summary()
aov(Indexvalue~Plant,data=alpha_frame[alpha_frame$Indexname=="Evenness",]) %>% summary()

leave_stucurure=structure_plot(taxobj =leave_obj_plan,taxlevel = "Genus",diagram = "stick")
p3=leave_stucurure$PCoA_Plot
leave_stucurure$PERMANOVA_statistics

leave_com=community_plot(taxobj =leave_obj_plan,taxlevel = "Phylum",n = 10)
p4=leave_com$mean_barplot

((p1|p2)/(p3|p4))+plot_layout(guides = "collect")
ggsave("./leave_profile.pdf",width = 8,height = 6)

library(dplyr)
phylum_long=leave_com$Grouped_Top10Phylum
phylum_stat=data.frame()
for(i in (unique(phylum_long$tax)%>% as.character())){
  subdata=phylum_long[phylum_long$tax==i,]
  stat=aov(rel_abundance~Plant,data=subdata) %>% summary() %>% .[[1]] %>% as.data.frame()
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
#write.table(phylum_stat,"./results/community/phylum_anova_stat.txt",row.names = F,sep="\t")
phylum_wide=leave_com$Top10Phylum
colSums(phylum_wide[1:5,-1])
