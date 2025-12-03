library(LorMe)
load("./beetle_community.rda")
load("./beetle_metabolome.rda")
load("./leave_obj_plan.rda")

##leave
leave_com=community_plot(taxobj = leave_obj_plan,taxlevel = "Phylum",n = 10)
leave_com_top_phylum=leave_com$Top10Phylum
sem_data=data.frame(metafile,t(leave_com_top_phylum[,-1]))
colnames(sem_data)[9:19]=paste0("leaf_",leave_com_top_phylum$Phylum)

leave_str=structure_plot(taxobj = leave_obj_plan,taxlevel = "Genus")

sem_data=data.frame(sem_data,leave_str$PCoA_coordinates[,1:2])
colnames(sem_data)[20:21]=c("leaf_PC1","leaf_PC2")

##beetle community
beetle_com=community_plot(taxobj =tax_obj_plan2 ,taxlevel = "Phylum",n = 10)
beetle_com_top_phylum=beetle_com$Top10Phylum
sem_data=data.frame(sem_data,t(beetle_com_top_phylum[,-1]))
colnames(sem_data)[22:32]=paste0("beetle_",beetle_com_top_phylum$Phylum)

beetle_str=structure_plot(taxobj = tax_obj_plan2,taxlevel = "Genus")
sem_data=data.frame(sem_data,beetle_str$PCoA_coordinates[,1:2])
colnames(sem_data)[33:34]=c("beetle_PC1","beetle_PC2")

##beetle cazy
sem_data=data.frame(sem_data,t(cate2_mean[,-1]))
colnames(sem_data)[35:43]=paste0("CAZy_",cate2_mean$Type)


sem_data=data.frame(sem_data,PCAframe[,1:2])
colnames(sem_data)[44:45]=c("CAZy_PC1","CAZy_PC2")

##beetle metabolome
sem_data=data.frame(sem_data,t(metabolite_summary[,-1]))
colnames(sem_data)[46:55]=paste0("Kegg_",metabolite_summary$Group.1)

sem_data=data.frame(sem_data,plsdaScore[,3:4])
colnames(sem_data)[56:57]=c("Kegg_PC1","Kegg_PC2")

temp1=circulation_lm(y = sem_data$Length,sem_data[,7:57],margin = 2) 
temp2=circulation_lm(y = sem_data$Weight,sem_data[,7:57],margin = 2) 
select_factor=c(temp1$ID[temp1$pvalue<0.05],temp2$ID[temp2$pvalue<0.05]) %>% unique()
##sem

library(lavaan)
library(semPlot)

df=scale(sem_data[,7:57],center=TRUE,scale=TRUE)
df[,"beetle_PC2"]=-(df[,"beetle_PC2"] %>% as.numeric())
model <- '
# 1.direct effect
Weight ~  beetle_p__Patescibacteria + Kegg_Organoheterocyclic.compounds + 
         Kegg_Organic.nitrogen.compounds

Length ~ beetle_p__Patescibacteria + Kegg_Organoheterocyclic.compounds  +beetle_PC2

# leaf -> beetle taxa
beetle_p__Patescibacteria~leaf_PC2
beetle_PC2~leaf_PC2
beetle_p__Spirochaetota~leaf_p__Cyanobacteriota

# leaf direct effect on CAZy



# beetle taxa -> CAZy 
#CAZy_PC1~beetle_p__Patescibacteria+beetle_p__Spirochaetota
CAZy_SLH~beetle_PC2+beetle_p__Spirochaetota

# CAZy -> Kegg
Kegg_Organoheterocyclic.compounds~CAZy_SLH+ beetle_p__Patescibacteria
Kegg_Organic.nitrogen.compounds~leaf_PC2
 
# 残差相关
Length ~~ Weight
beetle_p__Patescibacteria ~~beetle_p__Spirochaetota
beetle_p__Patescibacteria ~~beetle_PC2
beetle_p__Spirochaetota ~~beetle_PC2
leaf_p__Cyanobacteriota~~leaf_PC2
Kegg_Organoheterocyclic.compounds ~~   Kegg_Organic.nitrogen.compounds
'

fit_model<-sem(model,data=df)
#resid(fit_model, type = "cor")
#查看模型结果以及路径显著性
summary(fit_model,fit.measures=TRUE,standardized=TRUE,rsq=T)
fitmeasures(fit_model,c("chisq", "df", "gfi", "rmsea", "srmr","pvalue"))
##参数标准，检验模型是否合理,一般标准为
#df<5;gfi>0.8;rmsea<0.05;srmr<0.08
#具体标准可自行参考文献
modindices(fit_model,sort.=TRUE)
#查看模型系数
inspect(fit_model,what="std")
#查看SEM模型参数，用于评估可靠性


#对模型结果进行简单可视化
set.seed(333)
semPaths(fit_model,"std",edge.label.cex=0.8,#字体
         fade=F, layout = "spring",#布局
         optimizeLatRes = FALSE, residuals = FALSE)

