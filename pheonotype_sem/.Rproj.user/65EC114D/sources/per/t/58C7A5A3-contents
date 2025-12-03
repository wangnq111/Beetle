leave_genus=leave_obj_plan$Genus_percent
leave_genus_taxonomy=leave_obj_plan$Genus_taxonomy

leave_genus_taxonomy$Phylum[!leave_genus_taxonomy$Phylum %in% leave_com_top_phylum$Phylum[1:10] ] ="Other"

beetle_genus=tax_obj_plan2$Genus_percent
beetle_genus_taxonomy=tax_obj_plan2$Genus_taxonomy

beetle_genus_taxonomy$Phylum[!beetle_com_top_phylum$Phylum %in% beetle_com_top_phylum$Phylum[1:10] ] ="Other"

CAZY=read.table("./input//CAZY.txt",header = T,sep="\t")
CAZY_pct=sweep(CAZY[,-c(1:2)],colSums(CAZY[,-c(1:2)]),"/",MARGIN=2)

metabolite=read.csv("./input/metabolism.csv",header = T,encoding = "UTF-8",as.is = TRUE)
metabolite_pct=sweep(metabolite[,35:82],colSums(metabolite[,35:82]),"/",MARGIN=2)
kegg_anno=data.frame(KeggID=paste0("kegg",1:nrow(metabolite)),metabolite[,1:34])

combine_frame=data.frame(t(leave_genus[,-1]),t(beetle_genus[,-1]),t(CAZY_pct),t(metabolite_pct))
colnames(combine_frame)=c(paste0("leaf_",leave_genus_taxonomy$GenusID),
                          paste0("beetle_",beetle_genus_taxonomy$GenusID),
                          paste0("cazy_",CAZY$Subtype),
                          paste0("kegg_",kegg_anno$KeggID))


combine_anno=data.frame(ID=colnames(combine_frame),
                        cate=c(paste0("leaf_",leave_genus_taxonomy$Phylum),
                        paste0("beetle_",beetle_genus_taxonomy$Phylum),
                        paste0("cazy_",CAZY$Type),
                        paste0("kegg_",kegg_anno$ClassI)))



cate_name=c(paste0("leaf_",leave_com_top_phylum$Phylum[1:10]),
            paste0("beetle_",beetle_com_top_phylum$Phylum[1:10]),
            paste0("cazy_",unique(CAZY$Type)),
            paste0("kegg_",unique(kegg_anno$ClassI)))
cate_name=cate_name[!cate_name %in% c("kegg__")]
phylumlist=list()
for(i in cate_name){
  num=grep(i,combine_anno$cate)
  phylumlist=c(phylumlist,list(num))
  names(phylumlist)[[length(phylumlist)]]=i
}


library(linkET)
#mantel_results <- mantel_test(combine_frame, metafile[,c(7,8)], 
#                              spec_select = phylumlist,mantel_fun = "mantel") 
#r_frame=data.frame()
#p_frame=data.frame()
#for(i in names(phylumlist)){
#  r_list=c()
#  p_list=c()
#  for(j in names(phylumlist)){
#    matrix1=(combine_frame[,phylumlist[[i]]])
#    matrix2=(combine_frame[,phylumlist[[j]]])
#    
#    veg.dist <- vegdist(matrix1) # Bray-Curtis
#    env.dist <- vegdist(matrix2)
#    cor=mantel(veg.dist, env.dist, method="spear",na.rm = T)
#    r=cor$statistic
#    r_list=c(r_list,r)
#    p=cor$signif
#    p_list=c(p_list,p)
#    
#  }
#  if(sum(dim(r_frame)==0)){
#    r_frame=rbind(r_frame,data.frame(r_list))
#    colnames(r_frame)[ncol(r_frame)]=i
#    p_frame=rbind(p_frame,data.frame(p_list))
#    colnames(p_frame)[ncol(p_frame)]=i
#  }else{
#    r_frame=data.frame(r_frame,data.frame(r_list))
#    colnames(r_frame)[ncol(r_frame)]=i
#    p_frame=data.frame(p_frame,data.frame(p_list))
#    colnames(p_frame)[ncol(p_frame)]=i
#  }
# 
#}
#rownames(r_frame)=rownames(p_frame)=names(phylumlist)
#r_frame_copy=r_frame
#p_frame_copy=p_frame



###
select_factor_frame=sem_data[,select_factor[-c(1:2,9,12,14,15,18)]]


cor_obj=correlate(sem_data[,7:8],select_factor_frame,method="spearman")
cor_result=data.frame(cor_obj$r,spec = rownames(cor_obj$r))%>% gather("env","r",-spec)
cor_result$p=data.frame(cor_obj$p,spec = rownames(cor_obj$p))%>% gather("env","p",-spec) %$% p
cor_result$r[cor_result$p>0.05]=0
cor_result_clean=dplyr::mutate(cor_result,rd = cut(r, breaks = c(-Inf, 0.3, 0.4, Inf), 
                                      labels = c("Low", "Mid", "High")),
              pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf), 
                       labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
              effect=cut(r, breaks = c(-Inf, 0, Inf), 
                         labels = c("Neg", "Pos")))
library(ggplot2)
qcorrplot(correlate(select_factor_frame,method="spearman"),type = "lower", diag = FALSE ) +
  geom_square() +
  geom_couple(aes(colour = effect, size = rd), 
              data = cor_result_clean, 
              curvature = nice_curvature())+
  geom_mark(sep = '\n',size = 3, sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05, color = 'white') + ##添加显著性和相关性值
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#B2182B","#2166AC")) +
  guides(size = guide_legend(title = "Spearman's rho",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Effect",
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Spearman's rho", order = 3))
ggsave("allcor.pdf",width = 12,height = 8)


