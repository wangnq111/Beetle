library(LorMe)
library(tidyr)
library(dplyr)
library(magrittr)
library(igraph)
load("./tax_obj_plan1.rda")
load("./metabolite.rda")
metabolite_pct1=metabolite_pct
rownames(metabolite_pct1)=metabolite$Compound_ID
genus_pct=tax_obj_plan2$Genus_percent


micro_metabo_newtwork=network_analysis2(input =genus_pct,inputtype = 2,n = 24,
                  threshold = 0.7,method = "spearman",display = T,
                  input2 =metabolite_pct1,input2type = 3  )


network_adj=micro_metabo_newtwork$Adjacency_column_table
network_nodes=micro_metabo_newtwork$Nodes_info
network_nodes$group="Microbe"
network_nodes$group[grep("Com_",network_nodes$nodes_id)]="Metabolites" 


nodes_phylum=data.frame(Genus=network_nodes$nodes_id[network_nodes$group=="Microbe"]) %>% left_join(tax_obj_plan2$Genus_taxonomy) %$% Phylum
main_phylum_col=community_plot(taxobj = tax_obj_plan1,taxlevel = "Phylum",n = 10,palette = "Paired",nrow = 1) %$%
  filled_color
main_phylum_name= main_phylum_col%>% names()
nodes_phylum[!nodes_phylum %in% main_phylum_name]="Others"
network_nodes$Category[network_nodes$group=="Microbe"]=nodes_phylum
nodes_class=data.frame(Compound_ID=network_nodes$nodes_id[network_nodes$group=="Metabolites"]) %>% left_join(metabolite[,c(1,14)]) %$% ClassI
nodes_class[!nodes_class %in% names(classfy_col)]="Others_metabo"
network_nodes$Category[network_nodes$group=="Metabolites"]=nodes_class
#col
network_nodes$col[network_nodes$group=="Microbe"]=main_phylum_col[nodes_phylum] %>% as.character()
network_nodes$col[network_nodes$group=="Metabolites"]=classfy_col[nodes_class] %>% as.character()
network_nodes$col[network_nodes$Category=="Others_metabo"]="#B09C85FF"
#size
network_nodes$size=case_when(
  network_nodes$node_degree==1 ~ 1,
  network_nodes$node_degree>1&network_nodes$node_degree<=3 ~ 3,
  network_nodes$node_degree>3&network_nodes$node_degree<=6 ~ 5,
  network_nodes$node_degree>6&network_nodes$node_degree<=12 ~ 7,
  network_nodes$node_degree>12 ~9
)
#shape
network_nodes$shape=case_when(
  network_nodes$group=="Microbe" ~"circle",
  network_nodes$group=="Metabolites" ~"square",
)
network_nodes$name=NA
network_nodes$name[network_nodes$node_degree>=5&network_nodes$group=="Microbe"]=network_nodes$nodes_id[network_nodes$node_degree>=5&network_nodes$group=="Microbe"]
key_name=data.frame(Name=metabolite$Name,row.names = metabolite$Compound_ID) %>% .[network_nodes$nodes_id[network_nodes$node_degree>=5&network_nodes$group=="Metabolites"],]
network_nodes$name[network_nodes$node_degree>=5&network_nodes$group=="Metabolites"]=key_name



#visualization
adjfile=network_adj
verticefile=network_nodes
t_net_its <- graph_from_data_frame(adjfile,direct=F,vertices = verticefile)
#color
V(t_net_its)$color <- V(t_net_its)$col

#size
V(t_net_its)$size <- V(t_net_its)$size

#shape
V(t_net_its)$shape=V(t_net_its)$shape

#linetype
E(t_net_its)$lty <- ifelse(E(t_net_its)$value == 1, 1, 2) 


set.seed(1)
coords_t_its <- layout_(t_net_its,with_fr(niter=9999, grid="auto"))
pdf("./micro_metabo_network.pdf",width = 12,height = 8)
plot(t_net_its,vertex.label=V(t_net_its)$name,
     vertex.size=V(t_net_its)$size,
     layout=coords_t_its,
     vertex.shape=V(t_net_its)$shape,
     edge.lty = E(t_net_its)$lty)



legend("topright", legend=c("Microbe","Metabolite"),
       pch=c(21,22), pt.bg="grey70", pt.cex=2, bty="n", title="Group")

# 2) 节点大小
legend("topleft",
       legend=c("Degree=1","3","5","10","15"),
       pt.cex=c(1,1.5,2,2.5,3),
       pch=21, pt.bg="grey70", bty="n", title="Degree")

# 3) 节点颜色 (Phylum/Category)
# 假设你有一个颜色向量 main_phylum_col
legend("bottomleft",
       legend = names(main_phylum_col),
       inset = c(0, -0.1),
       pch = 21,                # 圆形
       pt.bg = main_phylum_col, # 填充色
       col = "black",           # 边框色
       pt.cex = 1.8,
       bty = "n",
       title = "Phylum")
legend("bottomright",
       legend = names(classfy_col),
       inset = c(-0.3, 0.2),
       pch = 22,                # 圆形
       pt.bg = classfy_col, # 填充色
       col = "black",           # 边框色
       pt.cex = 1.8,
       bty = "n",
       title = "Phylum")

# 4) 边的线型 (相关性)
legend("bottomright",
       legend=c("Positive","Negative"),
       lty=c(1,2), col="black", bty="n", title="Correlation")
dev.off()

