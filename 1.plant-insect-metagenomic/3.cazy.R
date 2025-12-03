library(factoextra)
library(FactoMineR)
library(magrittr)
library(dplyr)
library(LorMe)
library(vegan)
library(patchwork)

CAZY=read.table("./data/CAZY.txt",header = T,sep="\t")
CAZY_cato=data.frame(CAZYID=paste0("CAZY",1:nrow(CAZY)),CAZY[,1:2])

CAZY_pct=sweep(CAZY[,-c(1:2)],colSums(CAZY[,-c(1:2)]),"/",MARGIN=2)
metafile=read.table("./data/metafile.txt",header = T,sep="\t")
plant_col=color_scheme("Plan7",names=unique(metafile$Plant))
plant_order=c("Locust","Rose","Apricot","Seabuckthorn","Willow","Apple","Elm","Filbert")

res.pca <- PCA(t(CAZY_pct), 
               scale.unit = F,
               ncp = 11, 
               #ind.sup = 1, 
               quanti.sup = NULL,
               quali.sup = NULL,
               row.w = NULL, col.w = NULL, 
               graph = FALSE)

sam=data.frame(res.pca$ind$coord) %>% .[1:2] %>% data.frame(metafile)
cent = aggregate(sam[, 1:2], by = list(sam$Plant),FUN = mean)
colnames(cent) = c("Plant", "cent1", "cent2")
PCAframe = left_join(data.frame(sam, joint = sam$Plant), 
                     data.frame(cent[, 2:3], joint = cent$Plant)) %>% suppressMessages()

ggplot(PCAframe,aes(x=Dim.1,y=Dim.2,fill=factor(Plant)))+
  geom_point(size=1.5,alpha=.8,aes(shape=sex))+
  geom_hline(yintercept=0,linetype=4,color="grey") +
  geom_vline(xintercept=0,linetype=4,color="grey") +
  geom_segment(aes(xend = cent1, yend = cent2,color=Plant),show.legend = FALSE, size = 0.3, alpha = 0.6)+
  scale_shape_manual(values = c(21,22))+
  scale_color_manual(values = plant_col)+
  scale_fill_manual(values = plant_col)+
  labs(x=paste0("PC1:",get_eigenvalue(res.pca) %>% .[1,2]%>% round(.,2),"%"),
       y=paste0("PC1:",get_eigenvalue(res.pca) %>% .[2,2]%>% round(.,2),"%"),
       fill="")+
  # stat_ellipse(level=0.95)+
  theme_zg()+guides(color="none")#+
#theme(legend.position = c(.85,.7),legend.background = element_rect(color="black",linetype = 2,linewidth = .5))
adonis2(t(CAZY_pct)~metafile$Plant*metafile$sex,method="euclidean",by = "terms")
#ggsave("./results/CAZY_function/function_pcoa.pdf",width = 4.5,height = 2.8)

###community####
cate2_mean=aggregate(CAZY_pct,by=list(Type=CAZY_cato$Type),FUN=sum)
cate2_mean_long=gather(cate2_mean,"Sample","Rel",-c(Type))
cate2_mean_long=left_join(cate2_mean_long,metafile)
cate2_mean_long_mean=aggregate(cate2_mean_long$Rel,by=list(Type=cate2_mean_long$Type,Plant=cate2_mean_long$Plant),FUN=mean)
cate2_mean_long_mean$Plant=factor(cate2_mean_long_mean$Plant,levels = plant_order)
cate2_mean_long_mean$Type=factor(cate2_mean_long_mean$Type,levels =cate2_mean$Type[rowMeans(cate2_mean[,-1]) %>% order(decreasing = T)] )

cate2_mean[cate2_mean$Type %in% c("GH","GT","CE","CBM","PL"),-1] %>% colSums() %>% sort()


cazy_col=color_scheme("Plan9") %>% .[1:9]
names(cazy_col)=levels(cate2_mean_long_mean$Type)
ggplot(cate2_mean_long_mean,aes(x = Plant, y = x, fill = Type)) + 
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values =cazy_col ) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + 
  xlab("") + ylab("Relative abundance") + 
  labs(fill = "CAZy type") +
  #facet_wrap(~Plant, nrow = 2, scales = "free_x")+
  theme_zg() + theme(
    #axis.text.x = element_blank(),
    legend.text = element_text(size = 8, face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(color = NA),
    axis.line.y = element_line(size = 0.4, color = "black"),
    axis.ticks.length.y = unit(0.4, "lines"),
    axis.ticks.y = element_line(color = "black", size = 0.4)
  )
ggsave("./results/CAZY_function/composition_cazy.pdf",width = 5,height = 3)


cate2_stat=data.frame()
for(i in (unique(cate2_mean$Type)%>% as.character())){
  subdata=cate2_mean[cate2_mean$Type==i,-1] %>% as.numeric()
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


write.table(cate2_stat,"./results/CAZY_function/CAZY_anova.txt",row.names = F,sep="\t")


###core cazy####
cazy_mean_plant=aggregate(t(CAZY_pct),by=list(Plant=metafile$Plant),FUN=mean)
###exist species####
exist=function(x){
  which(x>0) %>% length() %>% return()
}
exist_cazy=apply(cazy_mean_plant[,-1], 1, exist)


###high abundant cazy####
threshold=1e-4
high_abundant=function(x){
  which(x>threshold) %>% length() %>% return()
}
high_abundant_cazy=apply(cazy_mean_plant[,-1], 1, high_abundant)

###core high abundant cazy####
core_high_cazy=intersect(
  colnames(cazy_mean_plant)[which(cazy_mean_plant[1,-1]>threshold)+1],
  colnames(cazy_mean_plant)[which(cazy_mean_plant[2,-1]>threshold)+1]
)
i=3
while(i<=8){
  core_high_cazy=intersect(
    core_high_cazy,
    colnames(cazy_mean_plant)[which(cazy_mean_plant[i,-1]>threshold)+1]
  )
  i=i+1
}

###flower visualizaiton####
core_stat=data.frame(Plant=cazy_mean_plant$Plant,
                     All=exist_cazy,
                     High=high_abundant_cazy,
                     core=length(core_high_cazy))
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
#ggsave("./results/CAZY_function//core_flower_cazy.pdf",width = 5,height = 4)

###indicator####
indic_data=as.data.frame(t(CAZY_pct))
colnames(indic_data)=paste0("CAZY",1:ncol(indic_data))

CAZY_indicator=indicspecies::multipatt(indic_data, 
                                       metafile$Plant, func = "r.g", control = permute::how(nperm = 1000)) %$% 
  as.data.frame(sign)
CAZY_indicator=data.frame(CAZYID=rownames(CAZY_indicator),CAZY_indicator)
taxonomy=data.frame(CAZYID=paste0("CAZY",1:ncol(indic_data))) %>%
  left_join(CAZY_cato)
CAZY_indicator_anno = left_join(CAZY_indicator, taxonomy)

pie_list=list()
for (i in unique(metafile$Plant)) {
  CAZY_indicator_anno_select=
    CAZY_indicator_anno[CAZY_indicator_anno[,paste0("s.", i) ]== 1 & CAZY_indicator_anno$p.value <  0.05&CAZY_indicator_anno$stat>0.5,]
  type_frame=table(CAZY_indicator_anno_select$Type) %>% sort(decreasing = TRUE) %>% as.data.frame()
  type_frame$Var1=factor(type_frame$Var1,levels = levels(cate2_mean_long_mean$Type))
  p=ggplot(type_frame, aes_string(x = "0",y = "Freq", fill = "Var1")) + 
    geom_bar(width = 1,  stat = "identity", show.legend = TRUE) + 
    coord_polar("y", start = 0) + 
    scale_fill_manual(values = color_scheme("Plan9",names = levels(cate2_mean_long_mean$Type))) + 
    labs(x = "", y = "", title = i, fill = "Type") + 
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), 
          axis.text = element_blank(),  
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 8))
  pie_list=c(pie_list,list(p))
  names(pie_list)[length(pie_list)]=i
  }

up_pie=(pie_list$Locust|pie_list$Rose|pie_list$Apricot|pie_list$Seabuckthorn)*guides(fill="none")
down_pie=(pie_list$Willow|pie_list$Apple|pie_list$Elm|pie_list$Filbert)*guides(fill="none")

up_pie/down_pie
#ggsave("./results/CAZY_function/indicator_pie.pdf",width = 8,height = 4)

cate2_select=cate2_stat[cate2_stat$`Pr(>F)`<0.05&cate2_stat$Item=="Plant      ",]
cate2_select=cate2_select[order(cate2_select$`Pr(>F)`,decreasing = F),]
select_cate=CAZY_indicator_anno$CAZYID[CAZY_indicator_anno$p.value<0.05]
########CAZY circos#############

cate3=data.frame(Subtype=CAZY$Subtype,CAZY_pct)
cate3_stat=data.frame()
for(i in (unique(cate3$Subtype)%>% as.character())){
  subdata=cate3[cate3$Subtype==i,-1] %>% as.numeric()
  stat=aov(subdata~Plant*sex,data=metafile) %>% summary() %>% .[[1]] %>% as.data.frame()
  stat$Item=rownames(stat)
  stat$signif=case_when(
    stat$`Pr(>F)`>0.05 ~"ns",
    stat$`Pr(>F)`>0.01&stat$`Pr(>F)`<0.05 ~"*",
    stat$`Pr(>F)`>0.001&stat$`Pr(>F)`<0.01 ~"**",
    stat$`Pr(>F)`<0.001 ~"***",
    is.na( stat$`Pr(>F)`)==T~NA
  )
  stat$Subtype=i
  cate3_stat=rbind(cate3_stat,stat)
}

cate3_stat=left_join(cate3_stat,CAZY[,1:2])

##All Cazy subtype####
#CAZY_pct_order=CAZY_pct[order(rowMeans(CAZY_pct),decreasing = T),] #%>% .[1:20,]
CAZY_anno=CAZY[,1:2]
#CAZY_anno_order=CAZY_anno[order(rowMeans(CAZY_pct),decreasing = T),] #%>% .[1:20,]
CAZY_anno$mean=rowMeans(CAZY_pct)

CAZY_cate3_mean=aggregate(t(CAZY_pct),by=list(Plant=metafile$Plant),FUN=mean)

CAZY_cate3_mean_heat=t(CAZY_cate3_mean[,-1]) %>% as.data.frame()
colnames(CAZY_cate3_mean_heat)=CAZY_cate3_mean$Plant
rownames(CAZY_cate3_mean_heat)=CAZY_anno$Subtype

threshold=1e-3
CAZY_cate3_mean_heat=CAZY_cate3_mean_heat[CAZY_anno$mean>threshold,]
CAZY_anno=CAZY_anno[CAZY_anno$mean>threshold,]


CAZY_anno_heat=data.frame(Type=CAZY_anno$Type,row.names =CAZY_anno$Subtype )
CAZY_col=color_scheme("Plan9",names = levels(cate2_mean_long_mean$Type)) %>% .[1:9]
ann_colors=list(Type=CAZY_col)

CAZY_cate3_mean_heat=CAZY_cate3_mean_heat[,plant_order]
CAZY_cate3_mean_heat_re=data.frame()
for(i in levels(cate2_mean_long_mean$Type)[1:7]){
  sub_frame=CAZY_cate3_mean_heat[CAZY_anno_heat$Type==i,]
  temp=apply(sub_frame, 1, scale) %>% t() %>% as.data.frame()
  sub_frame_order=sub_frame[order(temp[,1]),]
  CAZY_cate3_mean_heat_re=rbind(CAZY_cate3_mean_heat_re,sub_frame_order)
}

library(ggtreeExtra)
library(ggnewscale)
library(ggtree)
library(circlize)
#pdf("./results/indicator/indicator.pdf",width = 6,height = 6)
circos_data=apply(CAZY_cate3_mean_heat_re[-135,], 1, scale) %>% t()
col_fun1 = colorRamp2(c(-2.5, 0, 2.5), c("#4070AF", "white", "#D42C24"))
#CAZY_anno_heat1=data.frame(Type=CAZY_anno_heat[-135,])
top20_stat=cate3_stat[cate3_stat$Subtype %in% CAZY_anno$Subtype[-135],]
circos_anno=data.frame(Subtype=rownames(CAZY_cate3_mean_heat_re)) %>% left_join(CAZY[,1:2]) %>% .[-135,]


circos.par(gap.degree=c(1,1,1,1,1,30),start.degree=90)

#circos.initialize(indicator_anno_long1$tag, x=indicator_anno_long1$statvalue) #initialize

circos.heatmap(circos_data,split = circos_anno$Type,col=col_fun1,cluster = FALSE,
               bg.border = "black",
               bg.lwd = 1,
               cell.border = "white",
               cell.lwd = 0.5,
               rownames.side = "outside",
               rownames.cex = 0.5,
               track.height = 0.25)
## 热图标签
circos.track(
  track.index = get.current.track.index(),
  bg.border = NA,
  panel.fun = function(x,y){
    if(CELL_META$sector.numeric.index == length(unique(circos_anno$Type))) {
      cn <- colnames(CAZY_cate3_mean_heat_re[-135,]) %>% rev()
      n <- length(cn)
      cell_height <- (CELL_META$cell.ylim[2] - CELL_META$cell.ylim[1]) / n
      y_coords <- seq(CELL_META$cell.ylim[1] + cell_height / 2, 
                      CELL_META$cell.ylim[2] - cell_height / 2, 
                      length.out = n)
 
      for (i in 1:n) {
        circos.lines(
          c(CELL_META$cell.xlim[2], CELL_META$cell.xlim[2] + convert_x(1, "mm")),
          c(y_coords[i], y_coords[i]),
          col = "black",
          lwd = 1
        )
      } 
   
      circos.text(
        rep(CELL_META$cell.xlim[2], n) + convert_x(1.5, "mm"), # x坐标
        y_coords, # y坐标
        cn,
        cex = 0.8,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)


circos.track(
  ylim = c(0, 1), 
  track.height = 0.065, 
  bg.col = CAZY_col[unique(circos_anno$Type)],
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2] - 0.75,
      CELL_META$sector.index, 
      facing = "bending.inside", 
      cex = 0.5, 
      adj = c(0.5, 0)
    )
  }
)

# plant pt
circos.track(
  ylim = c(0, 1), 
  track.height = 0.05, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    #sector_data <- CAZY_anno_heat1#[CAZY_anno_heat1$Type == CELL_META$sector.index, ]
    sector_data <- top20_stat[top20_stat$Item=="Plant      ",]
    sector_data$ptcol=case_when(
      sector_data$signif=="***"~"#4070AF",
      sector_data$signif=="**"~"#7293B7",
      sector_data$signif=="*"~"#A5B7BF",
      sector_data$signif=="ns"~"white"
    )
    sector_data <- sector_data[sector_data$Type == CELL_META$sector.index, ]
    for (i in 1:nrow(sector_data))  {
      circos.points(
        CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (i - 0.5) / nrow(sector_data), 
        0.5, 
        pch = 21,
        cex = 0.5,
        bg =sector_data$ptcol[i] ,
        col=NA
      )
    }
    
    if(CELL_META$sector.numeric.index == length(unique(circos_anno$Type))) {
      circos.text(
        CELL_META$cell.xlim[2] + convert_x(1.5, "mm"), 
        0.5,
        "Plant",
        cex = 0.5,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)
# sex pt
circos.track(
  ylim = c(0, 1), 
  track.height = 0.05, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    #sector_data <- CAZY_anno_heat1#[CAZY_anno_heat1$Type == CELL_META$sector.index, ]
    sector_data <- top20_stat[top20_stat$Item=="sex        ",]
    sector_data$ptcol=case_when(
      sector_data$signif=="***"~"#4070AF",
      sector_data$signif=="**"~"#7293B7",
      sector_data$signif=="*"~"#A5B7BF",
      sector_data$signif=="ns"~"white"
    )
    sector_data <- sector_data[sector_data$Type == CELL_META$sector.index, ]
    for (i in 1:nrow(sector_data))  {
      circos.points(
        CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (i - 0.5) / nrow(sector_data), 
        0.5, 
        pch = 21,
        cex = 0.5,
        bg =sector_data$ptcol[i] ,
        col=NA
      )
    }
    
    if(CELL_META$sector.numeric.index == length(unique(circos_anno$Type))) {
      circos.text(
        CELL_META$cell.xlim[2] + convert_x(1.5, "mm"), 
        0.5,
        "Sex",
        cex = 0.5,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)
# interaction pt
circos.track(
  ylim = c(0, 1), 
  track.height = 0.05, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    #sector_data <- CAZY_anno_heat1#[CAZY_anno_heat1$Type == CELL_META$sector.index, ]
    sector_data <- top20_stat[top20_stat$Item=="Plant:sex  ",]
    sector_data$ptcol=case_when(
      sector_data$signif=="***"~"#4070AF",
      sector_data$signif=="**"~"#7293B7",
      sector_data$signif=="*"~"#A5B7BF",
      sector_data$signif=="ns"~"white"
    )
    sector_data <- sector_data[sector_data$Type == CELL_META$sector.index, ]
    for (i in 1:nrow(sector_data))  {
      circos.points(
        CELL_META$xlim[1] + (CELL_META$xlim[2] - CELL_META$xlim[1]) * (i - 0.5) / nrow(sector_data), 
        0.5, 
        pch = 21,
        cex = 0.5,
        bg =sector_data$ptcol[i] ,
        col=NA
      )
    }
    
    if(CELL_META$sector.numeric.index == length(unique(circos_anno$Type))) {
      circos.text(
        CELL_META$cell.xlim[2] + convert_x(1.5, "mm"), 
        0.5,
        "Interaction",
        cex = 0.5,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)
circos.clear()
#save 6x6


##cazy microbe####
load("./tax_obj_plan1.rda")
library(LorMe)
library(dplyr)
CAZY_pct1=CAZY_pct
rownames(CAZY_pct1)=CAZY_cato$Subtype
cazy_genus_network=network_analysis2(input =tax_obj_plan1$Genus_percent,inputtype = 2,n = 24,
                                     threshold = 0.9,method = "spearman",display = T,input2 = CAZY_pct1,input2type = 3)

adj_data=cazy_genus_network$Adjacency_column_table

colnames(adj_data)=c("Genus","Subtype")
adj_data_anno=left_join(adj_data,tax_obj_plan1$Genus_taxonomy) %>% left_join(CAZY_cato)

main_phylum_col=community_plot(taxobj = tax_obj_plan_plant,taxlevel = "Phylum",n = 10,palette = "Paired",nrow = 1) %$%
  filled_color
main_phylum_name= main_phylum_col%>% names()
adj_data_anno$Phylum[!adj_data_anno$Phylum %in%main_phylum_name ]="Others"

subtype_phylum_adj=table(adj_data_anno[,c(6,2)]) %>% as.data.frame()
subtype_phylum_adj=subtype_phylum_adj[subtype_phylum_adj$Freq>2,]

###visualization####
library(dplyr)
library(ggplot2)
library(ggsankeyfier)
library(RColorBrewer)

sankey_edges <- subtype_phylum_adj %>%
  rename(Phylum = Phylum, Subtype = Subtype, Freq = Freq) %>%
  group_by(Phylum, Subtype) %>%
  summarise(Freq = sum(Freq), .groups = "drop")

sankey_long <- pivot_stages_longer(
  data = sankey_edges,
  stages_from = c("Phylum", "Subtype"),
  values_from = "Freq",
  additional_aes_from = "Phylum"   
)

phyla <- unique(sankey_edges$Phylum)
nphy <- length(phyla)

node_cols=c(main_phylum_col,cazy_col)

pos <- position_sankey(order = "ascending", v_space = 10)
pos_text <- position_sankey(order = "ascending", v_space = 10, nudge_x = 0.08)

sankey_long=left_join(sankey_long,data.frame(node=CAZY_cato$Subtype,CAZY_cato))
sankey_long$Type[sankey_long$connector=="from"]=as.character(sankey_long$Phylum[sankey_long$connector=="from"])

sankey_long$Phylum=factor(sankey_long$Phylum,levels = names(main_phylum_col))

ggplot(sankey_long,
       aes(x = stage,
           y = Freq,
           group = node,
           connector = connector,
           edge_id = edge_id)) +
  geom_sankeyedge(aes(fill = Phylum),
                  position = pos,
                  alpha = 0.85,
                  show.legend = TRUE) +
  geom_sankeynode(aes(fill = Type),
                  position = pos,
                  show.legend = FALSE) +
  geom_text(aes(label = node),
            stat = "sankeynode",
            position = pos_text,
            hjust = 0, size = 3) +
  scale_fill_manual(values = node_cols) +
  scale_x_discrete(expand = expansion(add = c(0.2, 0.5)), position = "top") +
  theme_void() +
  theme(legend.position = "right",
        plot.margin = margin(8, 8, 8, 8))
ggsave("./cazy_microbe_alluvial.pdf",width =8,height = 5 )
