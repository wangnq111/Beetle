library(patchwork)
pheonotype=read.table("./input/pheono.txt",header = T,sep="\t")
source("./autoplot.R")

metafile=read.table("./input/metafile.txt",header = T,sep="\t")

plant_col=color_scheme("Plan7",names=unique(metafile$Plant))
plant_order=c("HS","YJ","XS","SJ","LS","PG","YS","SHT")
sex_col=c("F"="#FE5D5D",FM="#6376A0")
pheonotype$Plant=factor(pheonotype$Plant,levels = plant_order)



make_violin=function(inputdata,value_col){
  signif_results<- auto_signif_test(data =inputdata,treatment_col =1,value_col =value_col,prior = T,comparison_method = "SNK")
  letters<- data.frame(signif_results$comparison_letters,letterp=max(1.4*(signif_results$comparison_letters %>% .[,'Mean'])+1.4*(signif_results$comparison_letters %>% .[,'std'])))
  ggplot(inputdata,aes(x=as.factor(inputdata[,1]),y=inputdata[,value_col]))+
    scale_y_continuous(expand = c(0.05,0.01))+
    scale_color_manual(values=plant_col)+
    scale_fill_manual(values=plant_col)+
    theme_zg()+
    geom_violin(aes(color=factor(inputdata[,1])),alpha=0.8,width=0.8,linewidth=0.5,trim = F)+
    labs(x='',color='',y=colnames(inputdata)[value_col])+
    geom_quasirandom(aes(color=as.factor(inputdata[,1])),alpha=0.8,size=1,pch=16,show.legend = F)+
    geom_text(data=letters,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
    theme(legend.position = 'right')+guides(fill="none",color="none") %>% return()
}

(make_violin(pheonotype,3)|make_violin(pheonotype,4))
ggsave("./Figure/growth.pdf",width = 8,height = 3)

aov(Length~Sex*Plant,data=pheonotype) %>% summary()
aov(Weight~Sex*Plant,data=pheonotype) %>% summary()
len=auto_signif_test(data =pheonotype,treatment_col =1,value_col =3,prior = T,comparison_method = "SNK")
####Labels####
#
#compare Letters  type     Mean       std   n         se  Min  Max    Q25  Q50    Q75
#SHT     SHT       a Plant 12.86552 0.8038398  58 0.11121228 11.0 14.6 12.325 12.7 13.475
#PG       PG       a Plant 12.85263 0.7497293  95 0.08689703 11.5 15.4 12.500 12.8 13.350
#SJ       SJ       a Plant 12.80889 0.7255598  90 0.08927821 11.3 14.5 12.300 12.7 13.300
#LS       LS       a Plant 12.73699 0.8367442  73 0.09913004 11.0 14.6 12.000 12.7 13.500
#HS       HS       a Plant 12.71029 0.8969801 204 0.05929958 10.0 15.0 12.000 12.7 13.300
#XS       XS       a Plant 12.64297 0.8747268 128 0.07486205 10.0 15.0 12.000 12.5 13.125
#YJ       YJ      ab Plant 12.58627 0.9589619  51 0.11859917 10.9 14.6 12.000 12.4 13.250
#YS       YS       b Plant 12.31556 0.8862473  45 0.12625846  9.7 14.0 12.000 12.5 13.000
wei=auto_signif_test(data =pheonotype,treatment_col =1,value_col =4,prior = T,comparison_method = "SNK")
####Labels####
#
#compare Letters  type     Mean      std   n       se  Min   Max     Q25    Q50     Q75
#SHT     SHT       a Plant 156.5000 30.29375  58 4.331993 76.0 228.4 136.400 151.05 177.625
#LS       LS      ab Plant 152.1685 36.15135  73 3.861360 75.4 237.0 123.200 145.50 176.500
#HS       HS     abc Plant 147.8696 36.10817 204 2.309865 67.0 247.6 120.475 148.45 171.725
#SJ       SJ     abc Plant 146.9889 30.97030  90 3.477607 95.0 210.8 125.050 143.00 171.525
#XS       XS     abc Plant 146.6227 32.44399 128 2.916062 67.0 239.1 124.150 143.65 168.350
#PG       PG      bc Plant 137.6547 28.86415  95 3.384854 89.3 234.0 120.500 130.90 147.050
#YJ       YJ       c Plant 135.5588 33.88652  51 4.619731 77.7 208.1 106.250 134.50 159.150
#YS       YS       d Plant 124.3911 28.29793  45 4.918079 58.3 179.8 105.200 124.60 137.700