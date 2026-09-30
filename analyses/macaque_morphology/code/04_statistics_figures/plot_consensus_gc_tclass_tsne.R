#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel)})
base <- "C:/Users/53461/OneDrive/Documentos/Patch-seq_Mouse_Acb_MSN_T-type_Visualization/macaque_m/m18_tempfreeze_NPC5_HCK4_res2.3"
out <- file.path(base,"tSNE_consensus_GC_Tclass")
dir.create(out,showWarnings=FALSE)
d <- fread(file.path(base,"09_tSNE_optimized_coordinates.csv"))
a <- fread(file.path(base,"01_temp_frozen_assignments_126.csv"))[,.(cell_label,Subclass)]
d <- merge(d,a,by="cell_label",all.x=TRUE)
cols <- c("1"="#1F77B4","2"="#D9A400","3"="#8C564B","4"="#E377C2")
d[,`:=`(cx=mean(tSNE1),cy=mean(tSNE2)),by=HC_K4]
d[,`:=`(x=cx+.85*(tSNE1-cx),y=cy+.85*(tSNE2-cy))]
con <- d[concordant==TRUE]

style <- theme_void(base_size=12)+theme(plot.title=element_text(face="bold"),legend.position="right",plot.margin=margin(10,14,10,10))
save_plot <- function(p,name) {
  ggsave(file.path(out,paste0(name,".png")),p,width=5.4,height=4.7,dpi=600,bg="white")
  ggsave(file.path(out,paste0(name,".pdf")),p,width=5.4,height=4.7,bg="white")
}
no_text <- function(p) {
  p$layers <- Filter(function(x) !inherits(x$geom,"GeomLabelRepel"),p$layers)
  p+labs(title=NULL,subtitle=NULL,color=NULL)+theme(legend.position="none")
}

cent1 <- con[,.(x=median(x),y=median(y),N=.N),by=HC_K4]
p1 <- ggplot(d,aes(x,y))+
  stat_ellipse(data=con,aes(color=factor(HC_K4),group=HC_K4),level=.80,type="norm",linewidth=.7,linetype="dashed",show.legend=FALSE)+
  geom_point(data=d[concordant==FALSE],color="#4D4D4D",size=2.1,alpha=.90)+
  geom_point(data=con,aes(color=factor(HC_K4)),size=2.1,alpha=.88)+
  geom_label_repel(data=cent1,aes(label=sprintf("M%d (n=%d)",HC_K4,N),color=factor(HC_K4)),fill="white",label.size=.2,fontface="bold",size=3.6,show.legend=FALSE)+
  scale_color_manual(values=cols)+coord_equal()+labs(title="HC-GC consensus morphology classes",subtitle="117 consensus cells; 9 discordant cells in dark gray",color="Consensus M")+style
save_plot(p1,"01_consensus_M4_darkgray_discordant")
save_plot(no_text(p1),"01_consensus_M4_darkgray_discordant_no_text")

cent2 <- d[,.(x=median(x),y=median(y),N=.N),by=GC_merged_K4]
p2 <- ggplot(d,aes(x,y))+
  stat_ellipse(data=con,aes(color=factor(GC_merged_K4),group=GC_merged_K4),level=.80,type="norm",linewidth=.7,linetype="dashed",show.legend=FALSE)+
  geom_point(aes(color=factor(GC_merged_K4)),size=2.1,alpha=.88)+
  geom_label_repel(data=cent2,aes(label=sprintf("GC%d (n=%d)",GC_merged_K4,N),color=factor(GC_merged_K4)),fill="white",label.size=.2,fontface="bold",size=3.6,show.legend=FALSE)+
  scale_color_manual(values=cols)+coord_equal()+labs(title="Merged graph-clustering classes",subtitle="All 126 cells; 80% boundaries fitted from 117 consensus cells",color="Merged GC")+style
save_plot(p2,"02_merged_GC4_all_cells")
save_plot(no_text(p2),"02_merged_GC4_all_cells_no_text")

make_tclass <- function(subclass,label,file) {
  target <- con[Subclass==subclass]
  p <- ggplot(d,aes(x,y))+
    stat_ellipse(data=con,aes(color=factor(HC_K4),group=HC_K4),level=.80,type="norm",linewidth=.7,linetype="dashed",show.legend=FALSE)+
    geom_point(data=con[Subclass!=subclass],color="#D9D9D9",size=2.0,alpha=.70)+
    geom_point(data=d[concordant==FALSE],color="#4D4D4D",size=2.1,alpha=.95)+
    geom_point(data=target,aes(color=factor(HC_K4)),size=2.25,alpha=.92)+
    scale_color_manual(values=cols)+coord_equal()+
    labs(title=paste0(label," distribution"),subtitle=sprintf("%d of 117 consensus cells; discordant cells in dark gray",nrow(target)),color="Consensus M")+style
  save_plot(p,file)
  save_plot(no_text(p),paste0(file,"_no_text"))
}
make_tclass("STR D1 MSN","D1","03_D1_distribution")
make_tclass("STR D2 MSN","D2","04_D2_distribution")
make_tclass("STR Hybrid MSN","Hybrid","05_Hybrid_distribution")

stopifnot(nrow(d)==126,nrow(con)==117,identical(as.integer(table(con$HC_K4)),c(43L,42L,20L,12L)))
cat(normalizePath(out,winslash="/"),"\n")
