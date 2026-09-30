#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
base <- "C:/Users/53461/OneDrive/Documentos/Patch-seq_Mouse_Acb_MSN_T-type_Visualization/macaque_m/m18_tempfreeze_NPC5_HCK4_res2.3"
out <- file.path(base,"panels_A_I")
d <- fread(file.path(base,"01_temp_frozen_assignments_126.csv"))[concordant==TRUE]
d[,M:=factor(paste0("M",HC_K4),levels=paste0("M",1:4))]
d[,T_class:=factor(fcase(Subclass=="STR D1 MSN","D1",Subclass=="STR D2 MSN","D2",default="Hybrid"),levels=c("D1","D2","Hybrid"))]
counts <- d[,.(N=.N),by=.(M,T_class)]
counts[,Percent:=100*N/sum(N),by=M]
fwrite(counts,file.path(out,"K_consensusM4_Tclass_counts.csv"),bom=TRUE)
cols <- c(D1="#D95F02",D2="#008F7A",Hybrid="#6E7180")

p <- ggplot(counts,aes(M,Percent,fill=T_class))+
  geom_col(width=.68,position=position_stack(reverse=TRUE))+
  geom_text(aes(label=sprintf("%.1f%%",Percent)),position=position_stack(vjust=.5,reverse=TRUE),color="white",fontface="bold",size=2.4)+
  scale_fill_manual(values=cols)+scale_y_continuous(breaks=c(0,25,50,75,100),labels=function(x)paste0(x,"%"),expand=expansion(mult=c(0,.02)))+
  labs(x=NULL,y="Cell proportion",fill="T class")+
  theme_classic(base_size=12)+theme(axis.line.x=element_blank(),axis.ticks.x=element_blank(),legend.position="right")

p_blank <- ggplot(counts,aes(M,Percent,fill=T_class))+
  geom_col(width=.68,position=position_stack(reverse=TRUE))+
  scale_fill_manual(values=cols)+scale_y_continuous(expand=expansion(mult=c(0,.02)))+
  theme_void()+theme(legend.position="none",plot.margin=margin(4,4,4,4))

for (ext in c("png","pdf")) {
  ggsave(file.path(out,paste0("K_consensusM4_Tclass_stacked_percent_labeled.",ext)),p,width=3.8,height=3.5,dpi=600,bg="white")
  ggsave(file.path(out,paste0("K_consensusM4_Tclass_stacked_percent_no_labels.",ext)),p_blank,width=2.7,height=3.5,dpi=600,bg="white")
}
stopifnot(sum(counts$N)==117,all(abs(counts[,sum(Percent),by=M]$V1-100)<1e-10),identical(as.integer(d[,table(M)]),c(43L,42L,20L,12L)))
cat(file.path(out,"K_consensusM4_Tclass_stacked_percent_labeled.png"),"\n")
