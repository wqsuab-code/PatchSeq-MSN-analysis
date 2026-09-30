#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(mclust)})
base<-"C:/Users/53461/OneDrive/Documentos/Patch-seq_Mouse_Acb_MSN_T-type_Visualization/macaque_m";out<-file.path(base,"m18_tempfreeze_NPC5_HCK4_res2.3");dir.create(out,showWarnings=FALSE)
s<-fread(file.path(base,"m18_adaptive_pca126","04_pca_scores.csv"));meta<-fread(file.path(base,"m18_adaptive_pca126","06_PC1_PC5_outlier_diagnostics.csv"))[,.(cell_label,donor_label,ROI,Subclass)]
g<-fread(file.path(base,"m18_adaptive126_all_confusions","03_all_GC_assignments.csv"))[npcs==5 & abs(resolution-2.3)<1e-8]
x<-as.matrix(s[,paste0("PC",1:5),with=FALSE]);hc<-cutree(hclust(dist(x),method="ward.D2"),4);gc<-g$GC[match(s$cell_label,g$cell_label)];stopifnot(!anyNA(gc),uniqueN(gc)==13)
rawtab<-table(HC=factor(hc,1:4),GC=factor(gc,1:13));map<-apply(rawtab,2,which.max);mgc<-unname(map[gc]);tab<-table(HC=factor(hc,1:4),Merged_GC=factor(mgc,1:4));agreement<-mean(hc==mgc);ari<-adjustedRandIndex(hc,mgc)
a<-merge(data.table(cell_label=s$cell_label,HC_K4=hc,GC_raw_K13=gc,GC_merged_K4=mgc,concordant=hc==mgc),meta,by="cell_label");fwrite(a,file.path(out,"01_temp_frozen_assignments_126.csv"),bom=TRUE)
fwrite(data.table(GC_raw=1:13,GC_merged=map),file.path(out,"02_GC13_to_GC4_merge_map.csv"),bom=TRUE)
t<-as.data.table(tab);setnames(t,c("HC","Merged_GC","N"));t[,row_percent:=100*N/sum(N),by=HC];fwrite(t,file.path(out,"03_HC4_GC4_confusion.csv"),bom=TRUE)
p<-ggplot(t,aes(Merged_GC,HC,fill=N))+geom_tile(color="white",linewidth=1)+geom_text(aes(label=sprintf("%d\n%.1f%%",N,row_percent)),size=5)+scale_fill_gradient(low="white",high="#2166AC")+coord_equal()+labs(title="Temporary freeze: nPCS=5, HC K=4, GC res=2.3 (13 -> 4)",subtitle=sprintf("126 cells | agreement %.2f%% | ARI %.3f",100*agreement,ari),x="Merged GC class",y="Fixed HC class",fill="Cells")+theme_bw(base_size=12)
ggsave(file.path(out,"04_temp_frozen_crossplot.png"),p,width=8,height=6,dpi=400,bg="white");ggsave(file.path(out,"04_temp_frozen_crossplot.pdf"),p,width=8,height=6,bg="white")
cat(sprintf("agreement=%.8f ARI=%.8f\n",agreement,ari));print(map);print(tab)
