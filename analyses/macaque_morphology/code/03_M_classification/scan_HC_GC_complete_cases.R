#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Seurat);library(mclust);library(ggplot2);library(jsonlite)})
base_dir<-"C:/Users/53461/OneDrive/Documentos/Patch-seq_Mouse_Acb_MSN_T-type_Visualization/macaque_m"
args<-commandArgs(trailingOnly=TRUE);cohort_n<-if(length(args))as.integer(args[[1]]) else 117L;stopifnot(cohort_n%in%c(117L,126L))
pca_dir<-file.path(base_dir,sprintf("m18_adaptive_pca%d",cohort_n));out<-file.path(base_dir,sprintf("m18_adaptive%d_all_confusions",cohort_n));dir.create(out,recursive=TRUE,showWarnings=FALSE)
npcs_grid<-3:5;k_grid<-3:8;res_grid<-seq(.5,3,by=.1);knn<-20L;seed<-777L
s<-fread(file.path(pca_dir,"04_pca_scores.csv"));stopifnot(nrow(s)==cohort_n)

optimal_bijection<-function(tab){
  tab<-as.matrix(tab);k<-nrow(tab);if(ncol(tab)!=k)return(NULL);states<-2^k
  dp<-rep(-Inf,states);dp[1]<-0;pm<-matrix(NA_integer_,k,states);ph<-matrix(NA_integer_,k,states)
  for(j in seq_len(k)){nd<-rep(-Inf,states);for(mask in which(is.finite(dp))-1L)for(h in seq_len(k))if(bitwAnd(mask,bitwShiftL(1L,h-1L))==0){nm<-bitwOr(mask,bitwShiftL(1L,h-1L));pos<-nm+1;sc<-dp[mask+1]+tab[h,j];if(sc>nd[pos]){nd[pos]<-sc;pm[j,pos]<-mask;ph[j,pos]<-h}};dp<-nd}
  map<-integer(k);mask<-states-1L;for(j in k:1L){pos<-mask+1;map[j]<-ph[j,pos];mask<-pm[j,pos]};map
}

all_sum<-list();all_cells<-list();all_mats<-list();si<-ci<-mi<-0L
for(npcs in npcs_grid){
  x<-as.matrix(s[,paste0("PC",1:npcs),with=FALSE]);rownames(x)<-s$cell_label
  hc_all<-setNames(lapply(k_grid,function(k)cutree(hclust(dist(x),method="ward.D2"),k)),k_grid)
  snn<-FindNeighbors(x,k.param=knn,compute.SNN=TRUE,prune.SNN=1/15,verbose=FALSE)[["snn"]]
  for(res in res_grid){
    set.seed(seed);fit<-FindClusters(snn,algorithm=1,resolution=res,random.seed=seed,n.start=30,n.iter=30,verbose=FALSE)
    gc<-as.integer(as.factor(fit[,ncol(fit)]));rawK<-uniqueN(gc)
    ci<-ci+1L;all_cells[[ci]]<-data.table(cell_label=s$cell_label,npcs=npcs,resolution=res,GC=gc,GC_raw_K=rawK)
    for(k in k_grid){
      hc<-hc_all[[as.character(k)]];tab<-table(HC=factor(hc,levels=seq_len(k)),GC=factor(gc,levels=seq_len(rawK)))
      td<-as.data.table(tab);setnames(td,c("HC","GC","N"));td[,`:=`(npcs=npcs,resolution=res,HC_K=k,GC_raw_K=rawK)]
      setcolorder(td,c("npcs","resolution","HC_K","GC_raw_K","HC","GC","N"));mi<-mi+1L;all_mats[[mi]]<-td
      acc<-minrec<-minprec<-NA_real_;mapping<-""
      if(rawK==k){mp<-optimal_bijection(tab);mapped<-mp[gc];ct<-table(factor(hc,levels=seq_len(k)),factor(mapped,levels=seq_len(k)));acc<-mean(hc==mapped);minrec<-min(diag(ct)/rowSums(ct));minprec<-min(diag(ct)/colSums(ct));mapping<-paste(sprintf("GC%d->HC%d",seq_len(k),mp),collapse=";")}
      si<-si+1L;all_sum[[si]]<-data.table(npcs=npcs,resolution=res,HC_K=k,GC_raw_K=rawK,HC_sizes=paste(as.integer(table(hc)),collapse=";"),GC_sizes=paste(as.integer(table(gc)),collapse=";"),ARI=adjustedRandIndex(hc,gc),equal_K=rawK==k,optimal_agreement=acc,min_HC_recall=minrec,min_HC_precision=minprec,optimal_mapping=mapping)
    }
  }
}
summary<-rbindlist(all_sum);mats<-rbindlist(all_mats);cells<-rbindlist(all_cells)
fwrite(summary,file.path(out,"01_all_HC_GC_comparison_summary.csv"),bom=TRUE)
fwrite(mats,file.path(out,"02_all_cross_confusion_matrices_long.csv"),bom=TRUE)
fwrite(cells,file.path(out,"03_all_GC_assignments.csv"),bom=TRUE)
fwrite(summary[equal_K==TRUE][order(npcs,HC_K,-optimal_agreement)],file.path(out,"04_equalK_confusion_summary.csv"),bom=TRUE)

pdf(file.path(out,"05_all_cross_confusion_heatmaps.pdf"),width=8,height=6,onefile=TRUE)
for(i in seq_len(nrow(summary))){q<-summary[i];d<-mats[npcs==q$npcs & abs(resolution-q$resolution)<1e-8 & HC_K==q$HC_K]
  p<-ggplot(d,aes(GC,HC,fill=N))+geom_tile(color="white")+geom_text(aes(label=N),size=3)+scale_fill_gradient(low="white",high="#2166AC")+scale_x_discrete(drop=FALSE)+scale_y_discrete(drop=FALSE)+coord_equal()+labs(title=sprintf("NPC%d | HC K=%d vs GC res=%.1f (raw K=%d)",q$npcs,q$HC_K,q$resolution,q$GC_raw_K),subtitle=sprintf("ARI=%.3f%s",q$ARI,if(q$equal_K)sprintf(" | optimal agreement=%.1f%%",100*q$optimal_agreement)else ""),x="GC raw cluster",y="HC Ward.D2 cluster")+theme_bw(base_size=10)
  print(p)
}
dev.off()

write_json(list(analysis=sprintf("Macaque adaptive-transform %d-cell all HC-GC cross-confusions",cohort_n),input_PCA=file.path(pca_dir,"04_pca_scores.csv"),n=cohort_n,npcs=npcs_grid,HC_K=k_grid,HC="one Euclidean Ward.D2 tree per NPC, cut K=3:8",GC=list(kNN=knn,resolution=range(res_grid),step=.1,seed=seed,algorithm=1),HC_guided_GC_merge=FALSE,total_matrices=nrow(summary)),file.path(out,"00_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
print(summary[equal_K==TRUE][order(-optimal_agreement),.(npcs,resolution,HC_K,GC_raw_K,ARI,optimal_agreement,min_HC_recall,min_HC_precision,HC_sizes,GC_sizes)])
