#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite);library(digest);library(car)})
zip_path<-"C:/Users/53461/Downloads/Macaque-PatchSeq-BG.zip"
base_dir<-"C:/Users/53461/OneDrive/Documentos/Patch-seq_Mouse_Acb_MSN_T-type_Visualization/macaque_m"
args<-commandArgs(trailingOnly=TRUE);cohort_mode<-if(length(args))args[[1]] else "fixed117";stopifnot(cohort_mode%in%c("all126","fixed117"))
out<-file.path(base_dir,if(cohort_mode=="all126")"m18_adaptive_pca126" else "m18_adaptive_pca117");dir.create(out,recursive=TRUE,showWarnings=FALSE)
sha<-"8b0aeaed726e27658066230fe467bdbb765d641d1ea1e3a115e079e4ff0c13a6"
stopifnot(tolower(digest(file=zip_path,algo="sha256",serialize=FALSE))==sha)
features<-c("basal_dendrite_bias_dorsal","basal_dendrite_bias_medial","basal_dendrite_calculate_number_of_stems","basal_dendrite_extent_dorsal","basal_dendrite_extent_medial","basal_dendrite_max_branch_order","basal_dendrite_max_euclidean_distance","basal_dendrite_max_path_distance","basal_dendrite_mean_contraction","basal_dendrite_mean_diameter","basal_dendrite_num_branches","basal_dendrite_soma_percentile_dorsal","basal_dendrite_soma_percentile_medial","basal_dendrite_stem_exit_MedialLateral","basal_dendrite_stem_exit_dorsal","basal_dendrite_stem_exit_ventral","basal_dendrite_total_length","soma_surface_area")
meta<-as.data.table(read.csv(unz(zip_path,"Data/cell_metadata_AllCell.csv"),check.names=FALSE));morph<-as.data.table(read.csv(unz(zip_path,"Data/morphology_features.csv"),check.names=FALSE))
d<-merge(meta,morph,by="cell_label",sort=FALSE)[Lib_region_of_interest_label%chin%c("Ca","Pu","NAC") & Subclass_name%chin%c("STR D1 MSN","STR D2 MSN","STR Hybrid MSN")]
d<-d[complete.cases(d[,..features])];setorder(d,cell_label);stopifnot(nrow(d)==126L)
if(cohort_mode=="fixed117"){
  old_outliers<-fread(file.path(base_dir,"m18_kvote_fixed117","00_fixed_PC1_5_3SD_isolated_cells.csv"))$cell_label
  d<-d[!cell_label%chin%old_outliers]
}
stopifnot(nrow(d)==if(cohort_mode=="all126")126L else 117L)

adj_skew<-function(x){n<-length(x);m<-mean(x);s<-sd(x);n/((n-1)*(n-2))*sum(((x-m)/s)^3)}
yj_fit<-function(x){
  objective<-function(lambda){y<-car::yjPower(x,lambda); n<-length(y); -( -n/2*log(mean((y-mean(y))^2)) + (lambda-1)*sum(sign(x)*log(abs(x)+1)) )}
  optimize(objective,c(-5,5))$minimum
}
raw<-as.matrix(d[,..features]);rownames(raw)<-d$cell_label
raw_skew<-apply(raw,2,adj_skew);yj_features<-names(raw_skew)[abs(raw_skew)>=.5];direct_features<-setdiff(features,yj_features)
trans<-raw;lambdas<-setNames(rep(NA_real_,length(features)),features)
for(f in yj_features){lambdas[f]<-yj_fit(raw[,f]);trans[,f]<-car::yjPower(raw[,f],lambdas[f])}
z<-scale(trans);post_skew<-apply(trans,2,adj_skew)
pca<-prcomp(z,center=FALSE,scale.=FALSE);scores<-pca$x;pcz<-scale(scores[,1:5,drop=FALSE]);maxz<-apply(abs(pcz),1,max);flag<-maxz>3

audit<-data.table(feature=features,raw_skewness=raw_skew[features],rule=ifelse(features%chin%yj_features,"Yeo-Johnson + Z-score","direct Z-score"),YJ_lambda=lambdas[features],post_transform_skewness=post_skew[features],variance_after_Z=apply(z,2,var))
fwrite(audit,file.path(out,"01_transform_skewness_audit.csv"),bom=TRUE)
fwrite(data.table(cell_label=d$cell_label,z),file.path(out,"02_transformed_z_117.csv"),bom=TRUE)
var_dt<-data.table(PC=seq_along(pca$sdev),variance_percent=100*pca$sdev^2/sum(pca$sdev^2));var_dt[,cumulative_percent:=cumsum(variance_percent)]
fwrite(var_dt,file.path(out,"03_pca_variance.csv"),bom=TRUE)
fwrite(data.table(cell_label=d$cell_label,scores),file.path(out,"04_pca_scores.csv"),bom=TRUE)
fwrite(data.table(feature=rownames(pca$rotation),pca$rotation),file.path(out,"05_pca_loadings.csv"),bom=TRUE)
outdt<-data.table(cell_label=d$cell_label,donor_label=d$donor_label,ROI=d$Lib_region_of_interest_label,Subclass=d$Subclass_name,max_abs_PC1_5_z=maxz,PC_exceedances=apply(abs(pcz)>3,1,function(v)paste(colnames(pcz)[v],collapse=";")),flag_PC1_5_3SD=flag)
fwrite(outdt,file.path(out,"06_PC1_PC5_outlier_diagnostics.csv"),bom=TRUE)
fwrite(outdt[flag==TRUE][order(-max_abs_PC1_5_z)],file.path(out,"07_new_PC1_PC5_3SD_candidates.csv"),bom=TRUE)

p1<-ggplot(var_dt[PC<=10],aes(PC,variance_percent))+geom_col(fill="#4C78A8")+geom_line(aes(y=cumulative_percent),color="#E45756")+geom_point(aes(y=cumulative_percent),color="#E45756")+scale_x_continuous(breaks=1:10)+labs(title=sprintf("Adaptive-transform PCA: %d cells",nrow(d)),y="Variance / cumulative percent",x="PC")+theme_classic()
ggsave(file.path(out,"08_scree.png"),p1,width=7,height=4.5,dpi=400,bg="white")
plotdt<-data.table(PC1=scores[,1],PC2=scores[,2],outlier=flag,cell_label=d$cell_label)
p2<-ggplot(plotdt,aes(PC1,PC2,color=outlier))+geom_point(size=2)+scale_color_manual(values=c(`FALSE`="#4C78A8",`TRUE`="#E45756"))+labs(title="PC1-PC2 after adaptive transformation",color="PC1-5 >3 SD")+theme_classic()
ggsave(file.path(out,"09_PC1_PC2_outliers.png"),p2,width=6,height=5,dpi=400,bg="white")
write_json(list(analysis=paste("Macaque",cohort_mode,"adaptive feature transformation and PCA"),source_zip=normalizePath(zip_path,winslash="/"),source_zip_sha256=sha,mouse_morphology_results_read=FALSE,n=nrow(d),features=features,rule="abs raw adjusted skewness >=0.5: maximum-likelihood Yeo-Johnson; otherwise direct Z-score",PCA_outlier_rule="any of PC1-PC5 absolute standardized score >3",clustering_performed=FALSE),file.path(out,"00_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
print(audit);print(var_dt[1:5]);print(outdt[flag==TRUE][order(-max_abs_PC1_5_z)])
