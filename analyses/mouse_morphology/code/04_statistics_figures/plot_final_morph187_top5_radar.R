#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
options(stringsAsFactors = FALSE)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
input_file <- file.path(root, "outputs", "morph_qc", "final_morph187_NPC3_HCK4_GCres0p50_figures", "00_final_morph187_cell_assignments.csv")
output_dir <- file.path(root, "outputs", "morph_qc", "final_morph187_NPC3_HCK4_GCres0p50_figures", "20_top5_feature_radar")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

id <- "MSN_unique_ID"
metric_columns <- c(
  "M_soma_aspect_ratio",
  "M_soma_circularity_index",
  "M_total_number_of_neurites",
  "M_Number_of_bifurcation_points",
  "M_Total_neurite_length_(sections)"
)
metric_labels <- c("Soma AR", "Soma circ.", "Primary N", "Bifurcations", "Total length")
m_levels <- paste0("M", 1:4)
m_colours <- c(M1="#00468B", M2="#42B540", M3="#ED0000", M4="#0099B4")
cell_alpha <- 0.30

d <- read.csv(input_file, check.names=FALSE)
stopifnot(nrow(d)==187L, all(c(id,"M_class",metric_columns) %in% names(d)))
d$M_class <- factor(d$M_class, levels=m_levels)
raw <- as.matrix(d[,metric_columns,drop=FALSE]); storage.mode(raw) <- "double"
stopifnot(!anyNA(raw), all(is.finite(raw)))

# Same rule used for the E-class radar: winsorize each feature at Q2.5/Q97.5,
# then map to 0--1 without reversing biological direction.
lower <- apply(raw,2,quantile,probs=.025,na.rm=TRUE)
upper <- apply(raw,2,quantile,probs=.975,na.rm=TRUE)
scaled <- sweep(raw,2,lower,pmax)
scaled <- sweep(scaled,2,upper,pmin)
scaled <- sweep(scaled,2,lower,"-")
scaled <- sweep(scaled,2,pmax(upper-lower,.Machine$double.eps),"/")
scaled[scaled<0] <- 0; scaled[scaled>1] <- 1
colnames(scaled) <- metric_labels

cell_long <- data.frame(
  MSN_unique_ID=rep(d[[id]],each=length(metric_labels)),
  M_class=rep(d$M_class,each=length(metric_labels)),
  Metric=factor(rep(metric_labels,times=nrow(d)),levels=metric_labels),
  Axis=rep(seq_along(metric_labels),times=nrow(d)),
  Scaled_value=as.vector(t(scaled))
)
mean_rows <- do.call(rbind,lapply(m_levels,function(m){
  keep <- d$M_class==m
  data.frame(M_class=factor(m,levels=m_levels),Metric=factor(metric_labels,levels=metric_labels),
             Axis=seq_along(metric_labels),Mean_scaled_value=colMeans(scaled[keep,,drop=FALSE]))
}))
raw_summary <- do.call(rbind,lapply(m_levels,function(m){
  keep <- d$M_class==m
  data.frame(M_class=m,N_cells=sum(keep),Metric=metric_labels,Original_column=metric_columns,
             Mean_raw=colMeans(raw[keep,,drop=FALSE]),SD_raw=apply(raw[keep,,drop=FALSE],2,sd),
             Median_raw=apply(raw[keep,,drop=FALSE],2,median),
             Mean_scaled=colMeans(scaled[keep,,drop=FALSE]))
}))

rings <- do.call(rbind,lapply(c(.25,.5,.75,1),function(v)
  data.frame(Axis=seq_along(metric_labels),Value=v,Ring=factor(v))))
spokes <- data.frame(Axis=seq_along(metric_labels),Value=0,
                     Axis_end=seq_along(metric_labels),Value_end=1)

make_radar <- function(classes=m_levels,facet=TRUE,clean=FALSE){
  cells <- cell_long[as.character(cell_long$M_class)%in%classes,,drop=FALSE]
  means <- mean_rows[as.character(mean_rows$M_class)%in%classes,,drop=FALSE]
  p <- ggplot()+
    geom_polygon(data=rings,aes(Axis,Value,group=Ring),fill=NA,colour="grey72",linewidth=.25/ggplot2::.pt)+
    geom_segment(data=spokes,aes(Axis,Value,xend=Axis_end,yend=Value_end),colour="grey76",linewidth=.25/ggplot2::.pt)+
    geom_path(data=cells,aes(Axis,Scaled_value,group=interaction(M_class,MSN_unique_ID),colour=M_class),
              alpha=cell_alpha,linewidth=.20/ggplot2::.pt)+
    geom_path(data=means,aes(Axis,Mean_scaled_value,group=M_class),colour="black",linewidth=1.05/ggplot2::.pt)+
    geom_point(data=means,aes(Axis,Mean_scaled_value),colour="black",size=.65)+
    scale_colour_manual(values=m_colours,guide="none")+
    scale_x_continuous(breaks=seq_along(metric_labels),labels=if(clean) NULL else metric_labels,
                       limits=c(.5,length(metric_labels)+.5))+
    scale_y_continuous(limits=c(0,1),breaks=c(.25,.5,.75,1),labels=NULL)+
    coord_polar(start=-pi/2,clip="off")+
    theme_void(base_family="Arial")+
    theme(axis.text.x=element_text(size=if(facet) 4.2 else 5,colour="black",
                                  margin=margin(1.5,1.5,1.5,1.5,"pt")),
          strip.text=element_text(face="bold",size=6,colour="black"),
          panel.spacing=grid::unit(8,"pt"),plot.margin=margin(10,10,10,10,"pt"))
  if(facet){
    counts <- table(d$M_class)
    labs <- setNames(paste0(names(counts)," (n=",as.integer(counts),")"),names(counts))
    p <- p+facet_wrap(~M_class,nrow=1,labeller=as_labeller(labs))
  }
  p
}

save_pair <- function(stem,p,w,h){
  ggsave(file.path(output_dir,paste0(stem,".png")),p,width=w,height=h,units="in",dpi=900,
         device=if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else "png",bg="white")
  ggsave(file.path(output_dir,paste0(stem,".pdf")),p,width=w,height=h,units="in",device=cairo_pdf,bg="white")
  ggsave(file.path(output_dir,paste0(stem,".svg")),p,width=w,height=h,units="in",device=svg,bg="white")
}

combined <- make_radar()+
  labs(caption="Coloured lines: individual cells   |   Black line and points: within-class mean")+
  theme(plot.caption=element_text(family="Arial",size=5,colour="black",hjust=.5,margin=margin(t=2,unit="pt")))
save_pair("20_M1-M4_Top5_CellwiseAlpha0p30_PlusMean_Radar_W6_H1p8",combined,6,1.8)

clean <- make_radar(clean=TRUE)+theme(axis.text.x=element_blank(),strip.text=element_blank(),
                                      plot.margin=margin(1,1,1,1,"pt"))
save_pair("20_M1-M4_Top5_Radar_NoLabels_W5p2_H1p35",clean,5.2,1.35)

for(m in m_levels){
  p <- make_radar(m,facet=FALSE)+labs(title=paste0(m," morphology"))+
    theme(plot.title=element_text(face="bold",size=6,hjust=.5))
  save_pair(paste0("20_",m,"_Top5_CellwiseAlpha0p30_PlusMean_Radar_W1p6_H1p6"),p,1.6,1.6)
}

# Direct comparison of the four class means on one shared radar.
overlay <- ggplot()+
  geom_polygon(data=rings,aes(Axis,Value,group=Ring),fill=NA,colour="grey76",linewidth=.25/ggplot2::.pt)+
  geom_segment(data=spokes,aes(Axis,Value,xend=Axis_end,yend=Value_end),colour="grey80",linewidth=.25/ggplot2::.pt)+
  geom_path(data=mean_rows,aes(Axis,Mean_scaled_value,group=M_class,colour=M_class),linewidth=.8)+
  geom_point(data=mean_rows,aes(Axis,Mean_scaled_value,colour=M_class),size=.8)+
  scale_colour_manual(values=m_colours,name=NULL)+
  scale_x_continuous(breaks=seq_along(metric_labels),labels=metric_labels,limits=c(.5,5.5))+
  scale_y_continuous(limits=c(0,1),breaks=c(.25,.5,.75,1),labels=NULL)+
  coord_polar(start=-pi/2,clip="off")+theme_void(base_family="Arial")+
  theme(axis.text.x=element_text(size=5,colour="black",margin=margin(2,2,2,2,"pt")),
        legend.position="bottom",legend.text=element_text(size=5),legend.key.width=grid::unit(9,"pt"),
        plot.margin=margin(10,10,4,10,"pt"))
save_pair("20_M1-M4_Top5_GroupMean_Overlay_Radar_W2_H2",overlay,2,2)

write.csv(raw_summary,file.path(output_dir,"20_M1-M4_Top5_RawUnit_Group_Summary.csv"),row.names=FALSE)
write.csv(cell_long,file.path(output_dir,"20_M1-M4_Top5_Cellwise_RobustScaled_Values.csv"),row.names=FALSE)
write.csv(mean_rows,file.path(output_dir,"20_M1-M4_Top5_RobustScaled_Group_Means.csv"),row.names=FALSE)
writeLines(c("Final Morph n=187 top-five radar plots",
             "Scaling: winsorized at cohort Q2.5/Q97.5 and mapped to 0-1; no feature direction reversal.",
             "Coloured open paths: individual cells (alpha 0.30); black path and points: class mean.",
             "Top five selected by nested grouped-CV permutation importance.",
             paste0("Class counts: ",paste(names(table(d$M_class)),as.integer(table(d$M_class)),sep="=",collapse="; "))),
           file.path(output_dir,"20_M1-M4_Top5_Radar_run_log.txt"))
cat(normalizePath(output_dir,winslash="/"),"\n")
