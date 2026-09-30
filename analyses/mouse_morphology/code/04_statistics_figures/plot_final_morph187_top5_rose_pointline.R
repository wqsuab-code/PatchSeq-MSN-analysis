#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
options(stringsAsFactors=FALSE)

root <- normalizePath(".",winslash="/",mustWork=TRUE)
base <- file.path(root,"outputs","morph_qc","final_morph187_NPC3_HCK4_GCres0p50_figures")
radar_dir <- file.path(base,"20_top5_feature_radar")
out <- file.path(base,"21_top5_feature_rose_pointline")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

means <- read.csv(file.path(radar_dir,"20_M1-M4_Top5_RobustScaled_Group_Means.csv"),check.names=FALSE)
summary <- read.csv(file.path(radar_dir,"20_M1-M4_Top5_RawUnit_Group_Summary.csv"),check.names=FALSE)
levels_m <- paste0("M",1:4)
cols <- c(M1="#00468B",M2="#42B540",M3="#ED0000",M4="#0099B4")
labels <- c("Soma AR","Soma circ.","Primary N","Bifurcations","Total length")
means$M_class <- factor(means$M_class,levels=levels_m)
means$Metric <- factor(means$Metric,levels=labels)

rings <- do.call(rbind,lapply(c(.25,.5,.75,1),function(v)
  data.frame(Axis=seq_along(labels),Value=v,Ring=factor(v))))
spokes <- data.frame(Axis=seq_along(labels),Value=0,
                     Axis_end=seq_along(labels),Value_end=1)

make_rose <- function(classes=levels_m,facet=TRUE,clean=FALSE){
  d <- means[as.character(means$M_class)%in%classes,,drop=FALSE]
  p <- ggplot()+
    geom_polygon(data=rings,aes(Axis,Value,group=Ring),fill=NA,colour="grey76",linewidth=.25/ggplot2::.pt)+
    geom_segment(data=spokes,aes(Axis,Value,xend=Axis_end,yend=Value_end),colour="grey80",linewidth=.25/ggplot2::.pt)+
    geom_col(data=d,aes(Axis,Mean_scaled_value,fill=M_class),width=.72,alpha=.50,colour=NA)+
    geom_polygon(data=d,aes(Axis,Mean_scaled_value,group=M_class),fill=NA,colour="black",
                 linewidth=.25/ggplot2::.pt,linetype="solid")+
    geom_point(data=d,aes(Axis,Mean_scaled_value),shape=16,size=1,colour="black")+
    scale_fill_manual(values=cols,guide="none")+
    scale_x_continuous(breaks=seq_along(labels),labels=if(clean) NULL else labels,
                       limits=c(.5,length(labels)+.5))+
    scale_y_continuous(limits=c(0,1),breaks=c(.25,.5,.75,1),labels=NULL,expand=c(0,0))+
    coord_polar(start=-pi/2,clip="off")+
    theme_void(base_family="Arial")+
    theme(axis.text.x=element_text(size=if(facet)4.2 else 5,colour="black",
                                  margin=margin(1.5,1.5,1.5,1.5,"pt")),
          strip.text=element_text(face="bold",size=6,colour="black"),
          panel.spacing=grid::unit(8,"pt"),plot.margin=margin(10,10,10,10,"pt"))
  if(facet){
    counts <- tapply(summary$N_cells,summary$M_class,function(x)unique(x)[1])
    facet_labels <- setNames(paste0(levels_m," (n=",as.integer(counts[levels_m]),")"),levels_m)
    p <- p+facet_wrap(~M_class,nrow=1,labeller=as_labeller(facet_labels))
  }
  p
}

save3 <- function(stem,p,w,h){
  ggsave(file.path(out,paste0(stem,".png")),p,width=w,height=h,units="in",dpi=900,
         device=if(requireNamespace("ragg",quietly=TRUE))ragg::agg_png else "png",bg="white")
  ggsave(file.path(out,paste0(stem,".pdf")),p,width=w,height=h,units="in",device=cairo_pdf,bg="white")
  ggsave(file.path(out,paste0(stem,".svg")),p,width=w,height=h,units="in",device=svg,bg="white")
}

combined <- make_rose()+
  labs(caption="Bars: class mean, 50% opacity   |   Black points: size 1   |   Black line: 0.25 pt")+
  theme(plot.caption=element_text(family="Arial",size=5,colour="black",hjust=.5,margin=margin(t=2,unit="pt")))
save3("21_M1-M4_Top5_RoseBars50pct_BlackPoint1_Line0p25pt_W6_H1p8",combined,6,1.8)

clean <- make_rose(clean=TRUE)+theme(axis.text.x=element_blank(),strip.text=element_blank(),
                                     plot.margin=margin(1,1,1,1,"pt"))
save3("21_M1-M4_Top5_RoseBars50pct_BlackPoint1_Line0p25pt_NoLabels_W5p2_H1p35",clean,5.2,1.35)

for(m in levels_m){
  p <- make_rose(m,facet=FALSE)+labs(title=paste0(m," morphology"))+
    theme(plot.title=element_text(face="bold",size=6,hjust=.5))
  save3(paste0("21_",m,"_Top5_RoseBars50pct_BlackPoint1_Line0p25pt_W1p6_H1p6"),p,1.6,1.6)
}

write.csv(means,file.path(out,"21_M1-M4_Top5_Rose_Input_GroupMeans.csv"),row.names=FALSE)
writeLines(c("Final Morph n=187 top-five rose plus point-line plots",
             "Bars use the M-class colour at alpha=0.50.",
             "Mean points are solid black, ggplot2 size=1.",
             "Connecting line is solid black, 0.25 pt.",
             "Values are the same cohort-wide Q2.5/Q97.5 winsorized 0-1 scaling used in the E radar analysis."),
           file.path(out,"21_M1-M4_Top5_Rose_run_log.txt"))
cat(normalizePath(out,winslash="/"),"\n")
