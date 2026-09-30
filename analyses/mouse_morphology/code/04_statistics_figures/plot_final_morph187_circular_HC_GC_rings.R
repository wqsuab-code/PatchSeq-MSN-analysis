#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(circlize))

scores_file <- file.path("outputs", "morph_qc", "morph187_NPC3-4_HCK4-5_GC_allres_raw", "01_PCA_scores_187cells.csv")
assign_file <- file.path("outputs", "morph_qc", "final_morph187_NPC3_HCK4_GCres0p50_figures", "00_final_morph187_cell_assignments.csv")
out <- file.path("outputs", "morph_qc", "final_morph187_NPC3_HCK4_GCres0p50_figures")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

scores <- read.csv(scores_file, check.names=FALSE)
lab <- read.csv(assign_file, check.names=FALSE)
lab <- lab[match(scores$MSN_unique_ID, lab$MSN_unique_ID),]
stopifnot(nrow(scores)==187L, !anyNA(lab$MSN_unique_ID))
lab$HC_GC_consensus <- tolower(as.character(lab$HC_GC_consensus)) == "true"
pc <- as.matrix(scores[,paste0("PC",1:3)]); rownames(pc)<-scores$MSN_unique_ID
tree <- hclust(dist(pc), method="ward.D2"); dend <- as.dendrogram(tree)
ids <- labels(dend); lab <- lab[match(ids,lab$MSN_unique_ID),]

cols <- c(M1="#00468B",M2="#42B540",M3="#ED0000",M4="#0099B4")
hc_col <- unname(cols[lab$M_class]); gc_col <- unname(cols[lab$GC_as_M_class])
n <- length(ids); maxh <- attr(dend,"height")

leaf_count <- function(node){
  if(is.leaf(node)) return(1L)
  sum(vapply(node, leaf_count, integer(1)))
}

style_edges <- function(node,parent_height){
  if(is.leaf(node)) return(node)
  for(i in seq_along(node)){
    child<-node[[i]]; ep<-attr(child,"edgePar");if(is.null(ep))ep<-list()
    radial_level <- max(0,min(1,parent_height/maxh))
    subtree_support <- leaf_count(child)/n
    # Strong continuous hierarchy: central/high-support trunks are clearly
    # thicker, whereas distal terminal branches remain fine.
    hierarchy <- .72*radial_level^1.35 + .28*sqrt(subtree_support)
    ep$col <- "#888888"
    ep$lwd <- .25 + 2.35*hierarchy
    attr(child,"edgePar")<-ep;node[[i]]<-style_edges(child,attr(child,"height"))
  };node
}
dend<-style_edges(dend,maxh)

ring <- function(cc,h=.047){
  circos.trackPlotRegion(ylim=c(0,1),track.height=h,bg.border=NA,panel.fun=function(x,y){
    for(i in seq_len(n))circos.rect(i-1,0,i,1,col=cc[i],border=NA)
  })
}
draw_core <- function(){
  circos.clear();circos.par(start.degree=90,gap.degree=0,cell.padding=c(0,0,0,0),track.margin=c(.002,.002),circle.margin=c(.008,.008,.008,.008))
  circos.initialize(factors="cells",xlim=c(0,n));ring(gc_col);ring(hc_col)
  circos.trackPlotRegion(ylim=c(0,maxh),track.height=.80,bg.border=NA,panel.fun=function(x,y)circos.dendrogram(dend,facing="outside",max_height=maxh))
}
draw_labeled <- function(){
  old<-par(no.readonly=TRUE);on.exit({try(circos.clear(),silent=TRUE);par(old)},add=TRUE)
  par(fig=c(.04,.96,.20,.94),mar=c(0,0,0,0),family="sans",xpd=NA);draw_core()
  par(fig=c(0,1,.90,1),mar=c(0,0,0,0),new=TRUE,family="sans");plot.new();text(.5,.6,"Final morphology taxonomy",font=2,cex=.8)
  par(fig=c(0,1,0,.20),mar=c(0,0,0,0),new=TRUE,family="sans");plot.new()
  text(.5,.82,"Inner: HC K=4   Outer: GC resolution=0.50",cex=.55)
  text(.5,.61,sprintf("Agreement: %d/%d (%.2f%%)",sum(lab$HC_GC_consensus),n,100*mean(lab$HC_GC_consensus)),font=2,cex=.58)
  legend("bottom",legend=names(cols),fill=cols,border=NA,bty="n",ncol=4,cex=.48,x.intersp=.45)
}
draw_clean <- function(){on.exit(try(circos.clear(),silent=TRUE),add=TRUE);par(mar=c(0,0,0,0));draw_core()}

png1<-file.path(out,"17_final_HC_inner_GC_outer_circular_dendrogram_labeled.png")
pdf1<-file.path(out,"17_final_HC_inner_GC_outer_circular_dendrogram_labeled.pdf")
png2<-file.path(out,"18_final_HC_inner_GC_outer_circular_dendrogram_clean.png")
pdf2<-file.path(out,"18_final_HC_inner_GC_outer_circular_dendrogram_clean.pdf")
if(requireNamespace("ragg",quietly=TRUE))ragg::agg_png(png1,width=3.2,height=3.2,units="in",res=900,background="white") else png(png1,width=3.2,height=3.2,units="in",res=900,bg="white")
draw_labeled();dev.off();cairo_pdf(pdf1,width=3.2,height=3.2,bg="white");draw_labeled();dev.off()
if(requireNamespace("ragg",quietly=TRUE))ragg::agg_png(png2,width=1.5,height=1.5,units="in",res=900,background="white") else png(png2,width=1.5,height=1.5,units="in",res=900,bg="white")
draw_clean();dev.off();cairo_pdf(pdf2,width=1.5,height=1.5,bg="white");draw_clean();dev.off()

write.csv(data.frame(Circular_leaf_order=seq_len(n),MSN_unique_ID=ids,M_class=lab$M_class,GC_as_M_class=lab$GC_as_M_class,HC_GC_consensus=lab$HC_GC_consensus),file.path(out,"17_circular_dendrogram_leaf_order.csv"),row.names=FALSE)
cat("Circular dendrograms written to",out,"\n")
