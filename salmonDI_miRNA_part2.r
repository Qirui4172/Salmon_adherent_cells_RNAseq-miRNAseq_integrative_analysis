#!/usr/bin/env Rscript


#----------------------------------------------------------------------------------------------------
Usage<-function(){
	cat("\n\tUsage: Rscript salmonDI_miRNA_part2.r <readCounts.mx> <sampleInfo> <readCounts/gene> <nonZeroLib/gene> <controlGroup> <adjPvalue> <foldChange>","\n\n",

		"\tParameters (all required)","\n",
		"\t<readCounts.mx>			miRNA read counts matrix (\"miRNA_readcount.mx\", generated from \"salmonDI_miRNA_part1.sh\")","\n",
		"\t<sampleInfo>			Sample information file (\"sample.info\")","\n",
		"\t<readCounts/miRNA>		Minimal total read counts in all 12 libraries per miRNA, miRNAs with read counts less than this number will be filtered out (suggest: 100)","\n",
		"\t<nonZeroLib/miRNA>		Minimal non-zero libraries per miRNA, miRNAs with non-zero libraries less than this number will be filtered out (suggest: 5)","\n",
		"\t<controlGroup>			Control group, \"DI\" or \"HK\" (DI: distal intestine, HK: head kidney)","\n",
		"\t<adjPvalue>				BH adjusted pvalue for running DESeq2 (suggest: 0.05)","\n",
		"\t<foldChange>			Minimal fold change threshold for selecting DEGs (suggest: 2)","\n\n",

		"\tExample","\n",
        "\tRscript salmonDI_miRNA_part2.r miRNA_readcount.mx sample.info 100 5 DI 0.05 2","\n\n",

		"\tFunction","\n",
		"\tRun DESeq2, to find differentially expressed miRNAs and generate plots.","\n\n",

		"\tContact: Qirui Zhang (qirui.zhang@med.lu.se)","\n",
		"\tUpdated: 11-06-2020","\n\n"
	)
	quit()
}

args<-commandArgs(TRUE)
if (length(args)!=7){Usage()}

#----------------------------------------------------------------------------------------------------
# Load libraries
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Start analysis", "\n")
cat("Loading libraries ...", "\n")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pdf("salmonDI_miRNA_plots.pdf")

#----------------------------------------------------------------------------------------------------
cat("\n", "==================================================================================", "\n")
# read data
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Reading data ...", "\n")

readCounts<-read.table(args[1], header=T)
readCounts2<-as.data.frame(lapply(readCounts,as.integer))
rownames(readCounts2)<-rownames(readCounts)
readCounts<-readCounts2
sampleInfo<-read.table(args[2], header=T)
sampleInfo$Replicate<-as.factor(sampleInfo$Replicate)
colnames(readCounts)<-sampleInfo$Sample
readCounts<-readCounts[which(rowSums(readCounts) >= as.numeric(args[3])),]
readCounts<-readCounts[which(rowSums(readCounts != 0) >= as.numeric(args[4])),]

#----------------------------------------------------------------------------------------------------
cat("\n", "==================================================================================", "\n")
# run DESeq2
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Generating metadata, runing DESeq2, and normalizing read counts ...", "\n")

dds<-DESeqDataSetFromMatrix(countData=readCounts,colData=sampleInfo,design=~Tissue)
dds$Tissue<-relevel(dds$Tissue, args[5])
dds<-DESeq(dds)

normalized.counts<-counts(dds, normalized=TRUE)
write.table(as.data.frame(normalized.counts), "salmonDI_miRNA_baseMean.tsv", row.names=T, col.names=T, quote=F, sep="\t")

#----------------------------------------------------------------------------------------------------
cat("\n", "==================================================================================", "\n")
# extract comparison results
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Extracting comparison results ...", "\n")

adjP<-as.numeric(args[6])
FC<-as.numeric(args[7])
res<-results(dds,alpha=adjP)
summary(res)
demir.num<-sum(res$padj<adjP & abs(res$log2FoldChange)>=log2(FC), na.rm=TRUE)
cat("\t","Total DEmiRNA num: ",demir.num,"\n")
up.num<-sum(res$padj<adjP & res$log2FoldChange>=log2(FC), na.rm=TRUE)
cat("\t","Up-regulated DEmiRNA num: ",up.num,"\n")
down.num<-sum(res$padj<adjP & -(res$log2FoldChange)>=log2(FC), na.rm=TRUE)
cat("\t","Down-regulated DEmiRNA num: ",down.num,"\n\n")
demir.all<-res[which(res$padj<adjP & abs(res$log2FoldChange)>=log2(FC)),]
write.table(as.data.frame(demir.all), "salmonDI_DEmiRNA_list.tsv", row.names=T, col.names=T, quote=F, sep="\t")

#----------------------------------------------------------------------------------------------------
cat("\n", "==================================================================================", "\n")
# make sample plots
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Sample plots:", "\n")

# transform data
cat("Transforming data with DESeq2::rlog ...", "\n")
vst<-varianceStabilizingTransformation(dds, blind=FALSE)

# plot dispersion estimate
cat("Dispersion estimate plot ...", "\n")
plotDispEsts(dds)

# PCA plots
cat("PCA plot ...", "\n")
PCA_Plot<-function(vst, DI.color, HK.color){
	data<-plotPCA(vst, intgroup="Tissue", returnData=TRUE)
	percentVar<-round(100*attr(data, "percentVar"))
	DI.color<-DI.color
	HK.color<-HK.color
	ggplot(data, aes(PC1, PC2, color=Tissue))+geom_point(size=5)+scale_colour_manual(values=c("DI"={DI.color},"HK"={HK.color}))+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(size=1),text=element_text(size=20),axis.text=element_text(size=15),legend.text=element_text(size=15))+xlab(paste0("PC1: ",percentVar[1],"% variance"))+ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed(ratio=2)
}
PCA_Plot(vst, "salmon", "firebrick")

# correlation heatmap
cat("Sample correlation heatmap ...", "\n")
sampleDists<-dist(t(assay(vst)))
sampleDistMatrix<-as.matrix(sampleDists)
rownames(sampleDistMatrix)<-vst$Sample
colnames(sampleDistMatrix)<-NULL
colors<-colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)

#----------------------------------------------------------------------------------------------------
cat("\n", "==================================================================================", "\n")
# DEG plots
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "DEG plots:", "\n")

# MA-plot
cat("MA plot ...", "\n")
plotMA(res, main="MA plot", ylim=c(-10,10))

# volcano plot
cat("volcano plot ...", "\n")
VolcanoPlot<-function(res, down.color, up.color){
	# divide groups
	volcano<-as.data.frame(res)
	volcano$significant<-as.factor(ifelse(!is.na(volcano$padj) & volcano$padj < adjP & abs(volcano$log2FoldChange)>=log2(FC), ifelse(volcano$log2FoldChange>=log2(FC), "Up", "Down"), "No"))
	volcano$padj<-ifelse(is.na(volcano$padj),1,volcano$padj)
	volcano$group<-rep(0)
	volcano[which(volcano$padj >= adjP | (volcano$padj < adjP & abs(volcano$log2FoldChange) < log2(FC))),8]=rep(1)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-20 & volcano$log2FoldChange <= -1 & volcano$log2FoldChange >= -10),8]=rep(2)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-20 & volcano$log2FoldChange <= 10 & volcano$log2FoldChange >= 1),8]=rep(3)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-20 & volcano$log2FoldChange < -10),c("log2FoldChange", "group")]=list(log2FoldChange = -10, group=4)
	volcano[which(volcano$padj < 1e-20 & volcano$log2FoldChange <= -1 & volcano$log2FoldChange >= -10),c("padj","group")]=list(padj = 1e-20, group=4)
	volcano[which(volcano$padj < 1e-20 & volcano$log2FoldChange < -10),c("log2FoldChange", "padj", "group")]=list(log2FoldChange = -10, padj = 1e-20, group=4)
	volcano[which(volcano$padj < adjP & volcano$padj >= 1e-20 & volcano$log2FoldChange > 10),c("log2FoldChange", "group")]=list(log2FoldChange = 10, group=5)
	volcano[which(volcano$padj < 1e-20 & volcano$log2FoldChange <= 10 & volcano$log2FoldChange >= 1),c("padj", "group")]=list(padj = 1e-20, group=5)
	volcano[which(volcano$padj < 1e-20 & volcano$log2FoldChange > 10),c("padj", "log2FoldChange", "group")]=list(log2FoldChange = 10, padj = 1e-20, group=5)

	# plot
	down.color<-down.color
	up.color<-up.color
	p<-ggplot(volcano,aes(log2FoldChange,-log10(padj)))+geom_point(data=volcano[which(volcano$group==1),],color="gray",alpha=0.75)+geom_point(data=volcano[which(volcano$group==2),],color={down.color},alpha=0.75)+geom_point(data=volcano[which(volcano$group==3),],color={up.color},alpha=0.75)+geom_point(data=volcano[which(volcano$group==4),],shape=2,color={down.color},alpha=0.75)+geom_point(data=volcano[which(volcano$group==5),],shape=2,color={up.color},alpha=0.75)
	p+labs(title="Volcano plot",x="log2FoldChange",y="-log10(padj)")+geom_hline(yintercept=-log10(adjP),linetype=2,color="gray80")+geom_vline(xintercept=c(-log2(FC),log2(FC)),linetype=2,color="gray80")+xlim(-10,10)+ylim(0,20)+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(size=1))
}
VolcanoPlot(res, "firebrick", "salmon")

# DEmiRNA heatmap
cat("DEmiRNA heatmap ...", "\n")
HeatmapPlot<-function(demir.all, DI_anno.color, HK_anno.color, heatmap.color){
	vst.demir<-vst[rownames(demir.all),]
	anno.label<-data.frame(Tissue=c(rep("DI", 6),rep("HK", 6)))
	rownames(anno.label)<-rownames(colData(vst))
	anno.color<-list(Tissue=c(DI=DI_anno.color, HK=HK_anno.color))
	pheatmap(assay(vst.demir),scale="row",main="Heatmap of DEmiRNAs",color={heatmap.color},cluster_cols=F,show_rownames=T,annotation_col=anno.label,annotation_colors=anno.color,annotation_names_col=F,cellwidth=15,border_color=NA)
}
DI_anno.color<-"salmon"
HK_anno.color<-"firebrick"
heatmap.color<-colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
HeatmapPlot(demir.all, DI_anno.color, HK_anno.color, heatmap.color)

#----------------------------------------------------------------------------------------------------
cat("\n", "==================================================================================", "\n")
time<-format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Done with analysis!", "\n")
cat("\nGenerated files:\n")
cat("\"salmonDI_DEmiRNA_list.tsv\"\t\tDifferentitally expressed miRNAs\n")
cat("\"salmonDI_miRNA_baseMean.tsv\"\tDESeq2-normalized read counts\n")
cat("\"salmonDI_miRNA_plots.pdf\"\t\tPlot file\n\n")
cat("Oppa, we got a lot of DEmiRNAs ^(00)^\n\n")

