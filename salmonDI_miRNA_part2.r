#!/usr/bin/env Rscript

##===========================================================================================================
Usage <- function() {
	cat("\n\tUsage: Rscript salmonDI_miRNA_part2.r <readCounts.mx> <sampleInfo> <readCounts/miRNA> <nonZeroLib/miRNA> <controlGroup> <adjPvalue> <foldChange>", "\n\n",

		"\tParameters (all required)","\n",
		"\t<readCounts.mx>          miRNA read counts matrix (\"miRNA_readcount.mx\", generated from \"salmonDI_miRNA_part1.sh\")","\n",
		"\t<sampleInfo>             Sample information file (\"sample.info\")","\n",
		"\t<readCounts/miRNA>       Minimal total read counts in all 12 libraries per miRNA, miRNAs with read counts less than this number will be filtered out (suggest: 100)","\n",
		"\t<nonZeroLib/miRNA>       Minimal non-zero libraries per miRNA, miRNAs with non-zero libraries less than this number will be filtered out (suggest: 5)","\n",
		"\t<controlGroup>           Control group, \"DI\" or \"HK\" (DI: distal intestine, HK: head kidney)","\n",
		"\t<adjPvalue>              BH adjusted pvalue for running DESeq2 (suggest: 0.05)","\n",
		"\t<foldChange>             Minimal fold change threshold for selecting DEGs (suggest: 2)","\n\n",

		"\tExample","\n",
        "\tRscript salmonDI_miRNA_part2.r miRNA_readcount.mx sample.info 100 5 DI 0.05 2","\n\n",

		"\tFunction","\n",
		"\tRun DESeq2, to find differentially expressed miRNAs and generate plots.","\n\n",

		"\tContact: Qirui Zhang (qirui.zhang@med.lu.se)","\n",
		"\tUpdated: 11-06-2020", "\n\n"
	)
	quit()
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
	Usage()
}

read_count_threshold <- as.numeric(args[3])
non_zero_library_threshold <- as.numeric(args[4])
adjP <- as.numeric(args[6])
FC <- as.numeric(args[7])
if (any(!is.finite(c(read_count_threshold, non_zero_library_threshold, adjP, FC))) ||
	read_count_threshold < 0 || non_zero_library_threshold < 1 || adjP <= 0 || adjP >= 1 || FC <= 0) {
	stop("Numeric thresholds are invalid.")
}

##===========================================================================================================
## Load libraries
time <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
cat(time, "Start analysis", "\n")
cat("Loading libraries ...", "\n")
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
pdf("salmonDI_miRNA_plots.pdf")

##===========================================================================================================
## Read data
cat("\n", "==================================================================================", "\n")
time <- format(Sys.time(), format='%Y-%m-%d %H:%M:%S')
cat(time, "Reading data ...", "\n")

readCounts <- read.table(args[1], header = TRUE, row.names = 1, check.names = FALSE)
readCounts <- as.data.frame(lapply(readCounts, function(counts) as.integer(round(counts))))
sampleInfo <- read.table(args[2], header = TRUE, stringsAsFactors = FALSE)

if (ncol(readCounts) != nrow(sampleInfo)) {
	stop("The number of read-count columns does not match the sample metadata.")
}
if (anyDuplicated(sampleInfo$Sample) || anyDuplicated(colnames(readCounts))) {
	stop("Sample IDs must be unique in both input files.")
}
if (!setequal(colnames(readCounts), sampleInfo$Sample)) {
	stop("Sample IDs in the count matrix and metadata do not match.")
}
sampleInfo <- sampleInfo[match(colnames(readCounts), sampleInfo$Sample), , drop = FALSE]
rownames(sampleInfo) <- sampleInfo$Sample
sampleInfo$Tissue <- factor(sampleInfo$Tissue)
sampleInfo$Replicate <- factor(sampleInfo$Replicate)
readCounts <- readCounts[rowSums(readCounts) >= read_count_threshold, , drop = FALSE]
readCounts <- readCounts[rowSums(readCounts != 0) >= non_zero_library_threshold, , drop = FALSE]
if (nrow(readCounts) == 0) {
	stop("No miRNAs remain after filtering.")
}

##===========================================================================================================
## Run DESeq2
cat("\n", "==================================================================================", "\n")
time <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
cat(time, "Generating metadata, runing DESeq2, and normalizing read counts ...", "\n")

dds <- DESeqDataSetFromMatrix(countData = readCounts, colData = sampleInfo, design =~ Tissue)
dds$Tissue <- relevel(dds$Tissue, args[5])
dds <- DESeq(dds)

normalized.counts <- counts(dds, normalized = TRUE)
write.table(as.data.frame(normalized.counts), "salmonDI_miRNA_baseMean.tsv", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

##===========================================================================================================
## Extract comparison results
cat("\n", "==================================================================================", "\n")
time <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
cat(time, "Extracting comparison results ...", "\n")

res <- results(dds, alpha = adjP)
summary(res)
demir.num <- sum(res$padj < adjP & abs(res$log2FoldChange) >= log2(FC), na.rm = TRUE)
cat("\t","Total DEmiRNA num: ",demir.num,"\n")
up.num <- sum(res$padj < adjP & res$log2FoldChange >= log2(FC), na.rm = TRUE)
cat("\t","Up-regulated DEmiRNA num: ",up.num,"\n")
down.num <- sum(res$padj < adjP & -(res$log2FoldChange) >= log2(FC), na.rm = TRUE)
cat("\t","Down-regulated DEmiRNA num: ",down.num,"\n\n")
demir.all <- res[which(res$padj < adjP & abs(res$log2FoldChange) >= log2(FC)), , drop = FALSE]
write.table(as.data.frame(demir.all), "salmonDI_DEmiRNA_list.tsv", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

##===========================================================================================================
## Generate plots
cat("\n", "==================================================================================", "\n")
time <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
cat(time, "Sample plots:", "\n")

# transform data
cat("Transforming data with DESeq2::rlog ...", "\n")
vst <- varianceStabilizingTransformation(dds, blind = FALSE)

# plot dispersion estimate
cat("Dispersion estimate plot ...", "\n")
plotDispEsts(dds)

# PCA plots
cat("PCA plot ...", "\n")

PCA_Plot <- function(vst, DI.color, HK.color) {
	data <- plotPCA(vst, intgroup = "Tissue", returnData = TRUE)
	percentVar <- round(100 * attr(data, "percentVar"))
	p <- ggplot(data, aes(PC1, PC2, color = Tissue)) +
		geom_point(size = 5) +
		scale_colour_manual(values = c(DI = DI.color, HK = HK.color)) +
		theme_bw() +
		theme(panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.border = element_rect(size = 1), 
			text = element_text(size = 20),
			axis.text = element_text(size = 15), 
			legend.text = element_text(size = 15)) +
		xlab(paste0("PC1: ", percentVar[1], "% variance")) +
		ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		coord_fixed(ratio = 2)
	
	print(p)
	return(invisible(p))
}
PCA_Plot(vst, "indianred1", "turquoise3")
PCA_Plot(vst, "salmon", "firebrick")

# correlation heatmap
cat("Sample correlation heatmap ...", "\n")
sampleDists <- dist(t(assay(vst)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vst)
colnames(sampleDistMatrix) <- colnames(vst)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, 
	clustering_distance_rows = sampleDists,
	clustering_distance_cols = sampleDists, 
	col = colors)

##===========================================================================================================
## DEG plots
cat("\n", "==================================================================================", "\n")
time <- format(Sys.time(), format = '%Y-%m-%d %H:%M:%S')
cat(time, "DEG plots:", "\n")

# MA-plot
cat("MA plot ...", "\n")
plotMA(res, main = "MA plot", ylim = c(-10,10))

# volcano plot
cat("volcano plot version1 ...", "\n")

VolcanoPlotV1 <- function(res, down.color, up.color) {
	volcano <- as.data.frame(res)
	volcano$significant <- as.factor(ifelse(
		!is.na(volcano$padj) & volcano$padj < adjP & abs(volcano$log2FoldChange) >= log2(FC),
		ifelse(volcano$log2FoldChange >= log2(FC), "Up", "Down"), "No"
	))
	volcano$padj[is.na(volcano$padj)] <- 1
	volcano$group <- 1L
	volcano$group[volcano$padj < adjP & volcano$log2FoldChange <= -log2(FC)] <- 2L
	volcano$group[volcano$padj < adjP & volcano$log2FoldChange >= log2(FC)] <- 3L
	volcano$log2FoldChange[volcano$group == 2 & volcano$log2FoldChange < -10] <- -10
	volcano$log2FoldChange[volcano$group == 3 & volcano$log2FoldChange > 10] <- 10
	volcano$group[volcano$padj < 1e-20 & volcano$log2FoldChange <= -log2(FC)] <- 4L
	volcano$group[volcano$padj < 1e-20 & volcano$log2FoldChange >= log2(FC)] <- 5L
	volcano$padj[volcano$padj < 1e-20] <- 1e-20

	# plot
	p <- ggplot(volcano, aes(log2FoldChange, -log10(padj))) +
		geom_point(data = volcano[volcano$group == 1, ], color = "gray", alpha = 0.75) +
		geom_point(data = volcano[volcano$group == 2, ], color = down.color, alpha = 0.75) +
		geom_point(data = volcano[volcano$group == 3, ], color = up.color, alpha = 0.75) +
		geom_point(data = volcano[volcano$group == 4, ], shape = 2, color = down.color, alpha = 0.75) +
		geom_point(data = volcano[volcano$group == 5, ], shape = 2, color = up.color, alpha = 0.75)
	
	p <- p + 
		labs(title = "Volcano plot", x = "log2FoldChange", y = "-log10(padj)") +
		geom_hline(yintercept = -log10(adjP), linetype = 2, color = "gray80") +
		geom_vline(xintercept = c(-log2(FC), log2(FC)), linetype = 2, color = "gray80") +
		xlim(-10, 10) + 
		ylim(0, 20) + 
		theme_bw() +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(size = 1))
	
	print(p)
	return(invisible(p))
}

VolcanoPlotV1(res, "royalblue3", brewer.pal(11,"RdYlBu")[2])
VolcanoPlotV1(res, "firebrick", "salmon")

# DEmiRNA heatmap
cat("DEmiRNA heatmap ...", "\n")
HeatmapPlot <- function(demir.all, DI_anno.color, HK_anno.color, heatmap.color) {
	if (nrow(demir.all) == 0) {
		plot.new()
		text(0.5, 0.5, "No differentially expressed miRNAs")
		return(invisible(NULL))
	}
	vst.demir <- vst[rownames(demir.all), , drop = FALSE]
	anno.label <- data.frame(Tissue = colData(vst)$Tissue)
	rownames(anno.label) <- colnames(vst)
	anno.color <- list(Tissue = c(DI = DI_anno.color, HK = HK_anno.color))
	pheatmap(assay(vst.demir), 
		scale = "row", 
		main = "Heatmap of DEmiRNAs",
		color = heatmap.color, 
		cluster_cols = FALSE, 
		show_rownames = TRUE,
		annotation_col = anno.label, 
		annotation_colors = anno.color,
		annotation_names_col = FALSE, 
		cellwidth = 15, 
		border_color = NA)
}

DI_anno.color<-brewer.pal(8,"Dark2")[1]
HK_anno.color<-brewer.pal(8,"Dark2")[2]
heatmap.color<-colorRampPalette(c(rev(brewer.pal(9,"Blues")[c(2:9)]), brewer.pal(9, "OrRd")[c(2:9)]))(100)
HeatmapPlot(demir.all, DI_anno.color, HK_anno.color, heatmap.color)

DI_anno.color<-"salmon"
HK_anno.color<-"firebrick"
heatmap.color<-colorRampPalette(brewer.pal(9,"YlGnBu"))(100)
HeatmapPlot(demir.all, DI_anno.color, HK_anno.color, heatmap.color)
dev.off()

##===========================================================================================================
cat("\n", "==================================================================================", "\n")
time <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
cat(time, "Done with analysis!", "\n")
