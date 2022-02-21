
#==========================================================================================================================
# General information

# Brief introduction:
# There in-house scripts were wrote for integrative analysis of mRNA-seq and miRNA-seq of adherent cells from Salmon distal intestine and head kidney. It includes two parts, part1 "salmonDI_miRNA_part1.sh" for trimming adapters, aligning to genome, controlling the quality, and quantifying read counts; and part2 "salmonDI_miRNA_part2.r" for comparing differential transcripts and generating plots. For more detailed introductions please see below.

# Author:
# Qirui Zhang (qirui.zhang@med.lu.se)

# Citation:
# Y Park, Q Zhang, J M.O. Fernandes, and V Kiron. Transcriptome analysis of mRNA and miRNA in adherent intestinal cells of Atlantic salmon, Salmo salar. XXX, 2020, XXX


#==========================================================================================================================
# Experimental design

# Leukocyte cells were isolated from Atlantic salmon distal intestine (DI) and head kidney (HK) tissues, and were shortly cultured in vitro to obtain adherent cells. Six fish were used for generating DI and HK samples, resulting in 6 DI replicates and 6 HK replicates. Adherent cells were harvested and subjected to extract total RNA and small RNAs, construct libraries, and run RNA-seq and small RNA-seq on an Illumina NextSeq 500 sequencer. miRNA expression profiles were compared between DI and HK replicates, miRNA-seq and RNA-seq were integratively analysed to find significantly correlated miRNA-target transcript pairs. More details can be found in "Y Park, Q Zhang, J M.O. Fernandes, and V Kiron. XXX, 2020, XXX"


#==========================================================================================================================
# Running environment

# The script "salmonDI_miRNA_part1.sh" was originally run in author's private work directory "/home2/park/salmonDI_mirna", raw fastq files were stored in directory "/home2/park/rawData/flowcell3_mirna/fastqfiles", Salmon genome, miRNA hairpin and mature fasta files were stored in directory "/home2/park/genome/for_miRNA". Other users have to modify these THREE directory variables, "WorkDir", "RawReads", and "Genome", in line5-7 of script to their own directory for custom use. Barcode sequence are stored in line8 of the script, other users also have to modify it to their own barcode sequence.

# The script "salmonDI_miRNA_part2.r" can be run directly in other users' own directory, please keep the input files ("miRNA_readcount.mx", and "sample.info") in the same directory. The file "miRNA_readcount.mx" has to be firstly generated from script "salmonDI_miRNA_part1.sh".

#==========================================================================================================================
# Part1 script: "salmonDI_mRNA_part1.sh"


############## Usage of the script1 ################

# Usage: bash salmonDI_miRNA_part1.sh --step=All
#
#	Parameters
#	-s|--step	(required)	Which step(s) to run (please type step NUMBERS separated by comma, see "Steps" and "Examples" below)
#	-h|--help	(optional)	Return help information
#
#	Steps
#	1		Quality control of raw fastq
#	2		Trim adapter
#	3		Quality control of trimmed fastq
#	4		Align to salmon genome (ICSASG_v2)
#	5		Generate read counts matrix
#	All		Run all above steps
#
#	Examples
#	bash salmonDI_miRNA_part1.sh -s3,4,5 (or --step=3,4,5)
#	bash salmonDI_miRNA_part1.sh -sAll (or --step=All)
#	bash salmonDI_miRNA_part1.sh -h (for help information)
#
#	Function
#	Analyze salmon distal intestine and head kidney miRNA-seq data, including control the quality of raw fastq; trim adapter and control the quality again; align to genome; and generate read counts matrix.


########### Input and output of the script1 #############

# Input:
# Step "-s/--step" is required, users have to specify which step(s) to run.

# Output:
# Folders of alignment/, cleanReads/, qualityControl/, qualityControl/cleanReads/, qualityControl/rawReads/, and quantity/ are created under ${WorkDir}
# Clean reads are stored in ${WorkDir}/cleanReads/
# Alignments are stored in ${WorkDir}/alignment/
# Quality control results of raw reads and clean reads are stored in the corresponding subdirectories under ${WorkDir}/qualityControl/.
# "miRNA_readcount.mx" is the newly generated raw read counts matrix of all 12 samples (6 DI and 6 HK replicates), which is to be used in the script "salmonDI_miRNA_part2.r".



#==========================================================================================================================
# Part2 script: "salmonDI_miRNA_part2.r"


############## Usage of the script2 ################

# Usage: Rscript salmonDI_miRNA_part2.r <readCounts.mx> <sampleInfo> <readCounts/gene> <nonZeroLib/gene> <controlGroup> <adjPvalue> <foldChange>

#	Parameters (all required)
#	<readCounts.mx>			miRNA read counts matrix ("miRNA_readcount.mx", generated from "salmonDI_miRNA_part1.sh")
#	<sampleInfo>			Sample information file ("sample.info")
#	<readCounts/miRNA>		Minimal total read counts in all 12 libraries per miRNA, miRNAs with read counts less than this number will be filtered out (suggest: 100)
#	<nonZeroLib/miRNA>		Minimal non-zero libraries per miRNA, miRNAs with non-zero libraries less than this number will be filtered out (suggest: 5)
#	<controlGroup>			Control group, "DI" or "HK" (DI: distal intestine, HK: head kidney)
#	<adjPvalue>				BH adjusted pvalue for running DESeq2 (suggest: 0.05)
#	<foldChange>			Minimal fold change threshold for selecting DEGs (suggest: 2)
#
#	Example
#	Rscript salmonDI_miRNA_part2.r miRNA_readcount.mx sample.info 100 5 DI 0.05 2
#
#	Function
#	Run DESeq2, to find differentially expressed miRNAs and generate plots.



########### Input and output of the script2 #############

# Input:
# "miRNA_readcount.mx" is the RAW read counts matrix of miRNAs generated from "salmonDI_miRNA_part1.sh".
# "sample.info" is a manually-prepared sample information file.

# Output:
# "salmonDI_miRNA_baseMean.tsv" is the DESeq2-normalized read counts of all 12 samples (6 DI and 6 HK replicates).
# "salmonDI_DEmiRNA_list.tsv" is the full list of DEmiRNAs with an absolute fold change > foldChange, and BH adjusted < adjPvalue.
# "salmonDI_miRNA_plots.pdf" is the plot collection file, including sample correlation heatmap, PCA, MA plot, normalization dispersion estimate plot, DEmiRNA volcano plot, and DEmiRNA heatmap.


#==========================================================================================================================

