#!/bin/bash


Time=`date "+%Y-%m-%d %H:%M:%S"`
WorkDir="/home2/park/salmonDI_mirna"
RawReads="/home2/park/rawData/flowcell3_mirna/fastqfiles"
Genome="/home2/park/genome/for_miRNA"
SampleBarcode=("DI1:CCGTCC" "DI2:GTAGAG" "DI3:GTCCGC" "DI4:GTGAAA" "DI5:GTGGCC" "DI6:GTTTCG" "HK1:ATCACG" "HK2:GATCAG" "HK3:CTTGTA" "HK4:AGTCAA" "HK5:AGTTCC" "HK6:ATGTCA")

#-------------------------------------------------------------------------------------------------
Args=`getopt -o s:h --long step:,help -- "$@"`
eval set -- "$Args"

Usage()
{
	echo -e "\n\tUsage: bash salmonDI_miRNA_part1.sh --step=All"

	echo -e "\n\tParameters"
	echo -e "\t-s|--step	(required)	Which step(s) to run (please type step NUMBERS separated by comma, see \"Steps\" and \"Examples\" below)"
	echo -e "\t-h|--help	(optional)	Return help information\n"

	echo -e "\tSteps"
	echo -e "\t1		Quality control of raw fastq"
	echo -e "\t2		Trim adapter"
	echo -e "\t3		Quality control of trimmed fastq"
	echo -e "\t4		Align to salmon genome (ICSASG_v2)"
	echo -e "\t5		Generate read counts matrix"
	echo -e "\tAll		Run all above steps\n"

	echo -e "\tExamples"
	echo -e "\tbash salmonDI_miRNA_part1.sh -s3,4,5 (or --step=3,4,5)"
	echo -e "\tbash salmonDI_miRNA_part1.sh -sAll (or --step=All)"
	echo -e "\tbash salmonDI_miRNA_part1.sh -h (for help information)\n"

	echo -e "\tFunction"
	echo -e "\tAnalyze salmon distal intestine and head kidney miRNA-seq data, including control the quality of raw fastq; trim adapter and control the quality again; align to genome; and generate read counts matrix.\n"

	echo -e "\tContact: Qirui Zhang (qirui.zhang@med.lu.se)"
	echo -e "\tDate: 16-11-2019\n\n"
}

#-------------------------------------------------------------------------------------------------
# Inspect work folder
#if [[ $PWD != $WorkDir ]];then echo -e "\n\t*** WARNING: You're attempting to run the script in the wrong path, please run the script in the path: \"$WorkDir/\"!";Usage;exit;fi


if ([[ $@ != *"-s"* ]]);then	# didn't specify "-s"
	if ([[ $@ != *"-h"* ]] && [[ $# -ne 1 ]] );then # didn't specify "-s" or "-h" but input some other unclear values
		echo -e "\n\t*** WARNING: Incorrect parameters, please check again!\n"
	fi
	Usage;
	exit;
fi

# Prepare steps
echo -e "Checking steps ..."
#if [[ $# == 1 ]];then Usage;exit;fi
while [[ $# -ne 0 ]]
do
	case $1 in
		-s|--step)
		case $2 in
			"")
			echo -e "\n\t*** WARNING: Step numbers are missing! please specify step numbers after \"-s\" or \"--step=\"!\n"
			Usage
			shift
			exit
			;;
			All)
			step=$2
			shift 2
			;;
			[1-5]*)
			Steps=($(echo $2| sed 's/,/\t/g'))
			shift 2
			;;
		esac
		;;
        -h|--help)
        Usage
        shift
        exit
        ;;
        --)
        shift
        break
        ;;
    esac
done

mkdir -p alignment cleanReads qualityControl qualityControl/cleanReads qualityControl/rawReads quantity

#-------------------------------------------------------------------------------------------------
echo -e "\n======================================================================================"
echo -e "Step1: Quality control of raw fastq files\n"

# Step1: Quality control of raw fastq files
if [[ ${Steps[@]} =~ 1 ]] || [[ $step == "All" ]];then
	echo -e "$Time\tStep1 starts"

	for sample in ${Samples[@]}
	do
		echo -e "Controlling quality of raw fastq of sample ${sample} ..."
		fastqc -o qualityControl/rawReads -q $RawReads/${sample}*fastq.gz
	done

	echo -e "$Time\tStep1 done!"
else
	echo -e "Step1 skipped"
fi

#-------------------------------------------------------------------------------------------------
echo -e "\n======================================================================================"
echo -e "Step2: Trim adapters\n"

# Step2: Trim adapters
A1_1="AACTCCAGTCAC"
A1_2="ATCTCGTATGCCGTCTTCTGCTTG"
A2="TGGAATTCTCGGGTGCCAAGG"
G1="GTTCAGAGTTCTACAGTCCGACGATC"

if [[ ${Steps[@]} =~ 2 ]] || [[ $step == "All" ]];then
    echo -e "$Time\tStep2 starts"

	for samplebarcode in ${SampleBarcode[@]}
	do
		arr=(${samplebarcode//:/ })
		sample=${arr[0]}
		barcode=${arr[1]}
		echo -e "Trimming sample ${sample} ..."
		adapt3end1=${A1_1}${barcode}${A1_2}
		adapt3end2=${A2}
		adapt5end=${G1}
		rawreads=`echo ${RawReads}/${sample}*fastq`

		cutadapt -a ${adapt3end1} --no-indels -o cleanReads/${sample}.trim1.fastq ${rawreads}
		cutadapt -a ${adapt3end2} --no-indels -o cleanReads/${sample}.trim2.fastq cleanReads/${sample}.trim1.fastq
		cutadapt -g ${adapt5end} --no-indels -o cleanReads/${sample}.trim3.fastq cleanReads/${sample}.trim2.fastq
		cutadapt -u 4 -u -4 --minimum-length=15 --maximum-length=35 -o cleanReads/${sample}.trim4.fastq cleanReads/${sample}.trim3.fastq
		fastq_quality_filter -q 20 -p 85 -i cleanReads/${sample}.trim4.fastq  -o cleanReads/${sample}.fastq

		rm cleanReads/${sample}.trim*fastq
	done

	echo -e "$Time\tStep2 done!"
else
	echo -e "Step2 skipped"
fi

#-------------------------------------------------------------------------------------------------
echo -e "\n======================================================================================"
echo -e "Step3: Quality control of trimmed fastq\n"

# Step3: Quality control of trimmed fastq
if [[ ${Steps[@]} =~ 3 ]] || [[ $step == "All" ]];then
	echo -e "$Time\tStep3 starts"

	for sample in ${Samples[@]}
	do
		echo -e "Controlling quality of trimmed sample ${sample} ..."
		fastqc -o qualityControl/cleanReads -q cleanReads/${sample}.fastq
	done

	echo -e "$Time\tStep3 done!"
else
	echo -e "Step3 skipped"
fi

#-------------------------------------------------------------------------------------------------
echo -e "\n======================================================================================"
echo -e "Step4: Align to Atlantic salmon genome (ICSASG_v2)\n"

# Step4: Alignment
if [[ ${Steps[@]} =~ 4 ]] || [[ $step == "All" ]];then
	echo -e "$Time\tStep4 starts"
	i=0

	bowtie-build ${Genome}/ssa_genome.fa ${Genome}/ssa_genome

	for samplebarcode in ${SampleBarcode[@]}
	do
		arr=(${samplebarcode//:/ })
		sample=${arr[0]}
		i=$(($i+1))
		group=${sample:0:2}
		echo -e "${WorkDir}/cleanReads/${sample}.fastq\t${sample}" >> alignment/mapper_${group}.txt

		if [[ i==6 ]];then
			mapper.pl alignment/mapper_${group}.txt -e -d -h -i -j -l 18 -m -n -o 16 -p ${Genome}/ssa_genome -s alignment/${group}.pool.fa -t alignment/${group}.pool.arf
			i=0
		fi
	done

	echo -e "$Time\tStep4 done! Alignment statistics is in alignment/bowtie.log."
else
	echo -e "Step4 skipped"
fi

#-------------------------------------------------------------------------------------------------
echo -e "\n======================================================================================"
echo -e "Step5: Generate read counts matrix\n"

# Step5: Generate read counts matrix
if [[ ${Steps[@]} =~ 5 ]] || [[ $step == "All" ]];then
	echo -e "$Time\tStep5 Starts"

	for group in "DI" "HK"
	do
		echo -e "Generating read counts of group ${group} ..."
		mkdir -p quantity/${group}; cd quantity/${group}
		quantifier.pl -p ${Genome}/ssa_hairpin.fa -m ${Genome}/ssa_mature.fa -P -r alignment/${group}_pool.fa
		cut -f1,5-10 *csv|sort -k1,1|awk 'BEGIN {FS="\t";OFS="\t";mir="";rep1=0;rep2=0;rep3=0;rep4=0;rep5=0;rep6=0;i=0;} NR==1 {print $0;next;} NR==2 {mir=$1;rep1=$2;rep2=$3;rep3=$4;rep4=$5;rep5=$6;rep6=$7;i=1;next;} $1==mir {rep1+=$2;rep2+=$3;rep3+=$4;rep4+=$5;rep5+=$6;rep6+=$7;i+=1;next;} $1!=mir {print mir,rep1/i,rep2/i,rep3/i,rep4/i,rep5/i,rep6/i;mir=$1;rep1=$2;rep2=$3;rep3=$4;rep4=$5;rep5=$6;rep6=$7;i=1;next;} END {print mir,rep1,rep2,rep3,rep4,rep5,rep6;}' > ${group}_readcounts.txt
		cd ${WorkDir}/quantity
	done
	paste DI/DI_readcounts.txt HK/HK_readcounts.txt |cut -f1-7,9-14 > miRNA_readcounts.mx
	cd ${WorkDir}

	echo -e "$Time\tStep5 done!"
else
	echo -e "Step5 skipped"
fi

#-------------------------------------------------------------------------------------------------
echo -e "\n======================================================================================"
echo -e "$Time\tDone with the whole analysis!\n"

