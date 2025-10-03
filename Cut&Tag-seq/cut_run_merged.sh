#!/bin/bash
#ChIPseq.sh
#This is a pipeline for ChIP-seq analysis for merged data
#Requirement: conda activate daily
#Current version 1: J.Q.Z, 2021-10-27

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*                     PipeRiboseq: pipeline for ChIP-seq analysis.                         *"
  echo "*                            Version 2, 2021-10-28, Q.Z                                    *"
  echo "* Usage: `basename $0`                                                                     *"
  echo "*        Required (paired end data):                                                       *"
  echo "*                  -n [The prefix samplename]                                              *"
  echo "*                  -o [The outputpath]                                                     *"
  echo "*        Optional: -p [Number of CPUs, default=24]                                         *"
  echo "*                  -I [default=Not]                                                        *"
  echo "*                  -c [default=0]                                                          *"
  echo "* Inputs: The raw fastq files                                                              *"
  echo "* Run: Default to run trim_galore,fastqc,bowtie2,sam2bam,rmdup & deeptools                 *"
  echo "*      Figures will be generated in /plots folder, and bigWig files in /tracks folder      *"
  echo "* Outputs: All output files will be generated in the same folder as the pipeline submitted *"
  echo "* This pipeline requires 'conda activate daily'                                            *"
  echo "********************************************************************************************"
  exit 1
fi

#Get parameters
sample="unassigned"
out="unassigned"
CPU=24
input="Not"
class=0
for arg in "$@"
do
 if [[ $arg == "-n" ]]
  then
    sample=${array[$counter+1]}
    echo '   Sample name: '$sample
 elif [[ $arg == "-o" ]]
  then
    out=${array[$counter+1]}
    echo '   The output of all files: '$out
 elif [[ $arg == "-p" ]]
  then
    CPU=${array[$counter+1]}
    echo '   CPU: '$CPU
 elif [[ $arg == "-I" ]]
  then
    input=${array[$counter+1]}
    echo '   INPUT files: '$input
 elif [[ $arg == "-c" ]]
  then
    class=${array[$counter+1]}
    echo '   Peak calling class: '$class
 fi
  let counter=$counter+1
done

# downsample function
function SubSample {

## Calculate the sampling factor based on the intended number of reads:
FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$2 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

if [[ $FACTOR > 1 ]]
  then
  echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
fi

sambamba view -s $FACTOR -f bam -l 9 -t 24 $1

}
###
uni=$out/04_unique
peak=$out/06_peak
deep=$out/07_deeptools
###
echo "*****************************************************************************"
echo "1. downsample same reads"

uni_log=$uni/logs/${sample}.log
uni_out=$uni/$sample
mkdir -p $uni_out
mbam=$uni_out/${sample}_merge.bam
outbam=$uni_out/${sample}_extract.bam
#source activate TensorFlow
bam1=$uni/*${sample}-1/*no_chrYM.sort.rmdup.bam
bam2=$uni/*${sample}-2/*no_chrYM.sort.rmdup.bam
bam3=$uni/*${sample}-3/*no_chrYM.sort.rmdup.bam
# source activate /share/home/zhengzikai/anaconda3/envs/TensorFlow
samtools merge -@ 24 -f $uni_out/${sample}_merge.bam $bam1 $bam2 $bam3
samtools index -@ 24 $uni_out/${sample}_merge.bam
samtools flagstat $uni_out/${sample}_merge.bam > $uni_out/${sample}_merge.flagstat
echo Merge: $sample is Done

nn=`cat $uni/*/*merge.flagstat | grep "+ 0 mapped (100.00% : N/A)" | awk -F" " 'BEGIN {min = 1000000000} {if ($1+0 < min+0) min=$1} END {print min}'`
n=`expr $nn / 1000000 \* 1000000`

echo "Min reads = $n"
SubSample $mbam $n > $outbam
samtools index $outbam
sambamba flagstat -t 24 $outbam

echo "*****************************************************************************"
echo "2. bam2bw using deeptools"

for bin in {10,100,2000}
do

deep_log=$deep/logs/${sample}.bin${bin}.log
deep_out=$deep/$sample/bin$bin
mkdir -p $deep_out

bw=$deep_out/${sample}_bin${bin}.bw
bg=$deep_out/${sample}_bin${bin}.bedGraph

bamCoverage -p $CPU -b $outbam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrY 2>$deep_log 

done

echo "*****************************************************************************"
echo "3. Peak calling using MACS2 or SICER"


mkdir -p $peak/$sample

echo "Peak calling using MACS2 or SICER"

if (( $class == 0  )) ;then

echo 'The file is Input file, need to be control'

elif (( $class == 1  )) ;then

  echo 'Call peak using Macs2; For sharp signals'
  peak_log=$peak/logs/${sample}.macs2.log
  peak_out=$peak/$sample
  
  peaks_xls=$peak_out/${sample}_peaks.xls
  peaks_bed=$peak_out/${sample}.macs2.peaks.bed
  #-c $input
  macs2 callpeak -t $outbam -f BAM -g hs -B -q 0.05 -n ${sample} --outdir $peak_out 2>$peak_log
  #grep '^chr\S' $peaks_xls | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$8"\t""+""\t"$7}' >$peaks_bed
  grep '^[1-9XYM]' $peaks_xls | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$8"\t""+""\t"$7}' >$peaks_bed

elif (( $class == 2  )) ;then
  
  echo 'Call peak using SICER; For sharp signals with  default parameter'
  peak_out=$peak/$sample
  sw=/share/home/jiqianzhao/anaconda3/envs/macs2/bin/sicer
  
  bedt=$peak_out/${sample}.bed
  bedc=$peak_out/${sample}.Input.bed

  bedtools bamtobed -i $outbam > $bedt
  bedtools bamtobed -i $input > $bedc

  sicer -t $bedt -c $bedc -s hg19 -w 200 -fdr 0.01 -o $peak_out -cpu $CPU

elif (( $class == 3  )) ;then
  
  echo 'Call peak using SICER; For broad histone markers like H3K27me3'
  peak_out=$peak/$sample
  
  bedt=$peak_out/${sample}.bed
  bedc=$peak_out/${sample}.Input.bed

  # bedtools bamtobed -i $outbam > $bedt
  # bedtools bamtobed -i $input > $bedc

  sicer -t $bedt -c $bedc -s hg19 -w 500 -fdr 0.01 -g 2500 -o $peak_out -cpu $CPU  


elif (( $class == 4  )) ;then

  echo 'Call peak using SICER; For Super broad histone markers like H3K9me3'
  peak_out=$peak/$sample
  
  bedt=$peak_out/${sample}.bed
  bedc=$peak_out/${sample}.Input.bed

  # bedtools bamtobed -i $outbam > $bedt
  # bedtools bamtobed -i $input > $bedc

  sicer -t $bedt -c $bedc -s hg19 -w 1000 -fdr 0.01 -g 10000 -o $peak_out -cpu $CPU
fi
