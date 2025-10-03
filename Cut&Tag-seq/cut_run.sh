#!/bin/bash
#ChIPseq.sh
#This is a pipeline for Cut-seq analysis
#Requirement: conda activate daily
#Current version 1: J.Q.Z, 2021-10-27

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*                     PipeRiboseq: pipeline for Cut-seq analysis.                         *"
  echo "*                            Version 2, 2021-10-27, Q.Z                                    *"
  echo "* Usage: `basename $0`                                                                     *"
  echo "*        Required (paired end data):                                                       *"
  echo "*                  -F [fastq files1]                                                       *"
  echo "*                  -f [fastq files2]                                                       *"
  echo "*                  -n [The prefix samplename]                                              *"
  echo "*                  -o [The outputpath]                                                     *"
  echo "*        Optional: -p [Number of CPUs, default=24]                                         *"
  echo "*                  -g [The file of index files]                                            *"
  echo "* Inputs: The raw fastq files                                                              *"
  echo "* Run: Default to run trim_galore,fastqc,bowtie2,sam2bam,rmdup & deeptools                 *"
  echo "*      Figures will be generated in /plots folder, and bigWig files in /tracks folder      *"
  echo "* Outputs: All output files will be generated in the same folder as the pipeline submitted *"
  echo "* This pipeline requires 'conda activate daily'                                            *"
  echo "********************************************************************************************"
  exit 1
fi

#Get parameters
fq1="unassigned"
fq2="unassigned"
sample="unassigned"
out="unassigned"
CPU=24
index="unassigned"

for arg in "$@"
do
 if [[ $arg == "-F" ]]
  then
    fq1=${array[$counter+1]}
    echo '   raw fastq1: '$fq1
 elif [[ $arg == "-f" ]]
  then
    fq2=${array[$counter+1]}
    echo '   raw fastq2: '$fq2
 elif [[ $arg == "-n" ]]
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
 elif [[ $arg == "-g" ]]
  then
    index=${array[$counter+1]}
    echo '   Index: '$index
 fi
  let counter=$counter+1
done

echo "*****************************************************************************"
echo "1. Generate the all output files"

cd $out
mkdir -p 01_fastqc 02_trim 03_mapping 04_unique 06_peak 07_deeptools 08_plot

trim=$out/02_trim
map=$out/03_mapping
uni=$out/04_unique
bb=$out/05_bam2bed
deep=$out/07_deeptools
peak=$out/06_peak

for path in `ls $out | grep 0`
do
mkdir -p  $path/scripts
mkdir -p $path/logs
done

echo "*****************************************************************************"
echo "2. Run trim_galore"

trim_log=$trim/logs/${sample}.log
trim_out=$trim/$sample
mkdir -p $trim_out
trim_galore --fastqc --path_to_cutadapt /data/pingjiale/02_Software/01_anaconda/anaconda3/envs/trim_galore/bin/cutadapt --stringency 3 --paired --output_dir $trim_out $fq1 $fq2 2>$trim_log

echo "*****************************************************************************"
echo "3. Mapping using bowtie2"

map_log=$map/logs/${sample}.log
map_out=$map/$sample
mkdir -p $map_out

clean1=$trim_out/${sample}*1.fq.gz
clean2=$trim_out/${sample}*2.fq.gz

bowtie2 -p 24 -x $index --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -t -1 $clean1 -2 $clean2 -S $map_out/${sample}.sam 2>$map_log

echo "*****************************************************************************"
echo "4. Sam2bam and rmdup"

uni_log=$uni/logs/${sample}.log
uni_out=$uni/$sample
mkdir -p $uni_out

sam=$map_out/${sample}.sam
bam=$uni_out/${sample}_no_chrYM.bam
srt_bam=$uni_out/${sample}_no_chrYM.sort.bam
srt_rmdup_bam=$uni_out/${sample}_no_chrYM.sort.rmdup.bam
METRICS_FILE=$uni_out/${sample}.picard.matrix
idxstats_txt=$uni_out/${sample}.sort.rmdup.idxstats.txt
flagstat_txt=$uni_out/${sample}.sort.rmdup.flagstat.txt

#awk '$3!="chrM" && $3!="chrY"' $sam | samtools view -@ $CPU -S -b -q 10 > $bam 2>$uni_log &&
awk '$3!="M" && $3!="Y"' $sam | samtools view -@ $CPU -S -b -q 10 > $bam 2>$uni_log &&
samtools sort -@ $CPU -l 9 $bam -o $srt_bam 2>>$uni_log &&
picard MarkDuplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true INPUT=$srt_bam METRICS_FILE=$METRICS_FILE OUTPUT=$srt_rmdup_bam 2>>$uni_log &&
samtools index -@ $CPU $srt_rmdup_bam 2>>$uni_log &&
samtools flagstat -@ $CPU $srt_rmdup_bam > $flagstat_txt 2>>$uni_log

echo "*****************************************************************************"
echo "5. bam2bw using deeptools"

for bin in {10,100,2000,25000}
do

deep_log=$deep/logs/${sample}.bin${bin}.log
deep_out=$deep/$sample/bin$bin
mkdir -p $deep_out

bw=$deep_out/${sample}_bin${bin}.bw
bg=$deep_out/${sample}_bin${bin}.bedGraph

bamCoverage -p 24 -b $srt_rmdup_bam -o $bw --normalizeUsing RPKM --binSize $bin --ignoreForNormalization chrY 2>$deep_log

done


echo "Peak calling using MACS2 or SICER"
mkdir -p $peak/$sample
outbam=$srt_rmdup_bam
class=1
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
  cp $peak_out/${sample}-W200-G600-FDR0.01-island.bed $peak_out/${sample}_peak.bed

elif (( $class == 3  )) ;then
  
  echo 'Call peak using SICER; For broad histone markers like H3K27me3'
  peak_out=$peak/$sample
  
  bedt=$peak_out/${sample}.bed
  bedc=$peak_out/${sample}.Input.bed

  bedtools bamtobed -i $outbam > $bedt
  bedtools bamtobed -i $input > $bedc

  sicer -t $bedt -c $bedc -s hg19 -w 500 -fdr 0.01 -g 2500 -o $peak_out -cpu $CPU  
  cp $peak_out/${sample}-W500-G2500-FDR0.01-island.bed $peak_out/${sample}_peak.bed

elif (( $class == 4  )) ;then

  echo 'Call peak using SICER; For Super broad histone markers like H3K9me3'
  peak_out=$peak/$sample
  
  bedt=$peak_out/${sample}.bed
  bedc=$peak_out/${sample}.Input.bed

  bedtools bamtobed -i $outbam > $bedt
  bedtools bamtobed -i $input > $bedc

  sicer -t $bedt -c $bedc -s hg19 -w 1000 -fdr 0.01 -g 10000 -o $peak_out -cpu $CPU
  cp $peak_out/${sample}-W1000-G10000-FDR0.01-island.bed $peak_out/${sample}_peak.bed
fi
