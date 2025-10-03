#!/bin/bash

tar=/data/pingjiale/05_Result/01_qz_pbmc/01_bulkRNA
raw=/data/pingjiale/04_Raw_data/01_qz_pbmc/01_bulk/X101SC24114380-Z01-J003/01.RawData/
trim=$tar/02_trim

for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$trim/scripts/${sample}_trim.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Trimming:$sample is Starting

trim=$tar/02_trim
log=$trim/logs

fq1=$raw/$sample/*1.fq.gz 
fq2=$raw/$sample/*2.fq.gz

result=$trim/$sample
mkdir $result

trim_galore=/data/pingjiale/02_Software/01_anaconda/anaconda3/envs/trim_galore/bin/trim_galore
$trim_galore --fastqc --path_to_cutadapt /data/pingjiale/02_Software/01_anaconda/anaconda3/envs/trim_galore/bin/cutadapt --stringency 3 --paired --output_dir $result $fq1 $fq2 2>>$log/${sample}.log

echo Trimming has been Done'>>$trim/scripts/${sample}_trim.sh

done

echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
sh $i &
done'>$trim/run_trim.sh
