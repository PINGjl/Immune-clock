#!/bin/bash

tar=/data/pingjiale/05_Result/01_qz_pbmc/01_bulkRNA
raw=/data/pingjiale/04_Raw_data/01_qz_pbmc/01_bulk/X101SC24114380-Z01-J003/01.RawData/
map=$tar/03_map

for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$map/scripts/${sample}_map.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'
echo Mapping:$sample is Starting
trim=$tar/02_trim
map=$tar/03_map

clean1=$trim/$sample/*_val_1.fq.gz
clean2=$trim/$sample/*_val_2.fq.gz

index=/data/pingjiale/03_Database/01_hisat2_index/hg19/genome

result=$map/$sample
log=$map/logs
mkdir $result

hisat2=/data/pingjiale/02_Software/01_anaconda/anaconda3/envs/trim_galore/bin/hisat2
$hisat2 -p 24 -x $index -1 $clean1 -2 $clean2 --dta -S $result/${sample}.sam 2>$log/${sample}_map.log

echo Mapping:$sample is Done

bam=$result/${sample}.bam
sort=$result/${sample}.sort.bam
sort1=$result/${sample}.sort.name.bam

samtools=/data/pingjiale/02_Software/01_anaconda/anaconda3/envs/trim_galore/bin/samtools
$samtools view -bS -@ 24 -q 10 $result/${sample}.sam>$bam &&
$samtools sort -@ 24 $bam -o $sort &&
$samtools sort -@ 24 -n $bam -o $sort1 &&
$samtools index -@ 24 $sort


'>>$map/scripts/${sample}_map.sh

done


echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
sh $i &
done'>$map/run_map.sh

