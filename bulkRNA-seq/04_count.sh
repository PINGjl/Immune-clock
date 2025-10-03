#!/bin/bash
tar=/data/pingjiale/05_Result/01_qz_pbmc/01_bulkRNA
raw=/data/pingjiale/04_Raw_data/01_qz_pbmc/01_bulk/X101SC24114380-Z01-J003/01.RawData/
count=$tar/04_count

for sample in `ls $raw`
do

echo -e '#!/bin/bash'>$count/scripts/${sample}_count.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

gtf=/data/pingjiale/03_Database/02_ref_gtf/Homo_sapiens.GRCh37.87.final.gtf

uni=$tar/03_map
srt_bam=$uni/$sample/${sample}.sort.bam
count=$tar/04_count
result=$count/$sample
log=$count/logs
mkdir $result

echo HTseq:$sample is Starting

htseq=/data/pingjiale/02_Software/01_anaconda/anaconda3/envs/trim_galore/bin/htseq-count
$htseq -f bam -r name -s no -a 10 $srt_bam $gtf > $result/${sample}_all.txt 2>$log/${sample}_count.log

echo HTseq-count has been Done'>>$count/scripts/${sample}_count.sh

done


echo '#!/bin/bash
scripts=./scripts
cd $scripts
for i in `ls *.sh`
do
sh $i &
done'>$count/run_count.sh

