#!/bin/bash

omic=CUT_Tag
tar=/data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq
raw=/data/pingjiale/04_Raw_data/01_qz_pbmc/02_CUT_TAG_Seq/ANNO_XS01KF2022100242_PM-XS01KF2022100242-75/Rawdata
index=/data/pingjiale/03_Database/04_bowtie2_index/GRCh37/GRCh37
pip=/data/pingjiale/01_Script/01_qz_pbmc/03_CUT_TAG_Seq/cut_run.sh
each=$tar/scripts/each
mkdir -p $each
#echo $each
for sample in `ls $raw`
do

echo -e '#!/bin/bash' > $each/${sample}_${omic}.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw\\nomic=$omic\\nindex=$index\\npip=$pip\\neach=$each'

fq1=$raw/$sample/${sample}_R1.fq.gz
fq2=$raw/$sample/${sample}_R2.fq.gz

$pip -F $fq1 -f $fq2 -n $sample -o $tar -p 24 -g $index

'>> $each/${sample}_${omic}.sh
done

echo '#!/bin/bash
scripts=./each
cd $scripts
for i in `ls *.sh`
do
sh $i &
done'>$tar/scripts/run.sh
