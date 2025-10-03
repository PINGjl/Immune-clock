#!/bin/bash
omic=CUT_Tag
tar=/data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq
raw=/data/pingjiale/04_Raw_data/01_qz_pbmc/02_CUT_TAG_Seq/ANNO_XS01KF2022100242_PM-XS01KF2022100242-75/Rawdata
index=/data/pingjiale/03_Database/04_bowtie2_index/GRCh37/GRCh37
pip=/data/pingjiale/01_Script/01_qz_pbmc/03_CUT_TAG_Seq/cut_run_merged.sh

for sample in {CD8-G4,CD8-R1}
do
tar=/data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq
#mkdir -p $tar
each=$tar/scripts/each
#mkdir -p $each
out=$tar
echo -e '#!/bin/bash' > $each/${sample}_${omic}.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw\\nomic=$omic\\nindex=$index\\npip=$pip'

fq1=$raw/$sample*1.fq.gz
fq2=$raw/$sample*2.fq.gz

$pip -F $fq1 -f $fq2 -n $sample -o $tar -p 24 -i $index -c 1

'>> $each/${sample}_${omic}.sh

sh $each/${sample}_${omic}.sh
done