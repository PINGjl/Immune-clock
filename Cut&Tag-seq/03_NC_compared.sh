for sample in CD8-G4  CD8-R1;
do
macs2 callpeak -t /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/04_unique/$sample/${sample}_extract.bam -c /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/04_unique/00-CD8-NC/00-CD8-NC_no_chrYM.sort.rmdup.bam -f BAM -g hs -B -q 0.05 -n ${sample} --outdir /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/${sample}_noNC 2>/data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/${sample}_noNC/${sample}.macs2.log
done

for sample in CD8-G4  CD8-R1;
do
grep '^[1-9XYM]' /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/${sample}_noNC/${sample}_peaks.xls | awk '{print $1"\t"$2"\t"$3"\t"$10"\t"$8"\t""+""\t"$7}' > /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/${sample}_noNC/${sample}.macs2.peaks.bed
done

###比较大的region的
computeMatrix scale-regions -p 15 \
 -S /share/home/zhengzikai/03.result/lamin_EMD/${sample}/MSC/07_deeptools/${sample}/bin10/${sample}_bin10.bw /share/home/zhengzikai/03.result/lamin_EMD/${sample}/MSC_WT/07_deeptools/${sample}/bin10/${sample}_bin10.bw \
 -R /share/home/zhengzikai/03.result/lamin_EMD/DamID/bamCompare/iladswitch.bed \
 --beforeRegionStartLength 3000 \
 --regionBodyLength 5000 \
 --afterRegionStartLength 3000 \
 --skipZeros --missingDataAsZero -o /share/home/zhengzikai/03.result/lamin_EMD/${sample}/iladswitch.tab.gz

plotProfile -m /share/home/zhengzikai/03.result/lamin_EMD/DamID/bamCompare/union.HMMT.LAD.tab.gz --colors blue grey --perGroup \
              -out /share/home/zhengzikai/03.result/lamin_EMD/DamID/bamCompare/union.HMMT.LAD.PDF \
              --plotTitle "Test data profile"

###peak
computeMatrix reference-point -p 15 -S /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/07_deeptools/CD8-R1/bin10/CD8-R1_bin10.bw /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/07_deeptools/CD8-G4/bin10/CD8-G4_bin10.bw -R /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/result/Up.diffpeak.macs2.peaks.bed  -a 10000 -b 10000 --referencePoint center -o /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/result/union.LAD.tab.gz

plotProfile -m /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/result/union.LAD.tab.gz --colors blue grey --perGroup -out /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/result/new_LAD.pdf --plotTitle "Test data profile"

plotHeatmap -m /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/result/union.LAD.tab.gz -out /data/pingjiale/05_Result/01_qz_pbmc/03_CUT_TAG_Seq/06_peak/result/union.heat.iLAD.pdf --heatmapHeight 12 --colorList '#1a488e',"#97b2de","white","#de643f","#a73728"  --colorNumber 100 --whatToShow 'heatmap and colorbar'
