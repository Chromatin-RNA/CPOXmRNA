######################################################################################################################################
###                                                                                                                                ###
###                         THIS FILE CONTAINS CODE FOR EXTENDED FIGURE. 4A,C,D,F AND G                                            ###
###                                                                                                                                ###
######################################################################################################################################



###======================================================================================================================================
### 1. DOWNLOAD CHIP-SEQ FILES FROM ENCODE
###======================================================================================================================================

#mm10
# H3K27ac
wget https://www.encodeproject.org/files/ENCFF972YIT/@@download/ENCFF972YIT.bed.gz

wget https://www.encodeproject.org/files/ENCFF847MJA/@@download/ENCFF847MJA.bigWig

# H3K27me3
wget https://www.encodeproject.org/files/ENCFF932OHQ/@@download/ENCFF932OHQ.bed.gz

wget https://www.encodeproject.org/files/ENCFF336CRJ/@@download/ENCFF336CRJ.bigWig

#  EP300
wget https://www.encodeproject.org/files/ENCFF185UAQ/@@download/ENCFF185UAQ.bed.gz

wget https://www.encodeproject.org/files/ENCFF116DCE/@@download/ENCFF116DCE.bigWig


###======================================================================================================================================
### 2. Extended H3K27AC PEAK +/-1kb
###======================================================================================================================================

# h3k27ac 1kb extended peak

wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes

module load bedtools

bedtools slop  -i ENCFF972YIT.bed.gz -g mm10.chrom.sizes -b 1000 > MEL_H3K27ac_1k.bed



###======================================================================================================================================
### 3. BINARY HEATMAP FOR H3K27AC AND P300, GOT P200 ONLY PEAK AND P300/H3K27AC OVERLAPPING PEAK 
###======================================================================================================================================

######################################################### binary.heatmap.py#########################################################
import os
from matplotlib import pyplot as plt
from pybedtools.contrib import plotting
import pybedtools
import numpy as np


# set up the order in which to plot the columns of the binary heatmap
names, bts = zip(*[
    ('H3K27ac', 'MEL_H3K27ac_1k.bed'),
    ('p300', 'ENCFF185UAQ.bed.gz'),
])
bts = [pybedtools.BedTool(i).sort() for i in bts]

# set up the object by giving it a list of pybedtool.BedTool objects and a list
# of names to use.
b = plotting.BinaryHeatmap(
    bts=bts,
    names=names)

# plot it
b.plot()

# write out how many genomic location of each class were identified
with open('class_counts.txt', 'w') as fout:
    for cls, cnt in sorted(b.class_counts.items(), key=lambda x: x[1], reverse=True):
        fout.write('{0:>25}: {1:<15}\n'.format(cls, cnt))

# write out the actual intervals from each class
out_dir = 'intervals'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

for k, v in b.classified_intervals.items():
    label = k.replace(',', '_and_')
    v.cut([0, 1, 2]).saveas(os.path.join(out_dir, label + '.bed'), trackline='track name="%s"' % label)

# save the figure
fig = plt.gcf()
fig.tight_layout()
fig.savefig('binary_heatmap.png')
plt.show()

#########################################################

python binary.heatmap.py



###======================================================================================================================================
### 4. EXTENDED FIGURE. 4A
###======================================================================================================================================



#################
#plot heatmap

module load deeptools

computeMatrix reference-point -R /data/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/Binary.heatmap/intervals/H3K27ac_and_p300.bed /data/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/Binary.heatmap/intervals/p300.bed -S ENCFF116DCE.bigWig ENCFF847MJA.bigWig ENCFF336CRJ.bigWig \
  -a 5000 -b 5000 --skipZeros -o MEL_p300_H3K27ac_twoGroups_scaled2.gz --outFileNameMatrix MEL_p300_H3K27ac_twoGroups_scaled2.tab --outFileSortedRegions MEL_p300_H3K27ac_twoGroups_genes2.bed 

plotHeatmap -m MEL_p300_H3K27ac_twoGroups_scaled2.gz \
     -out MEL.p300.png \
     --colorMap RdBu_r \
	 --refPointLabel "p300" \
     --whatToShow 'heatmap and colorbar' \
     --zMin -3 --zMax 3 
	 
	 
	 
###======================================================================================================================================
### 5. EXTENDED FIGURE. 4C
###======================================================================================================================================

module load homer

findMotifsGenome.pl /data/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/Binary.heatmap/intervals/p300.bed mm10 MEL_p300_only_MotifOutput/ -size 200

###======================================================================================================================================
### 6. EXTENDED FIGURE. 4D
###======================================================================================================================================

#ELF1
wget https://www.encodeproject.org/files/ENCFF453OVN/@@download/ENCFF453OVN.bed.gz

wget https://www.encodeproject.org/files/ENCFF688GOC/@@download/ENCFF688GOC.bigWig

#GATA1
wget https://www.encodeproject.org/files/ENCFF841DLH/@@download/ENCFF841DLH.bed.gz

wget https://www.encodeproject.org/files/ENCFF222HZM/@@download/ENCFF222HZM.bigWig

#MYC
wget https://www.encodeproject.org/files/ENCFF152JNC/@@download/ENCFF152JNC.bed.gz

wget https://www.encodeproject.org/files/ENCFF944GOL/@@download/ENCFF944GOL.bigWig

module load deeptools

computeMatrix reference-point -R /data/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/Binary.heatmap/intervals/p300.bed -S ENCFF688GOC.bigWig ENCFF222HZM.bigWig ENCFF944GOL.bigWig \
  -a 5000 -b 5000 --skipZeros -o MEL_p300_only_TF_scaled.gz --outFileNameMatrix MEL_p300_only_TF_scaled.tab --outFileSortedRegions MEL_p300_only_TF.bed 
  
  
plotProfile -m MEL_p300_only_TF_scaled.gz -out MEL_p300_only_TF.png --perGroup --refPointLabel "p300 only" --samplesLabel "ELF1" "GATA1" "MYC" --plotTitle "MEL.p300.only.peak.overlapping.TFs.profile.2" --colors red blue green

###======================================================================================================================================
### 7. LIFTOVER LOOP FILES FROM MM9 TO MM10
###======================================================================================================================================

# APA analysis
#convert "p300 only all loops" file and "p300 only-H3K27me3 loops" file from mm9 to mm10

#P300.only.all loop
awk -v OFS='\t' '{print "chr"$1,$2,$3,"chr"$4,$5,$6}' /data/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/MEL.p300.ONLY.ALL.loop.txt > MEL.p300.ONLY.ALL.loop.mm9.bedpe 

# excel add loop name, score, strand1, strand2 to MEL.p300.ONLY.ALL.loop.mm9.bedpe, so that it becomes p300_only_loop.mm9.bedpe
# convert 
#https://github.com/cauyrd/liftOverBedpe/blob/main/liftOverBedpe.py

python bedpe/liftOverBedpe/liftOverBedpe.py --lift bedpe/liftOverBedpe/liftOver --chain mm9ToMm10.over.chain.gz --i p300_only_loop.mm9.bedpe --o MEL.p300.ONLY.ALL.loop.mm10.bedpe


#P300.only.H3K27me3 loop
awk -v OFS='\t' '{print "chr"$2,$3,$4,"chr"$5,$6,$7}' /data/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/MEL.p300.H3K27me3.loop.mm9.txt > MEL.p300.H3K27me3.loop.mm9.bedpe 
python bedpe/liftOverBedpe/liftOverBedpe.py --lift bedpe/liftOverBedpe/liftOver --chain mm9ToMm10.over.chain.gz --i p300.H3K27me3.loop.mm9.bedpe  --o MEL.p300.H3K27me3.loop.mm10.bedpe



###======================================================================================================================================
### 8. EXTENDED FIGURE. 4F
###======================================================================================================================================

module load juicer

juicer_tools apa  --threads 1 -w 6 -k VC_SQRT -u /data/xieb2/G1ER/4DNFI6H926RO.hic /gpfs/gsfs8/users/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/p300.only.loop.mm10.txt MEL.p300all.LOOP.APA_all_VC_SQRT_6

###======================================================================================================================================
### 9. EXTENDED FIGURE. 4G
###======================================================================================================================================

juicer_tools apa  --threads 1 -w 6 -k VC_SQRT -u /data/xieb2/G1ER/4DNFI6H926RO.hic /gpfs/gsfs8/users/xieb2/2023_MEL_p300_PRC2_loop/20230523.MEL.p300/MEL.p300.H3K27me3.loop.mm10.txt MEL.p300.H3K27me3.LOOP.APA_all_VC_SQRT_6



