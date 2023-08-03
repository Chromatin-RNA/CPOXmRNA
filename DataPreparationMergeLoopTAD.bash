
###################################################################################################################################################################
###                                                                                                                                                             ###
###                                THIS FILE CONTAINS CODE FOR MERGE LOOP AND TAD FROM CELL CYCLE G1ER HIC DATA                                                 ###
###                                                                                                                                                             ###
###################################################################################################################################################################



###====================================================================================================================================================================
###   MERGE LOOP FILES
###====================================================================================================================================================================



#Download data from GEO database
GSE129997_earlyG1_loops.tsv
GSE129997_lateG1_loops.tsv
GSE129997_midG1_loops.tsv
GSE129997_prometa_loops.tsv

# convert .tsv to .bedgraph
awk 'BEGIN{OFS="\t";} {print $6,$7,$8,$9,$10,$11}' GSE129997_earlyG1_loops.tsv | tail -n +2 > GSE129997_earlyG1_loops.bedgraph
awk 'BEGIN{OFS="\t";} {print $6,$7,$8,$9,$10,$11}' GSE129997_lateG1_loops.tsv | tail -n +2 > GSE129997_lateG1_loops.bedgraph
awk 'BEGIN{OFS="\t";} {print $6,$7,$8,$9,$10,$11}' GSE129997_midG1_loops.tsv | tail -n +2 > GSE129997_midG1_loops.bedgraph
awk 'BEGIN{OFS="\t";} {print $6,$7,$8,$9,$10,$11}' GSE129997_prometa_loops.tsv | tail -n +2 > GSE129997_prometa_loops.bedgraph

# merge cell cycle bedgraph files

module load hicexplorer
 
hicMergeLoops -i  GSE129997_earlyG1_loops.bedgraph GSE129997_lateG1_loops.bedgraph GSE129997_midG1_loops.bedgraph GSE129997_prometa_loops.bedgraph -o G1ER_merged_loop.bedgraph -r 5000



###====================================================================================================================================================================
### CREATE TAD DOMAIN BEDPE FILE WITH +/-5KB OF THE ANCHOR, THEN MERGE
###====================================================================================================================================================================

GSE129997_earlyG1_domains.bed GSE129997_lateG1_domains.bed GSE129997_midG1_domains.bed GSE129997_prometa_domains.bed

awk ‘{OFS=“} {print $1, $2-5000,$2+5000, $1, $3-5000, $3+5000,”G1ER_earlyG1_TAD_“NR}’ GSE129997_earlyG1_domains.bed > GSE129997_earlyG1_domains.bedgraph 
awk ‘{OFS=“} {print $1, $2-5000,$2+5000, $1, $3-5000, $3+5000,”G1ER_lateG1_TAD_“NR}’ GSE129997_lateG1_domains.bed > GSE129997_lateG1_domains.bed.bedgraph 
awk ‘{OFS=“} {print $1, $2-5000,$2+5000, $1, $3-5000, $3+5000,”G1ER_midG1_TAD_“NR}’ GSE129997_midG1_domains.bed > GSE129997_midG1_domains.bedgraph 
awk ‘{OFS=“} {print $1, $2-5000,$2+5000, $1, $3-5000, $3+5000,”G1ER_prometa_TAD_“NR}’ GSE129997_prometa_domains.bed > GSE129997_prometa_domains.bedgraph



hicMergeLoops -i GSE129997_earlyG1_domains.bedgraph GSE129997_lateG1_domains.bed.bedgraph GSE129997_midG1_domains.bedgraph GSE129997_prometa_domains.bedgraph -o G1ER_merged_TAD.bedgraph -r 5000