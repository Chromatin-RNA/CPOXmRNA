###############################################################################################
###                                                                                         ###
###       THIS FILE CONTAINS CODE FOR Extended Data Fig. 4 B,E,H                            ###
###                                                                                         ###
###############################################################################################


library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicInteractions)
library(GenomicFeatures)
library(ggpubr)




###================================================================================================
### 1. LOAD P300 ONLY PEAK (mm10)
###================================================================================================

MEL.P300.only.peak <- readPeakFile("MEL.p300.only.bed")
MEL.P300.only.peak



###================================================================================================
### 2. ANNOTATE P300 ONLY PEAKS WITH CHIPSEEKER
###================================================================================================

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakAnno <- annotatePeak(MEL.P300.only.peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")


###================================================================================================
### 3. Extended Data Fig. 4B
###================================================================================================

plotAnnoPie(peakAnno)

###================================================================================================
### 4. LOAD TAD AND LOOP FILES
###================================================================================================

# load G1ER merged TAD (mm9)
G1ER.TAD = makeGenomicInteractionsFromFile("G1ER_merged_TAD.bedgraph", 
                                type="bedpe", 
                                experiment_name="G1ER.TAD.5kb", 
                                description="G1ER.TAD.5kb")



#
# load G1ER merged loop (mm9)

G1ER.loop = makeGenomicInteractionsFromFile("G1ER_merged_loop.bedgraph", 
                                type="bedpe", 
                                experiment_name="G1ER.loop", 
                                description="G1ER.loop")

###================================================================================================
### 5. LOAD MM9 VERSION OF P300 ONLY PEAK FILE AND H3K27ME3 FILE
###================================================================================================

# p300 only peak mm9
MEL.P300.only.peak.mm9 <- readPeakFile("MEL.P300.only.mm9.bed")

# H3K27me3 peak mm9

MEL.H3K27me3.mm9 = readPeakFile("MEL.H3K27me3.mm9.bed")


###================================================================================================
### 6. LOAD mm9 GENE ANNOTATION
###================================================================================================

mm9.ens.db <- makeTxDbFromUCSC(genome="mm9", table="ensGene")
ens.genes = genes(mm9.ens.db)

###================================================================================================
### 7. ANNOTATE TAD
###================================================================================================

# annotate TAD
annotation.features.TAD = list(p300.only.all=MEL.P300.only.peak.mm9,
                               H3K27me3 = MEL.H3K27me3.mm9,
                           gene = ens.genes)
annotateInteractions(G1ER.TAD, annotation.features.TAD)
categoriseInteractions(G1ER.TAD)

# percentage
G1ER.TAD.fraction <- categoriseInteractions(G1ER.TAD)$count/sum(categoriseInteractions(G1ER.TAD)$count)*100
G1ER.TAD.fraction


 
###================================================================================================
### 8. LOAD ENCODE MEL RNA-SEQ GENE QUANTIFICATION FILE
###================================================================================================

# load gene quantification file

MEL_RNAseq_GeneQuant = read.csv("MEL_gene_quant.txt", h=T, sep = "\t")

head(MEL_RNAseq_GeneQuant)


###================================================================================================
### 9. FUNCTION---p300-H3K27me3 loop
###================================================================================================

# p300-H3K27me3 loop, expression level function
p300.only.H3K27me3.loop = function(x){
  x.1 = filter(x, x$node.class1 == "H3K27me3")


x.2 = filter(x, x$node.class2 == "H3K27me3")



x.gene = c(unlist(x.1$gene.id1), unlist(x.2$gene.id2))


x.gene.express = MEL_RNAseq_GeneQuant$TPM[match(x.gene, MEL_RNAseq_GeneQuant$gene_id)]
return(x.gene.express)
}

###================================================================================================
### 10. FUNCTION---ONLY.GENE
###================================================================================================

#only gene function

only.gene.exp = function(x){
  x.1 = filter(x, x$node.class1 == "gene")


x.2 = filter(x, x$node.class2 == "gene")



x.gene = c(unlist(x.1$gene.id1), unlist(x.2$gene.id2))


x.gene.express = MEL_RNAseq_GeneQuant$TPM[match(x.gene, MEL_RNAseq_GeneQuant$gene_id)]
return(x.gene.express)
}

###=====================================================================================================
### 11. TAD --- EXPRESSION VALUE OF P300 ONLY INTERACTING GENES (H3K27ME3 GENES VS. NON-H3K27ME3 GENES)
###=====================================================================================================


# H3K27me3-P300 TAD 

MEL_p300.only.H3K27me3.TAD = as.data.frame(G1ER.TAD[isInteractionType(G1ER.TAD, "H3K27me3", "p300.only.all")])


MEL_p300.only.H3K27me3.gene.express = p300.only.H3K27me3.loop(MEL_p300.only.H3K27me3.TAD)



# P300 only gene --TAD

MEL_p300.only = as.data.frame(G1ER.TAD[isInteractionType(G1ER.TAD, "p300.only.all", "gene")])

MEL_p300.only.exp = only.gene.exp(MEL_p300.only)

###================================================================================================
### 12. Extended Data Fig. 4H---LEFT PANEL
###================================================================================================
type.TAD.p300 = c(rep("non-H3K27me3 gene", length(MEL_p300.only.exp)),

rep("H3K27me3 gene", length(MEL_p300.only.H3K27me3.gene.express)))

value.TAD.p300 = c(MEL_p300.only.exp,MEL_p300.only.H3K27me3.gene.express)



value.TAD.p300.2 = log10(value.TAD.p300)

p300.H3K27me3.gene.exp.TAD.p300 = cbind.data.frame(type.TAD.p300, value.TAD.p300.2)


my_comparisons.p300 <- list( c("non-H3K27me3 gene", "H3K27me3 gene"))


ggviolin(p300.H3K27me3.gene.exp.TAD.p300, x = "type.TAD.p300", y = "value.TAD.p300.2", fill = "type.TAD.p300",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"),xlab = "TAD Boundary",
  ylab = "log10(TPM)")+
  stat_compare_means(comparisons = my_comparisons.p300,label = "p.signif", label.x = 1.5, label.y = 6)+ # Add significance levels
  stat_compare_means(label.y = 8) +
   font("xlab", size = 18)+
 font("ylab", size = 18)+
 font("xy.text", size = 18) 


###================================================================================================
### 13. ANNOTATE LOOPS
###================================================================================================

# annotate Loop
annotation.features.loop = list(p300.only.all=MEL.P300.only.peak.mm9,
                               H3K27me3 = MEL.H3K27me3.mm9,
                           gene = ens.genes)
annotateInteractions(G1ER.loop, annotation.features.loop)
categoriseInteractions(G1ER.loop)

# Percentage

G1ER.loop.fraction <- categoriseInteractions(G1ER.loop)$count/sum(categoriseInteractions(G1ER.loop)$count)*100
G1ER.loop.fraction



###================================================================================================
### 14. LOOP --- EXPRESSION VALUE OF P300 ONLY INTERACTING GENES (H3K27ME3 GENES VS. NON-H3K27ME3 GENES)
###================================================================================================

# H3K27me3-P300 loop 

MEL_p300.only.H3K27me3.loop = as.data.frame(G1ER.loop[isInteractionType(G1ER.loop, "H3K27me3", "p300.only.all")])


p300.only.H3K27me3.loop.exp = p300.only.H3K27me3.loop(MEL_p300.only.H3K27me3.loop)


# P300 only gene --loop

MEL_p300.only.loop = as.data.frame(G1ER.loop[isInteractionType(G1ER.loop, "p300.only.all", "gene")])

MEL_p300.only.loop.exp = only.gene.exp(MEL_p300.only.loop)

###================================================================================================
### 15. Extended Data Fig. 4H---RIGHT PANEL
###================================================================================================

type.loop.p300 = c(rep("non-H3K27me3 gene", length(MEL_p300.only.loop.exp)),

rep("H3K27me3 gene", length(p300.only.H3K27me3.loop.exp)))

value.loop.p300 = c(MEL_p300.only.loop.exp, p300.only.H3K27me3.loop.exp)



value.loop.p300.2 = log10(value.loop.p300)

p300.H3K27me3.gene.exp.loop.p300 = cbind.data.frame(type.loop.p300, value.loop.p300.2)



ggviolin(p300.H3K27me3.gene.exp.loop.p300, x = "type.loop.p300", y = "value.loop.p300.2", fill = "type.loop.p300",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"),xlab = "Chromatin loop",
  ylab = "log10(TPM)")+
  stat_compare_means(comparisons = my_comparisons.p300, label = "p.signif",label.y = 6, label.x = 1.5)+ # Add significance levels
  stat_compare_means(label.y = 8) +
   font("xlab", size = 18)+
 font("ylab", size = 18)+
 font("xy.text", size = 18) 

###================================================================================================
### 16. Extended Data Fig. 4E
###================================================================================================

# mirrored chart to show the interaction freq of TAD and loops
require("ggplot2")
require("gridExtra")

dataToPlot <- data.frame(
  "Category" = c("gene-gene", 
                 "gene-Intergenic", 
                 "gene-p300 only",
                 "gene-H3K27me3",
                 "Intergenic-Intergenic",
                 "Intergenic-p300 only",
                 "Intergenic-H3K27me3",
                 "p300 only-p300 only",
                 "p300 only-H3K27me3",
                 "H3K27me3-H3K27me3"
                 ),
  "TAD" = c(7169,			
	8136,			
	1448,			
	3192,			
	5434,			
	861,		
	1780,			
	102,			
	445,			
	710	),
  "loop" = c(6260,
             1849,
             4751,
             3672,
             468,
             832,
             686,
             1339,
             1993,
             2080))

plot1 <- ggdotchart(dataToPlot, x = "Category", y = "TAD",
           color = "#00AFBB",                                # Color by groups
           
          
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2),
           rotate = TRUE,                                # Rotate vertically
           group = "Category",                                # Order by groups
           dot.size = 6,                                 # Large dot size
          
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           ) +
  scale_y_continuous(trans = "reverse") +
  scale_x_discrete(position = "top") +
  theme(
    axis.text.y = element_blank()
  ) +
  labs(x = NULL)


plot2 <- ggdotchart(dataToPlot, x = "Category", y = "loop",
           color = "#FC4E07",                                # Color by groups
           
          
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2),
           rotate = TRUE,                                # Rotate vertically
           group = "Category",                                # Order by groups
           dot.size = 6,                                 # Large dot size
          
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           ) +
  
  theme(
    axis.text.y = element_text(size = 16, hjust = 0.5, margin=margin(r=30))
  ) +
  labs(x = "") 
  # Add dashed grids


TAD.loop.interaction.categories.plot = gridExtra::grid.arrange(plot1, plot2, ncol = 2, widths = c(1,2.5)) 


