##################################################################################################
###                                                                                            ###
###        THIS FILE CONTAINS EXAMPLE CODE FOR FIGURE 7 A&B AND SUPPLEMENTARY FIGURE 7 A&B     ###
###                                                                                            ###
##################################################################################################



rm(list=ls())
setwd("PATH TO WORKING DIRECTORY")

library(GenomicInteractions)
library(InteractionSet)
library(GenomicRanges)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(pheatmap)
library(ggpubr)
library(gginnards)
library(biomaRt)
library(GenomicFeatures)

###===============================================================================================
### 1. LOAD TAD FILE
###===============================================================================================

CELLTYPE.TAD = makeGenomicInteractionsFromFile("CELLTYPE.TAD.5K.bedpe", 
                                type="bedpe", 
                                experiment_name="CELLTYPE.TAD.5kb", 
                                description="CELLTYPE.TAD.5kb")
								
###===============================================================================================
### 2. GET THE LIST OF REFSEQ TRANSCRIPTS
###===============================================================================================								


###===============================================================================================
### 2. 1 MOUSE REFSEQ TRANSCRIPTS
###===============================================================================================	


mm9.refseq.db <- makeTxDbFromUCSC(genome="mm9", tablename = "refGene")
refseq.genes = genes(mm9.refseq.db)
refseq.transcripts = transcriptsBy(mm9.refseq.db, by="gene")
non_pseudogene = names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) 
refseq.transcripts = refseq.transcripts[non_pseudogene] 

###===============================================================================================
### 2. 2 HUMAN REFSEQ TRANSCRIPTS
###===============================================================================================
hg19.refseq.db <- makeTxDbFromUCSC(genome="hg19", table="refGene")
refseq.genes = genes(hg19.refseq.db)
refseq.transcripts = transcriptsBy(hg19.refseq.db, by="gene")
non_pseudogene = names(refseq.transcripts) %in% unlist(refseq.genes$gene_id) 
refseq.transcripts = refseq.transcripts[non_pseudogene] 


###===============================================================================================
### 3. Annotate promoter and terminator
###===============================================================================================

refseq.promoters = promoters(refseq.transcripts, upstream=1000, downstream=1000)
# unlist object so "strand" is one vector
refseq.transcripts.ul = unlist(refseq.transcripts) 
# terminators can be called as promoters with the strand reversed
strand(refseq.transcripts.ul) = ifelse(strand(refseq.transcripts.ul) == "+", "-", "+") 
refseq.terminators.ul = promoters(refseq.transcripts.ul, upstream=5000, downstream=5000) 
# change back to original strand
strand(refseq.terminators.ul) = ifelse(strand(refseq.terminators.ul) == "+", "-", "+") 
# `relist' maintains the original names and structure of the list
refseq.terminators = relist(refseq.terminators.ul, refseq.transcripts)


###===============================================================================================
### 4. Generate geneIDfile
###===============================================================================================

###===============================================================================================
### 4.1 MOUSE---Generate geneIDfile
###===============================================================================================


mm9.ensembl<-  useMart("ensembl", dataset="mmusculus_gene_ensembl")

mm9.values<- as.data.frame(refseq.transcripts)$tx_name

mm9.ensembl_IDs = getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = mm9.values, mart= mm9.ensembl)




mm9.refseq.transcript_id_name = data.frame(as.data.frame(refseq.transcripts)$tx_name, as.data.frame(refseq.transcripts)$group_name)

colnames(mm9.refseq.transcript_id_name) =c("refseq_mrna","gene.id")
head(mm9.refseq.transcript_id_name)

mm9.ensembl_IDs_groupname = merge(mm9.ensembl_IDs,mm9.refseq.transcript_id_name, by.x = "refseq_mrna" )

head(mm9.ensembl_IDs_groupname)

write.table(mm9.ensembl_IDs_groupname, "mm9_refseq_ensembl_hgnc_geneID.txt",
            quote = F,
            row.names = F,
            sep = "\t")
			
###===============================================================================================
### 4.2 HUMAN---Generate geneIDfile
###===============================================================================================


hg19.ensembl<-  useMart("ensembl", dataset="hsapiens_gene_ensembl")

hg19.values<- as.data.frame(refseq.transcripts)$tx_name


hg19.ensembl_IDs = getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), filters = "refseq_mrna", values = hg19.values, mart= hg19.ensembl)




hg19.refseq.transcript_id_name = data.frame(as.data.frame(refseq.transcripts)$tx_name, as.data.frame(refseq.transcripts)$group_name)

colnames(hg19.refseq.transcript_id_name) =c("refseq_mrna","gene.id")
head(hg19.refseq.transcript_id_name)

hg19.ensembl_IDs_groupname = merge(hg19.ensembl_IDs,hg19.refseq.transcript_id_name, by.x = "refseq_mrna" )

head(hg19.ensembl_IDs_groupname)

write.table(hg19.ensembl_IDs_groupname, "hg19_refseq_ensembl_hgnc_geneID.txt",
            quote = F,
            row.names = F,
            sep = "\t")


###===============================================================================================
### 5. ANNOTATE FEATURES
###===============================================================================================

annotation.features = list(promoter=refseq.promoters, 
                           gene.body=refseq.transcripts, 
                           terminator=refseq.terminators
                           
                           )
annotateInteractions(CELLTYPE.TAD, annotation.features)



###====================================================================================================================================
### 6. Bar plot shows the distribution of different types of interaction observed at TAD boundary. Fig.7A & Suppl.Fig.7A
###====================================================================================================================================

ggbarplot(categoriseInteractions(CELLTYPE.TAD), x = "category", y = "count",
          color = "#E7B800",            # Set bar border colors to white
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,
          title = "CELLTYPE.LX", # Rotate vertically x axis texts
          size = 2) +
  rremove("legend")+
  rremove("xlab")+
  font("ylab", size = 30)+
  font("title", size = 25)+
  font("xy.text", size = 25)

ggsave("CELLTYPE.TAD.category.barplot.5kb.png", width = 25, height = 20, units = "cm",dpi=300)


###===============================================================================================
### 7. load gene quantification file
###===============================================================================================

CELLTYPE_RNAseq_GeneQuant = read.csv("CELLTYPE_RNAseq_genequant.txt", h=T, sep = "\t")

head(CELLTYPE_RNAseq_GeneQuant)

###===============================================================================================
### 8. function---"TAD.boundary.gene.exp"
###===============================================================================================

TAD.boundary.gene.exp <- function(anchor.type.1,anchor.type.2,GIObject,genequant,geneIDfile){
   
   anchor.type.1.anchor.type.2.loop <- as.data.frame(GIObject[isInteractionType(GIObject, anchor.type.1,anchor.type.2)])
  
  anchor.type.1.1 = filter(anchor.type.1.anchor.type.2.loop, anchor.type.1.anchor.type.2.loop$node.class1 == anchor.type.1)


anchor.type.1.2 = filter(anchor.type.1.anchor.type.2.loop, anchor.type.1.anchor.type.2.loop$node.class2 == anchor.type.1)


col_name_1 <- paste0(anchor.type.1, ".id1")
  col_name_2 <- paste0(anchor.type.1, ".id2")

   anchor.type.1.gene <- c(unlist(get(col_name_1, anchor.type.1.1)), unlist(get(col_name_2, anchor.type.1.2)))

   
   anchor.type.1.gene.ensembl <- geneIDfile$ensembl_gene_id[match(anchor.type.1.gene, geneIDfile$gene.id)]

anchor.type.1.gene.express = genequant$TPM[match(anchor.type.1.gene.ensembl, genequant$gene_id)]

gene_file <- paste0(anchor.type.1,anchor.type.2, ".gene.txt")
express_file <- paste0(anchor.type.1,anchor.type.2, ".gene.express.txt")


# Combine gene names and expression values
   anchor.type.1.gene_data <- data.frame(Gene = anchor.type.1.gene.ensembl, Expression = anchor.type.1.gene.express)

   # Write the combined data to file
   write.table(anchor.type.1.gene_data, file = gene_file, sep = "\t", quote = FALSE, row.names = FALSE)




#
anchor.type.2.1 = filter(anchor.type.1.anchor.type.2.loop, anchor.type.1.anchor.type.2.loop$node.class1 == anchor.type.2)


anchor.type.2.2 = filter(anchor.type.1.anchor.type.2.loop, anchor.type.1.anchor.type.2.loop$node.class2 == anchor.type.2)



col_name_2_1 <- paste0(anchor.type.2, ".id1")
  col_name_2_2 <- paste0(anchor.type.2, ".id2")

   anchor.type.2.gene <- c(unlist(get(col_name_2_1, anchor.type.2.1)), unlist(get(col_name_2_2, anchor.type.2.2)))

   
anchor.type.2.gene.ensembl <- geneIDfile$ensembl_gene_id[match(anchor.type.2.gene, geneIDfile$gene.id)]
anchor.type.2.gene.express = genequant$TPM[match(anchor.type.2.gene.ensembl, genequant$gene_id)]
   
   

gene_file.2 <- paste0(anchor.type.2, anchor.type.1,".gene.txt")
express_file.2 <- paste0(anchor.type.2, anchor.type.1,".gene.express.txt")

# Combine gene names and expression values
   anchor.type.2.gene_data <- data.frame(Gene = anchor.type.2.gene.ensembl, Expression = anchor.type.2.gene.express)

   # Write the combined data to file
   write.table(anchor.type.2.gene_data, file = gene_file.2, sep = "\t", quote = FALSE, row.names = FALSE)
   






return(list(anchor.type.1.gene, anchor.type.1.gene.express,anchor.type.2.gene, anchor.type.2.gene.express))
}


###===============================================================================================
### 9. function---"TAD.boundary.sametype.gene.exp"
###===============================================================================================

TAD.boundary.sametype.gene.exp <- function(anchor.type,GIObject,genequant,geneIDfile){
   
   anchor.type.loop <- as.data.frame(GIObject[isInteractionType(GIObject, anchor.type,anchor.type)])
  
  col_name_1 <- paste0(anchor.type, ".id1")
  col_name_2 <- paste0(anchor.type, ".id2")

   
    anchor.type.gene.1 <- unlist(get(col_name_1, anchor.type.loop))
    anchor.type.gene.2 <- unlist(get(col_name_2, anchor.type.loop))

   
   anchor.type.gene.1.ensembl <- geneIDfile$ensembl_gene_id[match(anchor.type.gene.1, geneIDfile$gene.id)]

anchor.type.gene.1.express = genequant$TPM[match(anchor.type.gene.1.ensembl, genequant$gene_id)]

gene_file <- paste0(anchor.type,anchor.type, ".gene.1.txt")
express_file <- paste0(anchor.type, anchor.type,".gene.1.express.txt")

# Combine gene names and expression values
   anchor.type.gene.1_data <- data.frame(Gene = anchor.type.gene.1.ensembl, Expression = anchor.type.gene.1.express)

   # Write the combined data to file
   write.table(anchor.type.gene.1_data, file = gene_file, sep = "\t", quote = FALSE, row.names = FALSE)

##############
   #2
   anchor.type.gene.2.ensembl <- geneIDfile$ensembl_gene_id[match(anchor.type.gene.2, geneIDfile$gene.id)]

anchor.type.gene.2.express = genequant$TPM[match(anchor.type.gene.2.ensembl, genequant$gene_id)]

gene_file.2 <- paste0(anchor.type,anchor.type, ".gene.2.txt")
express_file.2 <- paste0(anchor.type, anchor.type,".gene.2.express.txt")

# Combine gene names and expression values
   anchor.type.gene.2_data <- data.frame(Gene = anchor.type.gene.2.ensembl, Expression = anchor.type.gene.2.express)

   # Write the combined data to file
   write.table(anchor.type.gene.2_data, file = gene_file.2, sep = "\t", quote = FALSE, row.names = FALSE)




return(list(anchor.type.gene.1.ensembl, anchor.type.gene.1.express,anchor.type.gene.2.ensembl, anchor.type.gene.2.express))
}


###===============================================================================================
### 10. 
###===============================================================================================


promoter_terminator_pair_gene = TAD.boundary.gene.exp("promoter","terminator",CH12.TAD,CH12_RNAseq_GeneQuant,mm9.ensembl_IDs_groupname)

promoter_genebody_pair_gene =
TAD.boundary.gene.exp("promoter","gene.body",CH12.TAD,CH12_RNAseq_GeneQuant,mm9.ensembl_IDs_groupname)




#
terminator_genebody_pair_gene =
TAD.boundary.gene.exp("terminator","gene.body",CH12.TAD,CH12_RNAseq_GeneQuant,mm9.ensembl_IDs_groupname)

#####

promoter_promoter_pair_gene <- TAD.boundary.sametype.gene.exp("promoter",CH12.TAD,CH12_RNAseq_GeneQuant,mm9.ensembl_IDs_groupname)

terminator_terminator_pair_gene <- TAD.boundary.sametype.gene.exp("terminator",CH12.TAD,CH12_RNAseq_GeneQuant,mm9.ensembl_IDs_groupname)

gene.body_gene.body_pair_gene <- TAD.boundary.sametype.gene.exp("gene.body",CH12.TAD,CH12_RNAseq_GeneQuant,mm9.ensembl_IDs_groupname)



###=================================================================================================================================================================
### 11. Violin plot shows the expression level at each TAD boundary for Promoter-Terminator (PT), Promoter-Genebody (PG) and Promoter-Promoter (PP) boundary pairs.
### Figure 7B & Suppl.Fig.7B
###==================================================================================================================================================================

PTG.type.TAD = c(
             rep("PT.Promoter", length(na.omit(promoter_terminator_pair_gene[[2]]))),
             rep("PT.Terminator", length(na.omit(promoter_terminator_pair_gene[[4]]))),
             rep("PG.Promoter", length(na.omit(promoter_genebody_pair_gene[[2]]))),
             rep("PG.GeneBody", length(na.omit(promoter_genebody_pair_gene[[4]]))),
             rep("Promoter1", length(na.omit(promoter_promoter_pair_gene[[2]]))),
             rep("Promoter2", length(na.omit(promoter_promoter_pair_gene[[4]])))
             )

PTG.value.TAD = c(
              na.omit(promoter_terminator_pair_gene[[2]]),
              na.omit(promoter_terminator_pair_gene[[4]]),
              na.omit(promoter_genebody_pair_gene[[2]]),
              na.omit(promoter_genebody_pair_gene[[4]]),
              na.omit(promoter_promoter_pair_gene[[2]]),
              na.omit(promoter_promoter_pair_gene[[4]])
              )



PTG.value.TAD.2 = log2(PTG.value.TAD+1)

PTG.gene.exp.TAD = data.frame(PTG.type.TAD, PTG.value.TAD.2)

write.table(PTG.gene.exp.TAD,
            "CELLTYPE.PTG.gene.exp.TAD.txt",
            quote= F,
            sep = "\t",
            row.names = F)

PTG.my_comparisons <- List(
                       c("PT.Promoter", "PT.Terminator"),
                       c("PG.Promoter", "PG.GeneBody"),
                       c("Promoter1", "Promoter2"),
                       
                       c("PG.GeneBody", "Promoter1"),
                       c("PG.GeneBody", "Promoter2"),
                       c("PT.Terminator", "Promoter1"),
                       c("PT.Terminator", "Promoter2")
)

ggviolin(PTG.gene.exp.TAD, x = "PTG.type.TAD", y = "PTG.value.TAD.2",
         
          color = "PTG.type.TAD", palette = "jco",
         add = "mean",
         add.params = list(size = 1.5),
         ylab = "log2(TPM+1)",title = "CELLTYPE")+
  stat_compare_means(comparisons = PTG.my_comparisons, label = "p.format", hide.ns = F, size = 5.5, vjust =0.2)+ # Add significance levels
  #stat_compare_means(label.y = 7) +
   font("xlab", size = 18)+
   font("ylab", size = 18)+
   font("xy.text", size = 18) +
  font("title", size = 18)+
   rremove("xlab")+
   rremove("legend")+
   rotate_x_text(angle = 45)

ggsave("CELLTYPE.TAD.Promoter_PTG.GeneExp.6groups.violin.pdf",  width = 18, height = 16, units = "cm",dpi=300)
