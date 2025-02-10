# Full length 16S rRNA gene sequences are imported of the SynCom strains to verify taxonomic assignments of ASVs, 
# present in SynCom samples. 

library(DECIPHER)
library(Biostrings)
fasta_16S<-readDNAStringSet("../Genomes/16S sequences SynCom.fasta")


library(phyloseq)
phyl_1<-readRDS("phyloseq_incl_meta.rds")

phyl_SC_ino<-subset_samples(phyl_1, Treatment %in% "SynCom" & Timepoint %in% "T0")
phyl_SC_ino<-prune_taxa(taxa_sums(phyl_SC_ino) > 0, phyl_SC_ino)
#29 ASVs

ASV_seqs<-DNAStringSet(refseq(phyl_SC_ino))
tax_table(phyl_SC_ino)
#Merge the two DNA string sets
combin_DNA<-c(fasta_16S,ASV_seqs)
#Align
align_DNA<-AlignSeqs(combin_DNA)


library(phangorn)

phang.align <- phyDat(as(align_DNA, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

library(ggtree)
library(tidytree)

ggtree(treeNJ)+geom_tiplab()

tax_table(phyl_SC_ino)

Para_data<-subset_taxa(phyl_1, Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" | Family == "Rhizobiaceae")

para_seqs<-DNAStringSet(refseq(Para_data))
para_align<-AlignSeqs(c(fasta_16S,para_seqs))
para.align <- phyDat(as(para_align, "matrix"), type="DNA")
dm_para <- dist.ml(para.align)
paratreeNJ <- NJ(dm_para)
# plot tree
ggtree(paratreeNJ)+ geom_text(aes(label=node), hjust=-.3)+geom_tiplab()
# a bit too crowded, so I subset based on a node
para_red<-tree_subset((paratreeNJ), node = 346, levels_back = 3)
ggtree(para_red)+geom_tiplab()

#Based on this the ASVs ASV25 ASV813, and ASV3186 are classified as Pararhizobium
# And ASV2, ASV34, ASV2711, ASV3038,ASV2719  as Agrobacterium

# ASV 4650 can be removed
tax_table(phyl_1)[,6][rownames(tax_table(phyl_1)) %in% c("ASV2","ASV34","ASV2711","ASV3038","ASV2719")]<-"Agrobacterium" 
tax_table(phyl_1)[,6][rownames(tax_table(phyl_1)) %in% c("ASV25","ASV813","ASV3186")]    <-"Pararhizobium"     
tax_table(phyl_1)[,6][rownames(tax_table(phyl_1)) %in% c("ASV13","ASV629")]<-"Peribacillus"          

saveRDS(phyl_1,"Phylo_upd_tax.rds")
