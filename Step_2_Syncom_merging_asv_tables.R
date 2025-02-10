# combine data sets, processed in step 1a-c.
# The amplicons were sequenced in three runs. Datasets are merged prior to assigning taxonomy.

#First the batch 1, "batch1_dada2.R"
setwd("..")
batch1_asv_table<-readRDS("batch1_syncom_data/asv_table_batch1.rds") 
head(batch1_asv_table)

# Field data, treated in field_dada2.R
field_asv_table<-readRDS("amplicons_field_SC/asv_table_field.rds")

head(field_asv_table)

# batc2 data, treated in batch2_syncom_16S.R
batch2_asv_table<-readRDS("asv_table.rds")

head(batch2_asv_table)

dim(field_asv_table)
# There are 2697 ASVs in batch2, 3398 ASVs in batch 1 and 9009 in field data set. 
length(setdiff(colnames(batch2_asv_table),colnames(batch1_asv_table))) # 1517 different
length(setdiff(colnames(batch1_asv_table),colnames(batch2_asv_table))) # 2218 different
length(intersect(colnames(batch1_asv_table),colnames(batch2_asv_table))) # 1180 shared
1180+1517
1180+2218

# merging them:
batch1_2_asvs<-dplyr::bind_rows(data.frame(batch2_asv_table),data.frame(batch1_asv_table))
1180+2218+1517 #this fits with the number of columns

#Difference with the field
length(intersect(colnames(batch1_2_asvs),colnames(field_asv_table))) #1221 shared
length(setdiff(colnames(batch1_2_asvs),colnames(field_asv_table))) # 3694 different
length(setdiff(colnames(field_asv_table),colnames(batch1_2_asvs))) # 7785 different

# merge them
total_asvs<-dplyr::bind_rows(batch1_2_asvs,data.frame(field_asv_table))
1221+3694+7785 # In agreement with above

library(dada2)

taxa <- assignTaxonomy(as.matrix(total_asvs), "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print[1:20,]

saveRDS(taxa,"Syncom_all/staxa_field.rds")

#Convert to a phyloseq object (easy to work with)
samples.out <- rownames(total_asvs)

library(phyloseq)
ps <- phyloseq(otu_table(total_asvs, taxa_are_rows=FALSE), 
               tax_table(taxa))

ps
dna <- Biostrings::DNAStringSet(taxa_names(ps)) #Here the DNA sequences are stored in a separate object, which is then put into the phyloseq object. It might be usefull later on in an anlysis.
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) #Change taxa names to ASV1, ASV2, etc.
ps

saveRDS(ps,"Syncom_all/phyloseq_all_SC.rds")
