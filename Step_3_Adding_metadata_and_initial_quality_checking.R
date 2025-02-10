# This script describes how metadata was added and quality filtered of samples of the final, combined data.

phyloseq_obj<-readRDS("Syncom_all/phyloseq_all_SC.rds")

# Importing metadata
library(xlsx)
setwd(".."); 
metadata<-read.xlsx("Field_amplicons/Field_metadata.xlsx", sheetIndex = 1)
rownames(metadata)<-metadata$Sample.name
colnames(metadata)[2]<-"Sample_ID"
metadata<-metadata[,-1]
metadata
library(tidyverse)
metadatab1_2<-read.xlsx("Merging_batch1_and2/Metadata.xlsx", sheetIndex = 1)
rownames(metadatab1_2)<-metadatab1_2$Sample_ID
metadatab1_2<-metadatab1_2 %>% select(c(Sample_ID,Pot_number,Treatment,Medium, Compartment, Timepoint, Days))
# Add batch to metadata for later investigation: 

metadatab1_2$Batch[metadatab1_2$Sample_ID %in% c("IC1","IC2","IC3","J3","K10","K16","K38","MOCK","SC11","SC12","SC13","SC20","SC29","SC31","SC9" )]<-"Batch2"

metadatab1_2$Batch[!metadatab1_2$Sample_ID %in% c("IC1","IC2","IC3","J3","K10","K16","K38","MOCK","SC11","SC12","SC13","SC20","SC29","SC31","SC9" )]<-"Batch1"

# Modify tables so they are better aligned:
colnames(metadata)[2]<-"Timepoint"
metadata$Batch<-"Field"
metadata$Treatment<-"Field"
metadata$Compartment<-"Rhizoplane"
metadata$Compartment[metadata$Timepoint %in% "T1"]<-"Soil"
metadata$Compartment[metadata$Timepoint %in% "Mock"]<-"Mock"

#Insert days 
metadata$Days[metadata$Timepoint %in% "T1"]<-0
metadata$Days[metadata$Timepoint %in% "T2"]<-9
metadata$Days[metadata$Timepoint %in% "T3"]<-16
metadata$Days[metadata$Timepoint %in% "T5"]<-23
metadata$Days[metadata$Timepoint %in% "T6"]<-30
metadata$Pot_number<-NA
metadata$Medium<-NA
metadata$Variety[metadata$Variety %in% "V2"]<-"Heerup"
metadata$Variety[metadata$Variety %in% "V1"]<-"Sheriff"

metadatab1_2$Plot<-NA
metadatab1_2$Variety<-"Sheriff"
#Combine the two sheets
metadata_all<-bind_rows(metadata,metadatab1_2)
library(phyloseq)
#The sample names need to be edited to make them match with the metadata
sample_names(phyloseq_obj)[95]<-"Field_Mock"
gsub("_lib*",sample_names(phyloseq_obj)[1:59]," ")
gsub("\._lib*","A_lib2","'_")

dummy_meta<-as.matrix(sample_names(phyloseq_obj))
colnames(dummy_meta)<-"Sample_ID"

dummy_sampl<-sample_data(data.frame(dummy_meta))
sample_names(dummy_sampl)<-sample_names(phyloseq_obj)
phylo_incl_dummy<-merge_phyloseq(phyloseq_obj,dummy_sampl)

# #remove the K25_lib663783_8373 due to very low read number 
phylo_red<-subset_samples(phylo_incl_dummy, !Sample_ID %in% c("K25_lib663793_8373","IC1_lib741308_10436","IC3_lib741310_10441",
                                                   "SC29_lib741305_10441", "SC20_lib741307_10436","K38_lib741303_10436"))

sample_names(phylo_red)[sample_names(phylo_red) %in% c("SC35_1_lib666015","SC35_2_lib663796")]<-c("SC35.1_lib666015_7","SC35.2_lib663796_8")
sample_names(phylo_red)

sample_names(phylo_red)[1:53]<-sapply(strsplit(basename(sample_names(phylo_red)[1:53]),"_"), function(x) paste(x[1],collapse = "_"))

# At this point there are 89 samples in the phyloseq object and 95 in the metadata. 
setdiff(sample_names(phylo_red),metadata_all$Sample_ID) # k19 is different
#Change to K19
sample_names(phylo_red)[sample_names(phylo_red) %in% c("k19")]<-"K19"
#The other way
setdiff(metadata_all$Sample_ID,sample_names(phylo_red)) #Field_NTC, G1,G2,G3,NTC1,NTC2 
# These have failed sequencing
metadata_red<-metadata_all[!metadata_all$Sample_ID %in% c("Field_NTC","G1","G2","G3","NTC1","NTC2"),]
# Add this metadata

sample_data(phylo_red)<-sample_data(metadata_red)

# In a previous r script called "Data_screening_ampvis_object_etc.R" No difference between sc35.1 and sc35.2, so the sc35.2 was deleted
#  SC34 removed as well
phylo_red2<-subset_samples(phylo_red, !Sample_ID %in% c("SC35.2","SC34"))

sample_names(phylo_red2)[51]<-"SC35"
sample_data(phylo_red2)$Sample_ID[sample_data(phylo_red2)$Sample_ID %in% "SC35.1"]<-"SC35"
saveRDS(phylo_red2,"Syncom_all/phyloseq_incl_meta.rds")
