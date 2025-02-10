# This script is a follow up on the "Comparing_full_genomes_with_ASVs_from_SC.R"
# Check of controls, starting inoculum and mock samples. 
# removal of mitochondria, chloroplast and non-bacterial ASVs.

library(phyloseq)
library(ampvis2)

phylo_obj<-readRDS("Phylo_upd_tax.rds")

otu_table(phylo_obj)<-t(otu_table(phylo_obj))
# due to merging, ASVs not found in one object, but another phyloseq object has gotten NA
otu_table(phylo_obj)[is.na(otu_table(phylo_obj))]<-0
sample_sums(phylo_obj) 

# Convert to ampvis2 object

otu_table(phylo_obj)

rownames(tax_table(phylo_obj))

df<-amp_load(otutable = otu_table(phylo_obj),
             taxonomy = tax_table(phylo_obj),
             metadata = data.frame(sample_data(phylo_obj))) #115 ASVs with 0 in all samples were removed
#First remove the chloroplast and mitochondria
#need to check why some are lost
any(taxa_sums(phylo_obj)==0) 
sum(taxa_sums(phylo_obj)==0) #115 are 0

total_amp<-amp_subset_taxa(df, tax_vector = c("Chloroplast","Mitochondria"),remove = TRUE) # 19 ASVs removed. All mitochondria
total_amp<-total_amp %>% amp_subset_taxa( tax_vector = c("Bacteria")) # Another 7 ASVs removed, this is now 12,559 ASVs
total_amp %>% amp_subset_taxa( tax_vector = c("Bacteria")) 
phylum_vec<-unique(total_amp$tax$Phylum) #some are unclassified
phylum_vec<-phylum_vec[phylum_vec != ""]
ladmigse<- total_amp %>% amp_subset_taxa(tax_vector = phylum_vec, remove = TRUE) 
ladmigse$abund # Either low abundant or found in negative controls.

total_amp<-total_amp %>% amp_subset_taxa(tax_vector = phylum_vec) #Another 106 ASVs are removed - total 12,453


total_amp$metadata$Treatment[total_amp$metadata$Compartment %in% "Mock"]<-"Mock"

#Have a look at the mock samples
total_amp %>% amp_subset_samples(Treatment %in% "Mock")%>%
  amp_heatmap(group_by = "Sample_ID",
              facet_by = "Batch",
              tax_show = 40,
              tax_aggregate = "Species",
              tax_add = "Genus")
# The Lactobacillus is a bit underrepresented compared to theoretical values (18.4%)
# Salmonella is split into phage and enterica.
# THe kosakonia is closely related to Enteroccous and Eschericiea (it was previously unknown family)
# This is good

# Gamma irradiated soil. Those called J1,J2 and J3
total_amp %>% amp_subset_samples(Treatment %in% "Soil_control")%>%
  amp_heatmap(group_by = "Sample_ID",
              facet_by = "Batch",
              tax_show = 40,
              tax_aggregate = "Species",
              tax_add = "Genus")
# This also looks really good (compared with first_exp_syncom.html)
total_amp %>% amp_subset_samples(Treatment %in% "Soil_control")%>%
  amp_heatmap(group_by = "Sample_ID",
              facet_by = "Batch",
              tax_show = 40,
              tax_aggregate = "Genus",
              tax_add = "Family")
# Still good

## Starting syncom called In1, IN2, IN3
total_amp %>% amp_subset_samples(Treatment %in% "SynCom")%>%
  amp_subset_samples(Compartment %in% "0")%>%
  amp_heatmap(group_by = "Sample_ID",
              facet_by = "Batch",
              tax_show = 40,
              tax_aggregate = "OTU",
              tax_add = "Genus")

# This also looks good. 

# save the ampvis object

saveRDS(total_amp,"total_ampvis.rds")
total_amp$metadata

# Also update the phyloseq object so it does not contain chloroplast and mitochondria
# start with 12,700 ASVs
phylo_minchl<-subset_taxa(phylo_obj,!Family =="Mitochondria"|is.na(Family)) # 20 ASVs removed
phylo_minchl<-subset_taxa(phylo_minchl,!Class =="Chloroplast"|is.na(Class)) # none as seen for the ampvis 2 object
phylo_mito<-subset_taxa(phylo_obj,Family =="Mitochondria") # 20 ASVs


taxphylo_upd<-subset_taxa(phylo_minchl,Kingdom =="Bacteria") # removes 7, which is as seen in ampvis2 object
phylo_fina<-prune_taxa(taxa_sums(taxphylo_upd) > 0, taxphylo_upd)

#Remove ASVs not classified at the phylum level
phyloseq_fin<-subset_taxa(phylo_fina, !Phylum =="")
phyloseq_fin<-prune_taxa(taxa_sums(phyloseq_fin)>0, phyloseq_fin)
#12,453 ASVs as in the ampvis
min(sample_sums(phyloseq_fin)) # 20,040 OK
max(sample_sums(phyloseq_fin)) # 166,116
mean(sample_sums(phyloseq_fin)) #77,063.16 ok

saveRDS(phyloseq_fin,"phyloseq_final.rds")

