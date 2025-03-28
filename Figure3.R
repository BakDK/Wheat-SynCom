# exploration of the data using the ampvis object. Following data quality check in step2
library(ampvis2)
library(tidyverse)
comb_amp<-readRDS("total_ampvis.rds")

# Start by extracting SynCom ASVs from the Inoculation samples and then plot only these
SC_start<-comb_amp %>% amp_subset_samples(Treatment %in% "SynCom" & Timepoint %in% "T0")

SC_Gen<-unique(SC_start$tax$Genus)
sync_asvs<-as.vector(SC_start$tax$OTU)
#Subset data based only on syncom ASVs
Amp_syncom_Gen<-comb_amp %>% amp_subset_taxa(tax_vector = SC_Gen, normalise = TRUE)


# Distribution of ASVs within these genera
## Subset the 13 genera, and then see how big all ASVs belonging to the 13 genera comprise
Amp_syncom_Gen$metadata

syncom_ASVs_df<-Amp_syncom_Gen %>%
  amp_subset_samples(Treatment %in% c("SynCom")) %>%
  amp_subset_samples(!Days %in% c("0")) %>%
  amp_heatmap(facet_by ="Timepoint", group_by ="Pot_number",
              tax_aggregate = "OTU", tax_show = 152, tax_add = "Genus",
              normalise = FALSE,
              textmap = TRUE) %>% t() %>%as.data.frame() 


#Get the metadata
syncom_asvs_met <-Amp_syncom_Gen%>%
  amp_subset_samples(Treatment %in% c("SynCom")) %>%
  amp_subset_samples(!Days %in% c("0"))  %>%
  .$metadata  %>% as.data.frame()
#Combine
syncom_ASVs_df2<-cbind(syncom_ASVs_df,Pot_Time=rownames(syncom_ASVs_df))

syncom_asvs_met$Pot_Time<-paste(syncom_asvs_met$Pot_number,syncom_asvs_met$Timepoint,sep = " ")
syncom_asv<-merge(syncom_ASVs_df2,syncom_asvs_met, by = "Pot_Time")

# pivot longer
syncom_asv_long_tibble <- syncom_asv %>% 
  select(-Pot_Time,-Batch, -Medium, -Compartment, -Pot_number, -Treatment, -Plot, -Variety) %>%
  pivot_longer(!c(Timepoint, Days,Sample_ID), names_to ="Tax_ID", values_to = "Relative_Abundance")

#Add catetogry to the different ASVs depending on whether they originate from SynCom or not
sync_ASV_df<-cbind(str_split_fixed(syncom_asv_long_tibble$Tax_ID,"; ",2),syncom_asv_long_tibble)
colnames(sync_ASV_df)[1:2]<-c("Genus","ASV")
(sync_ASV_df)
sync_ASV_df$SynCom[sync_ASV_df$ASV %in% sync_asvs]<-"SynCom"
sync_ASV_df$SynCom[!sync_ASV_df$ASV %in% sync_asvs]<-"Other_Source"

#Determine the mean relative abundance per day per genus for each category


#sum the relative abundance of each genera per sample
fordeling_ASV<-sync_ASV_df %>%
  group_by(Sample_ID, Genus,SynCom, Days) %>%
  mutate(Genus_per_sample = sum(Relative_Abundance))  %>% ungroup() 

distrib_ASVs <-fordeling_ASV %>%
  group_by( Genus,SynCom, Days) %>%
  mutate(Genus_per_day = mean(Genus_per_sample)) %>% ungroup()

fordel_asv<-distrib_ASVs %>% group_by(Genus, Genus_per_day) %>%
  distinct(Genus_per_day, .keep_all = TRUE)

# Change order of factors
fordel_asv$SynCom<-factor(fordel_asv$SynCom, levels = c("SynCom", "Other_Source"))
# PLot

SC_other_source_plot<-fordel_asv %>%  ggplot( aes(x = Days, y = Genus_per_day, fill = SynCom))+
  geom_bar(position = "stack", stat="identity")+
  facet_wrap(~Genus, scales = "free_y", ncol = 4) + 
  theme_bw()+
  theme(strip.text = element_text(size = 12,face = "italic"),
        panel.grid = element_blank(), strip.background = element_blank(),
        panel.background = element_blank(),
        legend.position.inside=  c(0.6,0.08))+
  guides(fill = guide_legend(direction = "horizontal", position = "inside"))+
  ylab("Mean Relative Abundance")+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(breaks = c(7,14,21,28))+
  scale_fill_manual(values = c("#1ab1cc","#7B0323"), name = "ASV Origin", labels = c("SynCom (Teal)", "Other Source (Red)"))
  
#Plot and export figure 3
x11(width = 8.2, height = 7)
SC_other_source_plot
ggsave("Plots/ASV_origin_fig3.png", dpi = 300)
dev.off()

