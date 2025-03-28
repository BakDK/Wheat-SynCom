# Plotting figure 2 and 
# Input files from "step_5_QC_and_filtering.R"
library(ampvis2)
library(tidyverse)
comb_amp<-readRDS("total_ampvis.rds")
colSums(comb_amp$abund)

# what is the range of observations in each sample?
colSums(comb_amp$abund) %>% range() # from 20k to 166k reads

# Start by extracting SynCom ASVs from the Inoculation samples and then plot only these
SC_start<-comb_amp %>% amp_subset_samples(Treatment %in% "SynCom" & Timepoint %in% "T0")
SC_start_heat<-SC_start %>% amp_heatmap(tax_show = 30,
                         tax_aggregate = "Genus",
                         round = 2)
ggsave("Plots/Syncom_starting.png",SC_start_heat, width = 7, height = 7) # Supplementary figure 7


SC_ASVs<-SC_start$tax$OTU

#Subset data based only on syncom ASVs
Amp_syncom_ASVs<-comb_amp %>% amp_subset_taxa(tax_vector = SC_ASVs, normalise = TRUE)

# retrieve relative abundances and then plot figure 2
AMP_textfile_ASVs<-Amp_syncom_ASVs %>%
  amp_subset_samples(Treatment %in% c("SynCom", "Soil_control")) %>%
  amp_subset_samples(Compartment %in% c("Rhizoplane","0")) %>%
  amp_heatmap(group_by ="Sample_ID",facet_by = "Days",
              tax_aggregate = "Genus", tax_show = 13, tax_add = "Family", round = 2,
              normalise = FALSE, textmap = TRUE)

AMP_df<-cbind(rownames(AMP_textfile_ASVs),AMP_textfile_ASVs[,-1:-3]) 
colnames(AMP_df)[1]<-"Taxa"
AMP_DF_Long<-AMP_df%>% pivot_longer(cols=!Taxa,names_to = "Sample", values_to = "Rel_abu")
AMP_DF_Long_com<-cbind(AMP_DF_Long,str_split_fixed(AMP_DF_Long$Sample, " ", 2))
colnames(AMP_DF_Long_com)[4:5]<-c("Sample_ID","Day")

AMP_DF_Long_com2<-cbind(AMP_DF_Long_com,str_split_fixed(AMP_DF_Long$Taxa, ";", 2))
colnames(AMP_DF_Long_com2)[6:7]<-c("Famiy","Genus")

AMP_DF_Long_com2_r<-AMP_DF_Long_com2 %>% filter(!Day %in% "0")
#Change order of days
levels(AMP_DF_Long_com2$Day)
AMP_DF_Long_com2_r$Day<-factor(AMP_DF_Long_com2_r$Day, levels =c("7","14","21","28"))

levels(AMP_DF_Long_com2_r$Genus)
AMP_DF_Long_com2_r$Genus<-gsub(" ","",AMP_DF_Long_com2_r$Genus)
AMP_DF_Long_com2_r$Genus<-factor(AMP_DF_Long_com2_r$Genus, levels = c("Agrobacterium", "Ensifer","Flavobacterium",
                                                                      "Microbacterium","Pararhizobium","Pedobacter",
                                                                      "Stenotrophomonas","Variovorax",
                                                                      "Arthrobacter","Bacillus","Pseudomonas","Paenarthrobacter",
                                                                      "Peribacillus"))

# Make a plot facetting based on taxa
first_ed_plot<-AMP_DF_Long_com2_r %>% 
  ggplot(aes(x = Day, y = Rel_abu))+geom_point(alpha = 0.2)+facet_wrap(~Genus, scales ="free_y", ncol = 4)+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        panel.background = element_rect(fill = "#fdf8ef"))+
  labs(y="Relative Abundance (%)")

AMP_DF_mean<-AMP_DF_Long_com2_r %>% 
  group_by(Genus,Day) %>%
  summarise( Mean = mean(Rel_abu))

sec_ed_pl<-first_ed_plot + geom_point(AMP_DF_mean, mapping= aes(x = Day, y = Mean),color = "Blue", shape = 4, size = 3)+
   theme(legend.position = "bottom",axis.text = element_text(size = 10),
         strip.text = element_text(size = 12,face = "italic"))

x11(height = 7, width = 8.2)
sec_ed_pl
ggsave("Plots/figure_2.png", dpi = 300) 
dev.off()

#Define the syncom members

Syncom_genera<-c("Microbacterium","Arthrobacter","Paenarthrobacter","Bacillus", "Peribacillus",
                 "Flavobacterium","Pedobacter","Agrobacterium","Ensifer","Pararhizobium", "Variovorax",
                 "Pseudomonas","Stenotrophomonas")


#Plotting the most abundant genera over time in the SynCom inoculated plants (Figure S9)
top_15_gen_sc_samples_heat<-comb_amp %>%
  amp_subset_samples(Treatment %in% c("SynCom")) %>%
  amp_subset_samples(Compartment %in% c("Rhizoplane")) %>%
  amp_heatmap(group_by ="Days",
              tax_aggregate = "Genus", tax_show = 15,  round = 2)
              
ggsave("Plots/syncom_heatmap_top_all_gen.png",top_15_gen_sc_samples_heat,width = 10, height = 8)

# Plot the Natcom for comparison as a suppl figure S10
top_15_gen_NC_samples_heat<-comb_amp %>%
  amp_subset_samples(Treatment %in% c("Control")) %>%
  amp_subset_samples(Compartment %in% c("Rhizoplane")) %>%
  amp_heatmap(group_by ="Days",
              tax_aggregate = "Genus", tax_show = 15,  round = 2)

ggsave("Plots/natcom_heatmap_top_all_gen.png",top_15_gen_NC_samples_heat,width = 10, height = 8) 

#Plot the starting relative abundance of each of the 13 genera in the starting inoculum of the NatCom, SynCom and soil
Amp_syncom_gene$metadata
Amp_syncom_gene$metadata$Plot<-paste(Amp_syncom_gene$metadata$Variety, Amp_syncom_gene$metadata$Treatment, sep = "_")

#Plot a supporting figure S12: 
AMP_start<-Amp_syncom_gene %>% 
  amp_subset_samples(Days %in% "0") %>%
  amp_subset_samples(!Treatment %in% "Soil_control")
AMP_start$metadata$Treatment[AMP_start$metadata$Treatment %in% "Control"]<-"NatCom"
Starting_SC_Members_all_samples<-AMP_start %>% 
  amp_heatmap(group_by = "Treatment",
              tax_show = 13,
              tax_aggregate = "Genus", 
              round = 2, normalise = FALSE)
ggsave("Plots/SC_gen_start_rel_abundance.png",Starting_SC_Members_all_samples, width = 8, height = 8 )



#Plot only the gamma irradiated soil ( Figure S5)
sterile_soil_plot<-comb_amp %>%
  amp_subset_samples(Treatment %in% "Soil_control") %>%
  amp_heatmap(tax_aggregate =  "Genus")

ggsave("Plots/Sterile_soil_heat.png",sterile_soil_plot,width = 5,height =4)


##### NB - might be deleted


#Use the same to get the total abundance over time
SC_Members_over_time_all_samples<-Amp_syncom_gene %>% 
  amp_subset_samples(!Days %in% "0") %>%
  amp_subset_samples(!Treatment %in% "Mock") %>%
   amp_heatmap(group_by = "Days",
               facet_by = "Plot",
              tax_show = 13,
              tax_aggregate = "Kingdom", 
              round = 1, normalise = FALSE)
ggsave("test_plots/SC_gen_total_over_time.png",SC_Members_over_time_all_samples, width = 7, height = 4) # This might not be included after all.


