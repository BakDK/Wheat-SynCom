# some import from step3_initial_plots.R is needed
library(ampvis2)
library(tidyverse)
comb_amp<-readRDS("total_ampvis.rds")

# Start by extracting SynCom ASVs from the Inoculation samples and then plot only these
SC_start<-comb_amp %>% amp_subset_samples(Treatment %in% "SynCom" & Timepoint %in% "T0")
SC_ASVs<-SC_start$tax$OTU

#Subset data based only on syncom ASVs
Amp_syncom_ASVs<-comb_amp %>% amp_subset_taxa(tax_vector = SC_ASVs, normalise = TRUE)

#Define the syncom members

Syncom_genera<-c("Microbacterium","Arthrobacter","Paenarthrobacter","Bacillus", "Peribacillus",
                 "Flavobacterium","Pedobacter","Agrobacterium","Ensifer","Pararhizobium", "Variovorax",
                 "Pseudomonas","Stenotrophomonas")
# Subsetting the full dataset based on the SynCom 
Amp_syncom_gene<-comb_amp %>% amp_subset_taxa(tax_vector = Syncom_genera, normalise = TRUE)


Amp_syncom_gene$metadata$Treatment[Amp_syncom_gene$metadata$Variety %in% "Heerup"]<-"Field_Heerup"
# Plot a time series of each bacteria

samples_genera_abund<-Amp_syncom_gene %>%
  amp_heatmap(group_by = "Sample_ID",
              tax_show = 13,
              tax_aggregate = "Genus",
              normalise = FALSE,
              textmap = TRUE) 

samples_inkl_meta<-samples_genera_abund  %>% t() %>% cbind(.,Sample_ID = rownames(.)) %>%
  merge(.,Amp_syncom_gene$metadata, by = "Sample_ID")

samples_inkl_meta_long<-samples_inkl_meta %>% select(-Batch, -Medium, -Compartment, -Pot_number,-Plot) %>%
  pivot_longer(!c(Treatment,Timepoint, Days,Sample_ID,Variety), names_to ="Genus", values_to = "Rel_abundance")

samples_inkl_meta_long$Rel_abundance<-as.numeric(samples_inkl_meta_long$Rel_abundance)


df3<-samples_inkl_meta_long%>% filter(!Treatment %in% c("Mock","Soil_control"))%>%
  filter(!Days %in% "0") 
  
df3$Treatment[df3$Treatment %in% "Control"]<-"NatCom"
df3$Treatment[df3$Treatment %in% "Field"]<-"Field Sheriff"
df3$Treatment[df3$Treatment %in% "Field_Heerup"]<-"Field Heerup"

df3$Treatment<-factor(df3$Treatment, levels =c("SynCom","NatCom","Field Sheriff","Field Heerup"))
levels(df3$Treatment)
saveRDS(df3,"Rel_abun_long_format_SC_gen.rds")

library(lemon)
plot_p4<-ggplot(df3,aes(x = Days, y = log10(Rel_abundance), color = Genus))+
  geom_point(alpha = 0.5,  size = 1)+
  facet_rep_grid(Genus~Treatment,scales = "free_y", repeat.tick.labels = FALSE,switch = "y",space = "free")+ 
  geom_vline(xintercept =c(10.5,17.5,24.5))+theme_bw()+
  scale_x_continuous(breaks = c(7,14,21,28), expand = c(0,0), limits = c(2,32))+ 
  scale_y_continuous(n.breaks = 4)+
  geom_smooth(data= subset(df3, Days>10),method = "lm", se = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background =  element_blank(),
        strip.background = element_rect(fill="white", color = "white"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "italic"),
        legend.position = "none")+
  xlab("Time (Days)")+
  ylab("Log10 Relative Abundance")



x11(height = 11, width =10)
plot_p4

ggsave("test_plots/fig5.png",plot_p4, height = 11, width = 10,units = "in",dpi = 300,device = "png")

