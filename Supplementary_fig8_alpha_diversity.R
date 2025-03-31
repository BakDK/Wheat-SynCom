# Supp fig S8.
# Determine the alpha diversity after rarefaction 100 times, as proposed by Schloss 2024 msphere.
library(phyloseq)
library(tidyverse)

phy_ob<-readRDS("phyloseq_final.rds")

sample_data(phy_ob)
min(sample_sums(phy_ob))
sample_sums(phy_ob) # 20,040 in SC29. 
#Calculating Shannon index 

# Create a matrix with all the samples as rows and 100 empty columns
shan_ra<-matrix(NA,nsamples(phy_ob),100)

# Then rarefy once and then calculate  Shannon diversity, repeat 100 times. Do not set rngseed as the values will the be the same
for(i in 1:100){
  {Full16_rare<-rarefy_even_depth(phy_ob, sample.size = min(sample_sums(phy_ob)))
  resultat<-plot_richness(Full16_rare ,x = "Sample_ID", measures = "Shannon")
  shan_ra[,i]<-resultat$data$value}
}

print(shan_ra)

#Calculate mean of the 100 rarefied Shannon indices
shan_ra_me<-apply(shan_ra,MARGIN =  1, mean)

shan_mean_met<-cbind(sample_data(phy_ob),Shannon =shan_ra_me)

#  Plot the data
shan_mean_met$Plot_var<-paste(shan_mean_met$Variety,shan_mean_met$Treatment, sep = "_")

#Modify the Plot_var variable
shan_mean_met$Plot_var[shan_mean_met$Plot_var %in% "Sheriff_Control"]<-"NatCom"
shan_mean_met$Plot_var[shan_mean_met$Plot_var %in% "Sheriff_SynCom"]<-"SynCom"
shan_mean_met$Plot_var[shan_mean_met$Plot_var %in% "Sheriff_Field"]<-"Field Sheriff"
shan_mean_met$Plot_var[shan_mean_met$Plot_var %in% "Heerup_Field"]<-"Field Heerup"
shan_mean_met$Plot_var[shan_mean_met$Plot_var %in% "Sheriff_Soil_control"]<-"Gamma-irradiated soil"


Shann_plot_all_points<-shan_mean_met %>% 
  filter(!Treatment %in% "Mock" & !Compartment %in% "Mock") %>% 
  ggplot(aes(x = Days, y = Shannon, color = Plot_var))+
  geom_point(position =position_dodge(width = 2.5))+theme_bw()+
  labs(color = "Sample")+theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 2))
#Please notice, the Sheriff_Field and Heerup_Field at days 0 are soil samples from the field  and are not specific for the cultivar.

x11(width = 6, height = 6)
Shann_plot_all_points
ggsave("Plots/Suppl_fig_Shannon_div.png")
dev.off()
