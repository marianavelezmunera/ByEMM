# France data clean-up
# Metadata

metadatos<-(ps_Chiarello_et_al2020@sam_data)

coral_metadata<-subset(metadatos,Host_subtype=="coral") #metadata just for coral samples (82 samples)
coral_metadata<-subset(coral_metadata,Host_species!="sp") #Eliminating samples that weren't identified to spp level
coral_metadata<-subset(coral_metadata,Host_species!="sp2")
coral_metadata<-subset(coral_metadata,Host_species!="sp1") #38 samples
coral_metadata<-subset(coral_metadata,Host_species!="branching") #35 samples
coral_metadata<-subset(coral_metadata,Host_species!="massive") #30 samples
coral_metadata$muestra<-row.names(coral_metadata) #sample ID column
metadatos_francia<-coral_metadata # Final object name 
metadatos_francia<-subset(metadatos_francia,Host_genus!="Heteractis") #28 samples

#remove unusefull objects 
rm(metadatos)
rm(coral_metadata)

# ASV table

ASV_francia<-(ps_Chiarello_et_al2020@otu_table) #All ASV
muestras_ASV_francia<-c(row.names(ASV_francia)) #Add sample names
ASV_francia<-cbind(ASV_francia,muestras_ASV_francia) #Add column with sample names
muestras_selec_francia<-c(metadatos_francia$muestra) #Filtering the selected samples
ASV_francia<-subset(ASV_francia,muestras_ASV_francia%in%muestras_selec_francia)
ASV_francia<-ASV_francia[,-43717] #Remove sample names column, final ASV table

# Remove non-present ASV

ASV_francia_1<-as.data.frame(ASV_francia) #Format change
ASV_francia_1<-apply(ASV_francia_1, 2, as.numeric) #Change to numeric 
row.names(ASV_francia_1)<-row.names(ASV_francia) # Row names
sumatoria_francia<-apply(ASV_francia_1,2,FUN = sum) # Sum ASV abundance in every sample
ASV_presentes_francia<-rbind(ASV_francia_1,sumatoria_francia)
sumatoria_francia<-as.data.frame(sumatoria_francia) #43K ASV
sumatoria_francia<-subset(sumatoria_francia,sumatoria_francia==0) #Non-present ASV
asv_eliminar_francia<-row.names(sumatoria_francia) #Vector with ASVs to remove
ASV_presentes_francia<-as.data.frame(ASV_presentes_francia)
ASV_presentes_francia<-ASV_presentes_francia[,!names(ASV_presentes_francia) %in% c(asv_eliminar_francia)] #4867 ASV 
ASV_francia<-ASV_presentes_francia #Final name
ASV_francia<-ASV_francia[-29,]

#remove unusefull objects 
rm(muestras_ASV_francia)
rm(muestras_selec_francia)
rm(ASV_francia_1)
rm(ASV_presentes_francia)
rm(asv_eliminar_francia)
rm(sumatoria_francia)

#Tax table

taxa_francia<-(ps_Chiarello_et_al2020@tax_table)
