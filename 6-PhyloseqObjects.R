# Phyloseq objects

OTU_francia<-otu_table(ASV_francia_final,taxa_are_rows = FALSE)
OTU_florida<-otu_table(ASV_florida_final,taxa_are_rows = FALSE)

taxa_florida2<-taxa_florida
rownames(taxa_florida2)<-taxa_florida2$taxcnat
taxa_florida2$taxcnat<-NULL
taxa_florida2$seqs.ps<-NULL
colnames(taxa_florida2)[1]
colnames(taxa_florida2)[1]<-"Domain"
taxa_florida2<-taxa_florida2[,-7]
taxa_florida2<-as.matrix(taxa_florida2)
taxa_florida_ps<-tax_table(as.matrix(taxa_florida2))

meta_francia<-sample_data(metadatos_francia_finales)
meta_florida<-sample_data(metadatos_florida_finales)

meta_florida<-meta_florida[,c(1,5)]
meta_florida$Lugar<-c(rep("Florida",nrow(meta_florida)))
colnames(meta_florida)<-c("Muestra","Especies","Lugar")
meta_florida<-sample_data(meta_florida)

meta_francia<-meta_francia[,7:8]
meta_francia$Lugar<-c(rep("Francia",nrow(meta_francia)))
colnames(meta_francia)<-c("Especies","Muestra","Lugar")
meta_francia_1<-as.matrix(meta_francia)
meta_francia_1<-as.data.frame(meta_francia_1)
meta_francia_1<-meta_francia_1 %>%relocate(Especies,.after=Muestra)

meta_francia<-sample_data(meta_francia_1)

#Phyloseq objects 
francia_ps<-phyloseq(OTU_francia, taxa_francia,meta_francia)
florida_ps<-phyloseq(OTU_florida,taxa_florida_ps,meta_florida)
View(francia_ps@tax_table)
View(florida_ps@tax_table)
florida_ps

taxa_francia<-as.data.frame(taxa_francia)
colnames(taxa_francia[1])
colnames(taxa_francia)<-c("Domain","Phylum","Class","Order","Family","Genus")
taxa_francia<-as.matrix(taxa_francia)
taxa_francia<-tax_table(taxa_francia)
rm(taxa_florida_ps,taxa_francia)

# Merged phyloseq
intento_merge<-merge_phyloseq(francia_ps,florida_ps)


