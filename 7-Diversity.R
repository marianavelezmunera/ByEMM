# Diversity
# Alpha diversity
#Florida

florida_alpha <-microbiome::alpha(t(ASV_florida_final), index = "diversity_shannon") #Alpha diversity
florida_rich<-richness(t(ASV_florida_final)) #Richness
metadatos_florida_finales$shannon<-florida_alpha
colnames(metadatos_florida_finales[,6])<-"shannon"

#Francia

francia_alpha <-microbiome::alpha(t(ASV_francia_final), index = "diversity_shannon") #Alpha diversity
francia_rich<-richness(t(ASV_francia_final)) #Richness
metadatos_francia_finales$shannon<-francia_alpha
colnames(metadatos_francia_finales[,9])<-"shannon"


# Beta diversity

matrixdist<-vegdist(otus_juntas,method="bray",na.rm = TRUE) #Distance matrix
matrixdist<-as.matrix(matrixdist)
PCoA2<-dudi.pco(d = matrixdist, scannf = FALSE, nf = 2)
