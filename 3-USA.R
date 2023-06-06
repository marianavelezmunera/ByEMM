#Florida data clean-up

#ASV table

ASV_florida_completo<-subset(ASV_florida_completo,Sample.Type=="Healthy"& Site.Type=="Healthy") # Only healthy individuals and sites
ASV_florida_completo<-subset(ASV_florida_completo,Species!="Water") # removing water samples (44 samples)
row.names(ASV_florida_completo)<-c(ASV_florida_completo$Sample.ID)

ASV_florida<-ASV_florida_completo[,-(c(1:5))] #Only ASV

####### Remove non-present ASV

sumatoria_florida<-apply(ASV_florida,2,FUN = sum) #Sum ASV abundance in every sample
ASV_presentes_florida<-rbind(ASV_florida,sumatoria_florida) #Joint sum as the last row
sumatoria_florida<-as.data.frame(sumatoria_florida) #Formating
sumatoria_florida<-subset(sumatoria_florida,sumatoria_florida==0) #Non-present ASV
asv_eliminar_florida<-row.names(sumatoria_florida) #ASV to remove
ASV_presentes_florida <- ASV_presentes_florida[,!names(ASV_presentes_florida) %in% c(asv_eliminar_florida)] #Only present ASV (4643)
ASV_florida<-ASV_presentes_florida #Final name
ASV_florida<-ASV_florida[-45,]

#remove unuseful objects
rm(ASV_presentes_florida)
rm(sumatoria_florida)
rm(asv_eliminar_florida)

#Metadata
#Changing spp names
cambios<-str_replace(ASV_florida_completo$Species,"CNAT","Colpophyllia natans")
cambios<-str_replace(cambios,"OFAV","Orbicella faveolata")
cambios<-str_replace(cambios,"PSTR","Pseudodiploria strigosa")
cambios<-str_replace(cambios,"MCAV","Montastraea cavernosa")
cambios<-str_replace(cambios,"SSID","Siderastrea siderea")
cambios<-str_replace(cambios," ","_") #vector with spp names
nombres_fil<-ASV_florida_completo$Sample.ID #Column with sample names
row.names(ASV_florida_completo)<-nombres_fil #row names
metadatos_florida<-ASV_florida_completo[,1:5] #Only 5 first columns are metadata
metadatos_florida$Especies<-cambios #new col with spp names
metadatos_florida<-metadatos_florida[,-3] #removing column with abreviated names

#remove unuseful objects
rm(ASV_florida_completo)
rm(cambios)
rm(nombres_fil)

#Tax table

taxa_florida <- read.csv("~/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL/ByEMM/florida/taxa_florida.csv")


