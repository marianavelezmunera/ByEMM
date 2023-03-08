#Cargar paquetes
library(phyloseq)
library(tidyverse)
#Establecer directorio de trabajo
setwd("C:/Users/Maria/OneDrive/Documents/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL")

######################### DATOS CORALES FRANCIA ###########################

load("SM3-ps_Chiarello_et_al2020.rData") #datos de los autores sin modificar

##################################### Metadatos ##################################

metadatos<-(ps_Chiarello_et_al2020@sam_data)
coral_metadata<-subset(metadatos,Host_subtype=="coral") #metadatos para samples de coral (82 muestras)
coral_metadata<-subset(coral_metadata,Host_species!="sp") #quitar las que no se sabe qué especie son
coral_metadata<-subset(coral_metadata,Host_species!="sp2")
coral_metadata<-subset(coral_metadata,Host_species!="sp1") #38 muestras
coral_metadata<-subset(coral_metadata,Host_species!="branching") #35 muestras
coral_metadata<-subset(coral_metadata,Host_species!="massive") #30 muestras
coral_metadata$muestra<-row.names(coral_metadata) #Nueva columna con el nombre de la muestra
metadatos_francia<-coral_metadata #Nombre final 
metadatos_francia<-subset(metadatos_francia,Host_genus!="Heteractis") #28 muestras
#quitar los objetos que ya no van a ser útiles
rm(metadatos)
rm(coral_metadata)
unique(metadatos_francia$Species_ID)

########################## Tabla ASV #####################################

ASV_francia<-(ps_Chiarello_et_al2020@otu_table) #OTUs completos
ASV_francia
muestras_ASV_francia<-c(row.names(ASV_francia)) #nombres muestras 
ASV_francia<-cbind(ASV_francia,muestras_ASV_francia) #Añadir una columna con el nombre de cada muestra a la tabla de ASV
muestras_selec_francia<-c(metadatos_francia$muestra) #muestras que se van a usar
ASV_francia<-subset(ASV_francia,muestras_ASV_francia%in%muestras_selec_francia)
ASV_francia<-ASV_francia[,-43717] #quitamos la columna de muestras, tabla ASV final
View(ASV_francia)
####### Remover ASV que no están presentes

ASV_francia_1<-as.data.frame(ASV_francia) #Cambio formato
ASV_francia_1<-apply(ASV_francia_1, 2, as.numeric) #cambio a numeros 
row.names(ASV_francia_1)<-row.names(ASV_francia) #nombre filas
sumatoria_francia<-apply(ASV_francia_1,2,FUN = sum) #sumatoria abundancia ASV en todas las muestras
ASV_presentes_francia<-rbind(ASV_francia_1,sumatoria_francia)
sumatoria_francia<-as.data.frame(sumatoria_francia) #(43K ASV)
sumatoria_francia<-subset(sumatoria_francia,sumatoria_francia==0) #las que no están presentes
asv_eliminar_francia<-row.names(sumatoria_francia) #vector con las asv que vamos a eliminar
ASV_presentes_francia<-as.data.frame(ASV_presentes_francia)
ASV_presentes_francia<-ASV_presentes_francia[,!names(ASV_presentes_francia) %in% c(asv_eliminar_francia)] #4867 ASV 
ASV_francia<-ASV_presentes_francia
View(sumatoria_francia)
ASV_francia<-ASV_francia[-29,]
rm(muestras_ASV_francia)
rm(muestras_selec_francia)
rm(ASV_francia_1)
rm(ASV_presentes_francia)
rm(asv_eliminar_francia)
rm(sumatoria_francia)
############################### Tabla de taxonomía ###################

taxa_francia<-(ps_Chiarello_et_al2020@tax_table)

##############################  CORALES FLORIDA ################################

############################### ASV FLORIDA ##############################

ASV_florida_completo<-read.csv("~/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL/florida/ASV_florida.csv") #todos los datos, 387 muestras
ASV_florida_completo<-subset(ASV_florida_completo,Sample.Type=="Healthy"& Site.Type=="Healthy") #dejar solo los sitios e individuos saludables, 53 muestras
ASV_florida_completo<-subset(ASV_florida_completo,Species!="Water") #quitar las muestras de agua, 44 muestras
row.names(ASV_florida_completo)<-c(ASV_florida_completo$Sample.ID)

ASV_florida<-ASV_florida_completo[,-(c(1:5))] #dejamos solo las ASV

####### Remover ASV que no están presentes

sumatoria_florida<-apply(ASV_florida,2,FUN = sum) #sumatoria abundancia ASV en todas las muestras
View(sumatoria_florida)
ASV_presentes_florida<-rbind(ASV_florida,sumatoria_florida) #juntar la suma como una última fila (16369 ASV)
sumatoria_florida<-as.data.frame(sumatoria_florida) #cambio de formato
sumatoria_florida<-subset(sumatoria_florida,sumatoria_florida==0) #dejar solo los que no están en ninguna muestra
asv_eliminar_florida<-row.names(sumatoria_florida) #cuales ASV se van a eliminar
ASV_presentes_florida <- ASV_presentes_florida[,!names(ASV_presentes_florida) %in% c(asv_eliminar_florida)] #Solo quedan las ASV que si están (4643)
ASV_florida<-ASV_presentes_florida
ASV_florida<-ASV_florida[-45,]
rm(ASV_presentes_florida)
rm(sumatoria_florida)
rm(asv_eliminar_florida)

############################### METADATOS FLORIDA #############################
#limpiar los nombres de las especies
cambios<-str_replace(ASV_florida_completo$Species,"CNAT","Colpophyllia natans")
cambios
cambios<-str_replace(cambios,"OFAV","Orbicella faveolata")
cambios<-str_replace(cambios,"PSTR","Pseudodiploria strigosa")
cambios<-str_replace(cambios,"MCAV","Montastraea cavernosa")
cambios<-str_replace(cambios,"SSID","Siderastrea siderea")
cambios<-str_replace(cambios," ","_") #vector listo con los nombres de las especies
nombres_fil<-ASV_florida_completo$Sample.ID #para ponerle nombre a las filas
row.names(ASV_florida_completo)<-nombres_fil #nombres de las filas
metadatos_florida<-ASV_florida_completo[,1:5] #Solo las 5 primeras columnas son metadatos
metadatos_florida$Especies<-cambios #nueva columna con el nombre de las especies bien
metadatos_florida<-metadatos_florida[,-3] #quitamos la col de los nombres viejos
rm(ASV_florida_completo)
rm(cambios)
rm(nombres_fil)

############################### TAXA FLORIDA ############################

taxa_florida<-read.csv("~/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL/florida/taxa_florida.csv") 
