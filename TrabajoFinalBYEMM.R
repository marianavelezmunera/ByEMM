#Cargar paquetes
library(phyloseq)
library(tidyverse)
#Establecer directorio de trabajo
setwd("C:/Users/Maria/OneDrive/Documents/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL/ByEMM")

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
muestras_ASV_francia<-c(row.names(ASV_francia)) #nombres muestras 
ASV_francia<-cbind(ASV_francia,muestras_ASV_francia) #Añadir una columna con el nombre de cada muestra a la tabla de ASV
muestras_selec_francia<-c(metadatos_francia$muestra) #muestras que se van a usar
ASV_francia<-subset(ASV_francia,muestras_ASV_francia%in%muestras_selec_francia)
ASV_francia<-ASV_francia[,-43717] #quitamos la columna de muestras, tabla ASV final
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

ASV_florida_completo<-read.csv("~/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL/ByEMM/florida/ASV_florida.csv") #todos los datos, 387 muestras
ASV_florida_completo<-subset(ASV_florida_completo,Sample.Type=="Healthy"& Site.Type=="Healthy") #dejar solo los sitios e individuos saludables, 53 muestras
ASV_florida_completo<-subset(ASV_florida_completo,Species!="Water") #quitar las muestras de agua, 44 muestras
row.names(ASV_florida_completo)<-c(ASV_florida_completo$Sample.ID)

ASV_florida<-ASV_florida_completo[,-(c(1:5))] #dejamos solo las ASV

####### Remover ASV que no están presentes

sumatoria_florida<-apply(ASV_florida,2,FUN = sum) #sumatoria abundancia ASV en todas las muestras
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

taxa_florida <- read.csv("~/Cosas de la maestría/PRIMER SEMESTRE/BIOGEOQUÍMICA Y ECOLOGÍA MICROBIANA MARINA/TRABAJO FINAL/ByEMM/florida/taxa_florida.csv")

##### SUBMUESTREO ESPECIES DE FRANCIA

#SP1: Herpolitha_limax, SP2: Fungia_fungites

acroporaa<-subset(metadatos_francia,Host_genus=="Acropora")
sample(unique(acroporaa$Species_ID),3)

#Especies Herpolitha_limax, Fungia_fungites, Acropora_intermedia, Acropora_cytherea y Acropora_elseyi

### SUBMUESTREO DE LAS MUESTRAS DE FRANCIA 

Acropora_intermedia<-subset(metadatos_francia,Species_ID=="Acropora_intermedia")
sample(unique(Acropora_intermedia$muestra),2) #Muestras de A. intermedia: "F2-I-05-M" "F2-I-04-M"

Acropora_cytherea<-subset(metadatos_francia,Species_ID=="Acropora_cytherea")
sample(unique(Acropora_cytherea$muestra),2)#Muestras de A. cytherea: "F2-I-16-M" "F2-I-09-M"

Acropora_elseyi<-subset(metadatos_francia,Species_ID=="Acropora_elseyi")
sample(unique(Acropora_elseyi$muestra),2)#Muestras de A.elseyi: B2-I-07" "B2-I-09"

muestras_francia<-c("F2-I-05-M","F2-I-04-M","F2-I-16-M","F2-I-09-M","B2-I-07","B2-I-09","	
F1-I-13","F1-I-14","F1-I-09","F1-I-10")

metadatos_francia_finales<-subset(metadatos_francia,muestra%in%muestras_francia)
metadatos_francia_finales<-rbind(metadatos_francia_finales,subset(metadatos_francia,muestra=="F1-I-13")) #5 especies, 10 muestras, con estos se trabaja

intento<-ASV_francia[c(muestras_francia),]
miremos<-ASV_francia["F1-I-13",]
intento<-rbind(intento,miremos)
intento<-intento[-7,]
ASV_francia_final<-intento
rm(intento)

##### SUBMUESTREO MUESTRAS DE FLORIDA

Colpophyllia_natans<-subset(metadatos_florida,Especies=="Colpophyllia_natans")
sample(unique(Colpophyllia_natans$Sample.ID),2) #Muestras de Colpophyllia_natans: "SWG500.H","SWG502.H"

Orbicella_faveolata<-subset(metadatos_florida,Especies=="Orbicella_faveolata")
sample(unique(Orbicella_faveolata$Sample.ID),2) #Muestras de Orbicella_faveolata: "SWG533.H", "SWG504.H"

Montastraea_cavernosa<-subset(metadatos_florida,Especies=="Montastraea_cavernosa")
sample(unique(Montastraea_cavernosa$Sample.ID),2) #Muestras de Montastraea_cavernosa: SWG508.H" "SWG506.H"

Siderastrea_siderea<-subset(metadatos_florida,Especies=="Siderastrea_siderea")
sample(unique(Siderastrea_siderea$Sample.ID),2) #Muestras de Siderastrea_siderea: "SWG510.H" "SWG509.H"

Pseudodiploria_strigosa<-subset(metadatos_florida,Especies=="Pseudodiploria_strigosa")
sample(unique(Pseudodiploria_strigosa$Sample.ID),2) #Muestras de Pseudodiploria_strigosa: "SWG514.H" "SWG527.H"

muestras_florida<-c("SWG500.H","SWG502.H","SWG533.H","SWG504.H","SWG508.H","SWG506.H", "SWG510.H" ,"SWG509.H","SWG514.H","SWG527.H")

metadatos_florida_finales<-subset(metadatos_florida,Sample.ID%in%muestras_florida)
ASV_florida_final<-ASV_florida[c(muestras_florida),]

#Borrar lo otro y solo dejar en el environment lo que sirve
rm(Acropora_cytherea,Acropora_elseyi,Acropora_intermedia,acroporaa)
rm(ASV_florida,ASV_francia)
rm(Colpophyllia_natans,metadatos_florida)
rm(metadatos_francia)
rm(Montastraea_cavernosa,Orbicella_faveolata,Pseudodiploria_strigosa)
rm(muestras_florida,muestras_francia)
rm(ps_Chiarello_et_al2020)
rm(Siderastrea_siderea)

### ELIMINAR LAS ASV que no están en ninguna muestra (AGAIN)

##Florida

sumatoria_florida_2<-apply(ASV_florida_final,2,FUN = sum) #sumatoria abundancia ASV en todas las muestras
ASV_presentes_florida_2<-rbind(ASV_florida_final,sumatoria_florida_2) #juntar la suma como una última fila (4643 ASV)
sumatoria_florida_2<-as.data.frame(sumatoria_florida_2) #cambio de formato
sumatoria_florida_2<-subset(sumatoria_florida_2,sumatoria_florida_2==0) #dejar solo los que no están en ninguna muestra
asv_eliminar_florida_2<-row.names(sumatoria_florida_2) #cuales ASV se van a eliminar
ASV_presentes_florida_2 <- ASV_presentes_florida_2[,!names(ASV_presentes_florida_2) %in% c(asv_eliminar_florida_2)] #Solo quedan las ASV que si están (2264)
ASV_florida_final<-ASV_presentes_florida_2


## Francia

sumatoria_francia_2<-apply(ASV_francia_final,2,FUN = sum) #sumatoria abundancia ASV en todas las muestras
ASV_presentes_francia_2<-rbind(ASV_francia_final,sumatoria_francia_2) #juntar la suma como una última fila (4643 ASV)
sumatoria_francia_2<-as.data.frame(sumatoria_francia_2) #cambio de formato
sumatoria_francia_2<-subset(sumatoria_francia_2,sumatoria_francia_2==0) #dejar solo los que no están en ninguna muestra
asv_eliminar_francia_2<-row.names(sumatoria_francia_2) #cuales ASV se van a eliminar
ASV_presentes_francia_2 <- ASV_presentes_francia_2[,!names(ASV_presentes_francia_2) %in% c(asv_eliminar_francia_2)] #Solo quedan las ASV que si están (2898)
ASV_francia_final<-ASV_presentes_francia_2

#Quitar del ENV lo que no voy a usar (AGAIN)

rm(ASV_presentes_florida_2,ASV_presentes_francia_2)
rm(sumatoria_florida_2,sumatoria_francia_2)
rm(asv_eliminar_florida_2)
rm(asv_eliminar_francia_2)
rm(miremos)

#Quitar la última fila que son las sumas

ASV_florida_final<-ASV_florida_final[-11,]
ASV_francia_final<-ASV_francia_final[-11,]

#### ANALISIS ALPHA DIVERSITY
#florida

florida_alpha <-microbiome::alpha(t(ASV_florida_final), index = "diversity_shannon")
florida_rich<-richness(t(ASV_florida_final))
metadatos_florida_finales$shannon<-florida_alpha
colnames(metadatos_florida_finales[,6])<-"shannon"

#francia

francia_alpha <-microbiome::alpha(t(ASV_francia_final), index = "diversity_shannon")
francia_rich<-richness(t(ASV_francia_final))
metadatos_francia_finales$shannon<-francia_alpha
colnames(metadatos_francia_finales[,9])<-"shannon"

#Data frame junto
francia_meta<-metadatos_francia_finales[,7:9]
francia_meta$lugar<-c(rep("Francia",nrow(francia_meta)))
colnames(francia_meta)<-c("Especies","Muestra","Shannon","Lugar")

florida_meta<-metadatos_florida_finales[,c(1,5,6)]
florida_meta<-florida_meta %>%
  relocate(Especies,.before = Sample.ID)
florida_meta$lugar<-c(rep("Florida",nrow(florida_meta)))
colnames(florida_meta)<-c("Especies","Muestra","Shannon","Lugar")


## Diferencia diversidad alfa
metadatos_ambos<-rbind(francia_meta,florida_meta)
test_sitio<-metadatos_ambos[,3:4]
spl <- split(metadatos_ambos$Shannon, metadatos_ambos$Lugar)
pv <- ks.test(spl$Francia, spl$Florida) #p value 


anova_especie<-aov(Shannon~Especies,data=metadatos_ambos)
summary(anova_especie)


graph_alpha_lugar<-ggplot(data = metadatos_ambos,aes(x=Lugar,y=Shannon))+
  geom_boxplot()
graph_alpha #Organizarla luego

graph_alpha_espe<-ggplot(data = metadatos_ambos,aes(x=Especies,y=Shannon))+
  geom_boxplot()
graph_alpha_espe #Organizar luego
#NO HAY DIFERENCIAS NI POR SITIO NI POR ESPECIE

## BETA DIVERSIDAD 
#Pasar todo a phyloseq
OTU_francia<-otu_table(ASV_francia_final,taxa_are_rows = FALSE)
OTU_florida<-otu_table(ASV_florida_final,taxa_are_rows = FALSE)

taxa_florida2<-taxa_florida
rownames(taxa_florida2)<-taxa_florida2$taxcnat
taxa_florida2$taxcnat<-NULL
taxa_florida2$seqs.ps<-NULL
colnames(taxa_florida2)[1]<-"Domain"
taxa_florida2<-taxa_florida2[,-7]
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
#Objetos phyloseq listos
francia_ps<-phyloseq(OTU_francia, taxa_francia,meta_francia)
florida_ps<-phyloseq(OTU_florida,taxa_florida_ps,meta_florida)
francia_ps
florida_ps
florida_ps@tax_table$Species<-NULL

rm(taxa_florida_ps,taxa_francia)
# Hagamos un phyloseq con las dos cosas
intento_merge<-merge_phyloseq(francia_ps,florida_ps)
intento_merge
View(intento_merge@sam_data)
View(intento_merge@tax_table)
#POR FIN DIO ESTA MIERDA
View(taxa_francia)

otus_juntas<-as.data.frame(intento_merge@otu_table)
library(ade4)
matrixdist<-vegdist(otus_juntas,method="bray",na.rm = TRUE)
matrixdist<-as.matrix(matrixdist)

PCoA<-pcoa(matrixdist, correction="none", rn=NULL)
biplot(PCoA) #BetaDiversidad


### CORE MICROBIOME 

core <- microbiome::transform(intento_merge, "compositional")
taxa_names(core)[1:2]

core.taxa.standard <- core_members(core, detection = 0.0001, prevalence = 50/100)
core.taxa.standard #Este es el core (57 taxa)

tax_table(core)[tax_table(core) == "k__"] <- NA
tax_table(core)[tax_table(core) == "p__"] <- NA
tax_table(core)[tax_table(core) == "c__"] <- NA
tax_table(core)[tax_table(core) == "o__"] <- NA
tax_table(core)[tax_table(core) == "f__"] <- NA
tax_table(core)[tax_table(core) == "g__"] <- NA
tax_table(core)[tax_table(core) == "s__"] <- NA

tax_table(core)[, colnames(tax_table(core))] <- gsub(tax_table(core)[, colnames(tax_table(core))],  pattern = "[a-z]__", replacement = "")

core.f <- microbiome::add_besthit(core)

taxa_names(core.f)[1:10] 

miremos<-as.data.frame(core@tax_table)
colnames(miremos)[1]<-"Domain"

miremos<-as.matrix(miremos)
miremos<-tax_table(miremos)

core<-phyloseq(core@otu_table,miremos,core@sam_data)


core.taxa.standard <- core_members(core.f, detection = 0.0001, prevalence = 50/100)
core.taxa.standard

core.core <- core(core.f, detection = 0.0001, prevalence = .5)
core.taxa <- taxa(core.core)

tax.mat <- tax_table(core.core)
tax.df <- as.data.frame(tax.mat)

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
View((head(core.taxa.class))) #TABLA CON EL CORE MICROBIOME
