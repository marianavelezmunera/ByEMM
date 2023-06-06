# Subsampling

#France
#SP1: Herpolitha_limax, SP2: Fungia_fungites

acroporaa<-subset(metadatos_francia,Host_genus=="Acropora")
sample(unique(acroporaa$Species_ID),3) #Selecting 3 random Acropora species

#Sspecies: Herpolitha_limax, Fungia_fungites, Acropora_intermedia, Acropora_cytherea y Acropora_elseyi

#2 samples per species

Acropora_intermedia<-subset(metadatos_francia,Species_ID=="Acropora_intermedia")
sample(unique(Acropora_intermedia$muestra),2) #A. intermedia samples: "F2-I-05-M" "F2-I-04-M"

Acropora_cytherea<-subset(metadatos_francia,Species_ID=="Acropora_cytherea")
sample(unique(Acropora_cytherea$muestra),2)#A. cytherea samples: "F2-I-16-M" "F2-I-09-M"

Acropora_elseyi<-subset(metadatos_francia,Species_ID=="Acropora_elseyi")
sample(unique(Acropora_elseyi$muestra),2)#A.elseyi samples: B2-I-07" "B2-I-09"

muestras_francia<-c("F2-I-05-M","F2-I-04-M","F2-I-16-M","F2-I-09-M","B2-I-07","B2-I-09","	
F1-I-13","F1-I-14","F1-I-09","F1-I-10") #France samples

metadatos_francia_finales<-subset(metadatos_francia,muestra%in%muestras_francia) #Just metadata for the selected samples
metadatos_francia_finales<-rbind(metadatos_francia_finales,subset(metadatos_francia,muestra=="F1-I-13")) #5 spp, 10 samples

#Florida

Colpophyllia_natans<-subset(metadatos_florida,Especies=="Colpophyllia_natans")
sample(unique(Colpophyllia_natans$Sample.ID),2) #Colpophyllia_natans samples: "SWG500.H","SWG502.H"

Orbicella_faveolata<-subset(metadatos_florida,Especies=="Orbicella_faveolata")
sample(unique(Orbicella_faveolata$Sample.ID),2) #Orbicella_faveolata samples: "SWG533.H", "SWG504.H"

Montastraea_cavernosa<-subset(metadatos_florida,Especies=="Montastraea_cavernosa")
sample(unique(Montastraea_cavernosa$Sample.ID),2) #Montastraea_cavernosa samples: SWG508.H" "SWG506.H"

Siderastrea_siderea<-subset(metadatos_florida,Especies=="Siderastrea_siderea")
sample(unique(Siderastrea_siderea$Sample.ID),2) #Siderastrea_siderea samples: "SWG510.H" "SWG509.H"

Pseudodiploria_strigosa<-subset(metadatos_florida,Especies=="Pseudodiploria_strigosa")
sample(unique(Pseudodiploria_strigosa$Sample.ID),2) #Pseudodiploria_strigosa samples: "SWG514.H" "SWG527.H"

muestras_florida<-c("SWG500.H","SWG502.H","SWG533.H","SWG504.H","SWG508.H","SWG506.H", "SWG510.H" ,"SWG509.H","SWG514.H","SWG527.H") #Final samples

metadatos_florida_finales<-subset(metadatos_florida,Sample.ID%in%muestras_florida) #Just metadata for the selected samples
ASV_florida_final<-ASV_florida[c(muestras_florida),]

#remove unuseful 
rm(Acropora_cytherea,Acropora_elseyi,Acropora_intermedia,acroporaa)
rm(ASV_florida,ASV_francia)
rm(Colpophyllia_natans,metadatos_florida)
rm(metadatos_francia)
rm(Montastraea_cavernosa,Orbicella_faveolata,Pseudodiploria_strigosa)
rm(muestras_florida,muestras_francia)
rm(ps_Chiarello_et_al2020)
rm(Siderastrea_siderea)