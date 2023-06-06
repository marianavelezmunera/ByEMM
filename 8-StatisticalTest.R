# Statistical testing

#Full metadata, including place and shannon index
francia_meta<-metadatos_francia_finales[,7:9]
francia_meta$lugar<-c(rep("Francia",nrow(francia_meta)))
colnames(francia_meta)<-c("Especies","Muestra","Shannon","Lugar")

florida_meta<-metadatos_florida_finales[,c(1,5,6)]
florida_meta<-florida_meta %>%
  relocate(Especies,.before = Sample.ID)
florida_meta$lugar<-c(rep("Florida",nrow(florida_meta)))
colnames(florida_meta)<-c("Especies","Muestra","Shannon","Lugar")
metadatos_ambos<-rbind(francia_meta,florida_meta) #Joined info

## Kolmogorov-Smirnov test for diversities

test_sitio<-metadatos_ambos[,3:4]
spl <- split(metadatos_ambos$Shannon, metadatos_ambos$Lugar)
pv <- ks.test(spl$Francia, spl$Florida) #p value 

## ANOSIM

ANOSIM<-anosim(matrixdist, metadatos_ambos$Lugar, distance="bray",permutations=9999)
ANOSIM