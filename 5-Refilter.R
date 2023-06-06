# Re filtering data
#Removing AGAIN non-present ASV

#France

sumatoria_francia_2<-apply(ASV_francia_final,2,FUN = sum) #Sum ASV abundance in every sample
ASV_presentes_francia_2<-rbind(ASV_francia_final,sumatoria_francia_2) #Joint sum as the last row (4643 ASV)
sumatoria_francia_2<-as.data.frame(sumatoria_francia_2) ##Formating
sumatoria_francia_2<-subset(sumatoria_francia_2,sumatoria_francia_2==0) #dejar solo los que no están en ninguna muestra
asv_eliminar_francia_2<-row.names(sumatoria_francia_2) #Non-present ASV
ASV_presentes_francia_2 <- ASV_presentes_francia_2[,!names(ASV_presentes_francia_2) %in% c(asv_eliminar_francia_2)] #Only present ASV (2898)
ASV_francia_final<-ASV_presentes_francia_2 #Final name

# Florida

sumatoria_florida_2<-apply(ASV_florida_final,2,FUN = sum) 
ASV_presentes_florida_2<-rbind(ASV_florida_final,sumatoria_florida_2) 
sumatoria_florida_2<-as.data.frame(sumatoria_florida_2) 
sumatoria_florida_2<-subset(sumatoria_florida_2,sumatoria_florida_2==0) 
asv_eliminar_florida_2<-row.names(sumatoria_florida_2) 
ASV_presentes_florida_2 <- ASV_presentes_florida_2[,!names(ASV_presentes_florida_2) %in% c(asv_eliminar_florida_2)] 
ASV_florida_final<-ASV_presentes_florida_2

#remove unuseful objects

rm(ASV_presentes_florida_2,ASV_presentes_francia_2)
rm(sumatoria_florida_2,sumatoria_francia_2)
rm(asv_eliminar_florida_2)
rm(asv_eliminar_francia_2)

#Removing last rows (sums)

ASV_florida_final<-ASV_florida_final[-11,]
ASV_francia_final<-ASV_francia_final[-11,]



