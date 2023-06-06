# Plots
# Colors
colores<-met.brewer("Cross",24)
colores 

#Alpha diversity plots

graph_alpha_lugar<-ggplot(data = metadatos_ambos,aes(x=Lugar,y=Shannon,fill=Lugar))+
  geom_boxplot()+
  theme_classic()+
  theme(legend.position = "none")+
  scale_fill_manual(values=c(colours[1],colours[7]))+
  theme(axis.title = element_text(face = "bold"))
graph_alpha_lugar
ggsave(filename="boxplot_alpha.png",graph_alpha_lugar) #saving plot

# Betadiversity biplot
biplot(PCoA2)
colours = met.brewer("Cassatt1",type="discrete")
s.class(PCoA2$li, fac = as.factor(metadatos_ambos$Lugar) , col=c(colours[1],colours[2]))
png("intento.png") #saving plot

# Heatmap core microbiome
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)

p <- plot_core(core.f, plot.type = "heatmap", 
               min.prevalence = 0.6, 
               prevalences = prevalences, 
               detections = detections) +
  xlab("Detection Threshold (Relative Abundance (%))")
print(p)
colores2<-colores<-met.brewer("Cross",type = "continuous")
core.fam <- aggregate_taxa(core, "Phylum")
any(taxa_names(core.fam) == "Unknown")
core.fam <- subset_taxa(core.fam, Phylum!="Unknown")
core.fam <- microbiome::transform(core.fam, "compositional")

p1 <- plot_core(core.fam, 
                plot.type = "heatmap", 
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5,
                colours = met.brewer("Cassatt1",type="continuous")) +
  xlab("Detection Threshold (Relative Abundance (%))")
p1 <- p1 + theme_bw() + ylab("ASVs") + xlab("Umbral de detección (Abundancia relativa (%)")
print(p1)
ggsave("core_prev.png",plot = last_plot())

#Venn diagrams per site

table(meta(phyloseq_completo)$Lugar, useNA = "always")
pseq.rel <- microbiome::transform(phyloseq_completo, "compositional")

lugares_lista<-unique(as.character(meta(phyloseq_completo)$Lugar))
print(lugares_lista)
list_core <- c()


for (n in lugares_lista){
  ps.sub <- subset_samples(pseq.rel, Lugar == n)
  core_m <- core_members(ps.sub, 
                         detection = 0.0001, prevalence = 50/100)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m)))
  list_core[[n]] <- core_m
  print(as.data.frame(core_m))
}
list_core
mycols <- c(Francia=colours[3], Florida=colours[6]) 
plot(venn(list_core),
     fills = mycols)


# Core microbiome table

intersect(list_core$Francia,list_core$Florida)

toda_asv<-as.data.frame(intento_merge2@tax_table)

toda_asv$asv_name<-rownames(toda_asv)
toda_asv<-subset(toda_asv,asv_name %in% intersect(list_core$Francia,list_core$Florida))

rownames(toda_asv)<-NULL
toda_asv<-toda_asv[,4:7]
toda_asv<-toda_asv %>%  relocate(asv_name,.before=Order)
colnames(toda_asv)<-c("ASV","Orden","Familia","Género")

tabla_final<-toda_asv %>% 
  gt() %>%
  tab_source_note(md("Tabla 1. ASV del core microbiome compartido entre ambos sitios")) %>% 
  opt_row_striping() %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything()))
gtsave(tabla_final,"tabla_core.rtf")

# Abundances

plot_abundance = function(physeq, ylabn = "",
                          Facet = "Phylum",
                          Color = "Phylum",filos){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  mphyseq <- subset(mphyseq,Phylum%in%filos)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "Lugar", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet,ncol = 3) + ylab(ylabn) +
    scale_y_log10()+
    theme_classic()+
    theme(legend.position = "none")
  
}
filos<-c("Acidobacteria","Actinobacteria","Bacteroidetes","Chlamydiae","Chlroflexi","Cyanobacteria","Epsilonbacteraeota","Euryarchaeota","Firmicutes","Fusobacteria","Planctomycetes","Proteobacteria","Thaumarchaeota","Verrucomicrobia")
plot_abundance(phyloseq_completo,"Abundancias",filos=filos)
ggsave("abundancias.png",plot = last_plot())


