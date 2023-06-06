# Core microbiome

core<-as.data.frame(intento_merge@tax_table)
colnames(core)[1]<-"Domain"
core<-as.matrix(core)
core<-tax_table(core)
intento_merge2<-phyloseq(intento_merge@otu_table,core,intento_merge@sam_data)
intento_merge2
phyloseq_completo<-intento_merge2

core <- microbiome::transform(phyloseq_completo, "compositional")
taxa_names(core)

core.taxa.standard <- core_members(core, detection = 0.0001, prevalence = 50/100)
core.taxa.standard #core (60 taxa)
length(core.taxa.standard)
tax_table(core)[tax_table(core) == "k__"] <- NA
tax_table(core)[tax_table(core) == "p__"] <- NA
tax_table(core)[tax_table(core) == "c__"] <- NA
tax_table(core)[tax_table(core) == "o__"] <- NA
tax_table(core)[tax_table(core) == "f__"] <- NA
tax_table(core)[tax_table(core) == "g__"] <- NA
tax_table(core)[tax_table(core) == "s__"] <- NA

tax_table(core)[, colnames(tax_table(core))] <- gsub(tax_table(core)[, colnames(tax_table(core))],  pattern = "[a-z]__", replacement = "")

core.f <- microbiome::add_besthit(core)

taxa_names(core.f) 

core.taxa.standard <- core_members(core.f, detection = 0.0001, prevalence = 50/100)
length(core.taxa.standard)

core.core <- core(core.f, detection = 0.0001, prevalence = 50/100)
core.taxa <- taxa(core.core)
core.taxa
tax.mat <- tax_table(core.core)
tax.df <- as.data.frame(tax.mat)

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
