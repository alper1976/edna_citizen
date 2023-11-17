
library(phyloseq)
phyloseq_object_fish = readRDS("~/Desktop/Master/Analyses/phyloseq_mifish_lca_right.rds")

sample_names(phyloseq_object_fish) = sapply(strsplit(sample_names(phyloseq_object_fish),"_"),'[',1) # corrects the rownames, so compatible with sample_data


#fjernet identity and coverage på slutten
c=read.table("/Users/lonekv/blast_results_mifish_final_lca_03_98_out",sep = "\t",fill=T) #lca output for mifish sample and scandifish output.
colnames(c) = c("Query","LCA_rank","LCA_taxon","Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Method") #creates correct colnames.



control <- prune_samples(sample_names(phyloseq_object_fish) != "N" & sample_names(phyloseq_object_fish) != "P", phyloseq_object_fish)

sample_names(control)


seal_fish_raw_pruned = prune_taxa(taxa_sums(control)>0, control) #remove taxa with 0 reads
sample_names(seal_fish_raw_pruned)
plot(estimate_richness(seal_fish_raw_pruned, measures="Observed"))

#plot(estimate_richness(seal_fish_raw_pruned, measures = "ACE"))

otu_table(seal_fish_raw_pruned)


tab <- otu_table(seal_fish_raw_pruned)
class(tab) <- "matrix" # as.matrix() will do nothing
## you get a warning here, but this is what we need to have
tab <- (tab) # transpose observations to rows
library(vegan)
png(filename="ASV_1.png", width=900, bg="white")
par(cex.axis = 1.5)
rarecurve(tab, step=100, lwd=2, ylab="ASVs", cex.lab = 1.5, label=F)
dev.off()

tab.data <- as.data.frame(tab)
seal_fish_raw_pruned = rarefy_even_depth(seal_fish_raw_pruned, rngseed=1, sample.size=2000, replace=F) # rarifies at 2000, replace smaller samples
sample_names(seal_fish_raw_pruned)
richness_fish =estimate_richness(seal_fish_raw_pruned)
richness_fish$diff_ace = richness_fish$Observed - richness_fish$ACE

mean(richness_fish$diff_ace,na.rm=T)
estimate_richness(seal_fish_raw_pruned)[,"ACE"]

rarecurve(t(otu_table(seal_fish_raw_pruned)), step=50, cex=0.5) #rarefaction curve




#subset_samples(phyloseq_object,Sample!=)


spec_accum_random_S001 <- specaccum(t(otu_table(seal_fish_raw_pruned)), method = "random", permutations = 100,
                                    conditioned =TRUE, gamma = "jack1")

mod_random_S001 <- fitspecaccum(spec_accum_random_S001, "arrh")

png(filename="ASV.png", width=900, bg="white")
plot(mod_random_S001, col="hotpink", xlab = "Number of sites", ylab = "ASVs", cex.lab = 1.5, cex.axis = 1.5)
boxplot(spec_accum_random_S001, col = "yellow", border = "blue", lty=1, cex=0.5, add= TRUE)
dev.off()



#sample_sums(seal_fish_raw_pruned)

# koden må gå gjennom

seal_fish_phylo=plot_bar(seal_fish_raw_pruned,fill="Species")
seal_fish_phylo_gg =seal_fish_phylo$data #convert to gg object

sppec <- unique(seal_fish_phylo_gg$Species)

no.id <- subset(seal_fish_phylo_gg, seal_fish_phylo_gg$LCA_taxon == "no identification")
no.id.un <- subset(no.id, no.id$OTU == unique(no.id$OTU))


sppec <- as.data.frame(sppec)
colnames(sppec) <- "species"
ble <- rbind(my_species, sppec)
ble <- as.data.frame(table(ble))
ble <- subset(ble,Freq==1)

klokk <- match_df(ble, my_species)