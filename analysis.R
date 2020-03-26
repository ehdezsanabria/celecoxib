library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("SpiecEasi")
library("ggplot2")
library("dplyr")
library("Phenoflow")
library("dplyr")
library("grid")
library("vegan")
library("igraph")
library("DESeq2")
source("functions.R")

# Import data
df_phy <- import_mothur(mothur_shared_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pi....shared",
              mothur_constaxonomy_file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pi....taxonomy")
colnames(tax_table(df_phy)) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus")

# import metadata - make sure rownames are sample names
meta <- read.csv2("metadata.csv")

# Create sample data object
metadata <- sample_data(meta)
rownames(metadata) <- metadata$Sample
metadata$Donor <- gsub(pattern = "a", replacement = "", metadata$Donor)
metadata$Donor <- gsub(pattern = "b", replacement = "", metadata$Donor)
metadata$Donor <- factor(metadata$Donor)
metadata$merge_factor <- interaction(metadata$Donor, metadata$TPT, metadata$TRT)
metadata$Sample <- as.factor(metadata$Sample)

################################################################################
### Alpha diversity analysis
################################################################################

# Calculate diversity from rescaled OTU table
diversity_results <- Diversity_16S(scale_reads(df_phy), brea = FALSE, thresh = 1000, R = 100)
diversity_results <- data.frame(Sample = rownames(diversity_results), diversity_results)

# Merge diversity results with metadata
diversity_results <- left_join(diversity_results, metadata, by = "Sample")

# Plot diversity results
## Overall difference between treatments - D2 (Inverse Simpson)
p.div1_D2 <- ggplot(aes(x = TRT, y = D2, fill = TRT), data = diversity_results)+
  geom_jitter(shape = 21, size = 4, width = 0.2)+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=2,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-D2-OVERALL.png", width = 7, height = 6, res = 500, units = "in")
print(p.div1_D2)
dev.off()

## Overall difference between treatments - D0 (Observed Richness)
p.div1_D0 <- ggplot(aes(x = TRT, y = D0, fill = TRT), data = diversity_results)+
  geom_jitter(shape = 21, size = 4, width = 0.2)+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=2,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-D0-OVERALL.png", width = 7, height = 6, res = 500, units = "in")
print(p.div1_D0)
dev.off()

## Overall difference between treatments - Chao index (extrapolated Richness)
p.div1_Chao <- ggplot(aes(x = TRT, y = D0.chao, fill = TRT), data = diversity_results)+
  geom_jitter(shape = 21, size = 4, width = 0.2)+
  geom_boxplot(alpha = 0.7, outlier.shape = NA)+
  # facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=2,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)
  # geom_errorbar(aes(ymin = D0.chao- sd.D0.chao, ymax = D0.chao+sd.D0.chao), 
                # width= 0.015)

png(file = "./Figures_alphadiv/DIV-Chao-OVERALL.png", width = 7, height = 6, res = 500, units = "in")
print(p.div1_Chao)
dev.off()


## Overall difference between treatments for each donor and split for timepoints
p.div2_D2 <- ggplot(aes(x = TRT, y = D2, fill = TRT), data = diversity_results)+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.7)+
  facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=12,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-D2-DONOR.png", width = 14, height = 7, res = 500, units = "in")
print(p.div2_D2)
dev.off()

p.div2_D0 <- ggplot(aes(x = TRT, y = D0, fill = TRT), data = diversity_results)+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.7)+
  facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=12,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-D0-DONOR.png", width = 14, height = 7, res = 500, units = "in")
print(p.div2_D0)
dev.off()

p.div2_Chao <- ggplot(aes(x = TRT, y = D0.chao, fill = TRT), data = diversity_results)+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.7)+
  facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=12,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-Chao-DONOR.png", width = 14, height = 7, res = 500, units = "in")
print(p.div2_Chao)
dev.off()

# The same as the above but only for timepoint 16h
p.div2_D2_T16 <- diversity_results %>% dplyr::filter(TPT == "16h") %>% 
  ggplot(aes(x = TRT, y = D2, fill = TRT))+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.7)+
  facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=12,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-D2-DONOR-T16.png", width = 14, height = 5, res = 500, units = "in")
print(p.div2_D2_T16)
dev.off()

p.div2_D0_T16 <- diversity_results %>% dplyr::filter(TPT == "16h") %>%  
  ggplot(aes(x = TRT, y = D0, fill = TRT))+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.7)+
  facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=12,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-D0-DONOR-T16.png", width = 14, height = 5, res = 500, units = "in")
print(p.div2_D0_T16)
dev.off()

p.div2_Chao_T16 <- diversity_results %>% dplyr::filter(TPT == "16h") %>% 
  ggplot(aes(x = TRT, y = D0.chao, fill = TRT))+
  geom_point(shape = 21, size = 4)+
  geom_boxplot(alpha = 0.7)+
  facet_grid(TPT~Donor)+
  ylab("Alpha diversity")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=12,"Accent"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))+
  guides(fill=FALSE)

png(file = "./Figures_alphadiv/DIV-Chao-DONOR-T16.png", width = 14, height = 5, res = 500, units = "in")
print(p.div2_Chao_T16)
dev.off()

# Export diversity values with metadata
write.csv(file = "./outputfiles/diversity_CX.csv", diversity_results)

################################################################################
### Beta diversity analysis
################################################################################

# Add metadata to phyloseq object
sample_data(df_phy) <- metadata

# Rescale OTU table to account for library size differences
df_phy_scaled <- scale_reads(df_phy)

# Start making the PCoA
pcoa <- ordinate(
  physeq = df_phy_scaled, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)

# Now lets transform this to dataframe
pcoa.df <- data.frame(pcoa$vectors, sample_data(df_phy_scaled))

# And calculate the variance explained
var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

# Lets run an exploratory permanova to see suggested effect sizes.
# Similar variances across treatments: moe or less...
dist.seq <- vegan::vegdist(t(otu_table(df_phy_scaled)))
disper.seq <- vegan::betadisper(dist.seq, group = sample_data(df_phy_scaled)$Donor)
anova(disper.seq)
print(disper.seq)

# Permutations are constrained within each sampling year 
perm_results <- vegan::adonis(dist.seq ~ Donor + TRT + TPT, 
                              data = data.frame(sample_data(df_phy_scaled)))

# Add this information on plots
my_grob = grid::grobTree(textGrob(bquote(paste(r[Donor]^2 == 
                                           .(round(100 * perm_results$aov.tab[1, 5], 1)), 
                                         "%")), x = 0.7, y = 0.95, hjust = 0, gp = gpar(col = "black", 
                                                                                        fontsize = 14, fontface = "italic")))
my_grob2 = grid::grobTree(textGrob(bquote(paste(r[Treatment]^2 == 
                                            .(format(round(100 * perm_results$aov.tab[2, 5], 
                                                           1), nsmall = 1)), "%")), x = 0.7, y = 0.87, 
                             hjust = 0, gp = gpar(col = "black", fontsize = 14, 
                                                  fontface = "italic")))
my_grob3 = grobTree(textGrob(bquote(paste(r[Time]^2 == 
                                            .(round(100 * perm_results$aov.tab[3, 5], 1)), 
                                          "%")), x = 0.7, y = 0.79, hjust = 0, gp = gpar(col = "black", 
                                                                                         fontsize = 14, fontface = "italic")))
# Now we can plot the beta diversity plot
p_beta_donor <- ggplot(data = pcoa.df, aes(x=Axis.1, y=Axis.2, shape = TRT))+
  geom_point(alpha=0.7, size=7, aes(fill=Donor))+
  scale_shape_manual(values = c(21,24))+
  theme_bw()+
  scale_fill_brewer(palette = "Paired")+
  # scale_fill_manual(values = c('red', "blue"))+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)),
         shape  = guide_legend(title = "Treatment"))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)

png(file = "./Figures_betadiv/BETA-DONOR.png", width = 8, height = 6, res = 500, units = "in")
print(p_beta_donor)
dev.off()


p_beta_time <- ggplot(data = pcoa.df, aes(x=Axis.1, y=Axis.2, shape = TRT))+
  geom_point(alpha=0.7, size=7, aes(fill=TPT))+
  scale_shape_manual(values = c(21,24))+
  theme_bw()+
  scale_fill_brewer(palette = "Accent")+
  # scale_fill_manual(values = c('red', "blue"))+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), 
       y = paste0("PCoA axis 2 (",var[2], "%)"), 
       colour="")+  
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22), title = "Timepoint"),
         shape  = guide_legend(title = "Treatment"))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)

png(file = "./Figures_betadiv/BETA-TIME.png", width = 8, height = 6, res = 500, units = "in")
print(p_beta_time)
dev.off()

################################################################################
### Community composition overview figures
################################################################################

# Calculate relative abundances
physeq_rel_otu_reps <- transform_sample_counts(df_phy_scaled, function(x) 100*x/sum(x))

# Export relative OTU table
write.csv(file="./outputfiles/relative_otu-table_reps.csv", otu_table(physeq_rel_otu_reps))

# Pool at genus level
df_phy_genus = tax_glom(df_phy_scaled, "Genus")

# Calculate relative abundances
physeq_rel_reps <- transform_sample_counts(df_phy_genus, function(x) 100*x/sum(x))

# Export relative OTU table
write.csv(file="./outputfiles/relative_genus-table_reps.csv", otu_table(physeq_rel_reps))

# Merge replicate samples for visualization purposes
df_phy_genus <- merge_samples(df_phy_genus, "merge_factor")

# Repair replaced factors
sample_data(df_phy_genus)$TRT <- factor(sample_data(df_phy_genus)$TRT)
levels(sample_data(df_phy_genus)$TRT) <- c("CX", "PEG")

sample_data(df_phy_genus)$TPT <- factor(sample_data(df_phy_genus)$TPT)
levels(sample_data(df_phy_genus)$TPT) <- c("0h", "16h")

sample_data(df_phy_genus)$Donor <- factor(sample_data(df_phy_genus)$Donor)
levels(sample_data(df_phy_genus)$Donor) <- c("D1","D2","D3","D4","D5","D6","D7","D8")

# Calculate relative abundances
physeq_rel <- transform_sample_counts(df_phy_genus, function(x) 100*x/sum(x))

# Export relative OTU table
write.csv(file="./outputfiles/relative_genus-table.csv", otu_table(physeq_rel))

# Convert phyloseq to dataframe
df_physeq_rel <- psmelt(physeq_rel)

# Select top 12 OTUs
topN = 12
most_abundant_taxa = sort(taxa_sums(physeq_rel), TRUE)[1:topN]
print(most_abundant_taxa)
physeq_rel_top12 <- prune_taxa(names(most_abundant_taxa), physeq_rel)
physeq_rel_top12_df <- psmelt(physeq_rel_top12)

# Make barplot figure
p_barplot <- ggplot(aes(x = TRT, y = Abundance, fill = Genus), data = physeq_rel_top12_df)+
  geom_bar(stat='identity')+
  facet_grid(TPT~Donor)+
  ylab("Relative Abundance (%)")+ xlab("Treatment")+
  scale_fill_manual(values=brewer.pal(n=12,"Paired"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 14))

png(file = "./Figures_composition/BARPLOT-DONOR-TPT.png", width = 14, height = 8, res = 500, units = "in")
print(p_barplot)
dev.off()

write.csv(file = "./outputfiles/top12taxa_CX.csv", physeq_rel_top12_df)

################################################################################
### Network analysis with SparCC and SPIEC-EASI
################################################################################

# Lets do network stuff on non-normalized otu table (df_phy)
# Remove taxa not seen more than 3 times in at least 20% of the samples
df_phy_filtered <- filter_taxa(df_phy, function(x) sum(x > 20) > (0.5*length(x)), TRUE)

# Split up the data into the two treatments
df_phy_filtered_CX <- prune_samples(sample_data(df_phy_filtered)$TRT == "CX", df_phy_filtered)
df_phy_filtered_PEG <- prune_samples(sample_data(df_phy_filtered)$TRT == "PEG", df_phy_filtered)

# SPIEC-EASI network analysis for CEX
set.seed(777)
sp_easi_CX <- spiec.easi(df_phy_filtered_CX, method='mb', lambda.min.ratio=1e-2,
                           nlambda=20, icov.select.params=list(rep.num=20))
ig2.mb_CX <- adj2igraph(sp_easi_CX$refit,  vertex.attr=list(name=taxa_names(df_phy_filtered_CX)))

vsize_CX <- rowMeans(clr(t(otu_table(df_phy_filtered_CX)), 1))+10
Family_CX <- tax_table(df_phy_filtered_CX)[,"Family"]
vweights_CX <- summary(symBeta(getOptBeta(sp_easi_CX), mode='maxabs'))

png(file = "./Figures_network/NETWORK-CX-C20-A50.png", width = 8, height = 8, res = 500, units = "in")
plot_network_custom(ig2.mb_CX, df_phy_filtered_CX, type='taxa',
             line_weight = 2, hjust = +0.25, color = "Family", label_size = 3,
             label = "Genus")+
  scale_color_brewer(palette = "Paired")+
  geom_point(aes(size = vsize_CX, color = Family_CX))+
  guides(size = FALSE,
         color  = guide_legend(title = "Family", nrow = 3))+
  geom_text(aes_string(label = "Genus"), size = 3, 
            hjust = +0.25, na.rm = TRUE)+
  theme(legend.position="bottom", legend.text=element_text(size=12),
        text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

# SPIEC-EASI network analysis for PEG
vsize_PEG <- rowMeans(clr(t(otu_table(df_phy_filtered_PEG)), 1))+10
Family_PEG <- tax_table(df_phy_filtered_PEG)[,"Family"]
vweights_PEG <- summary(symBeta(getOptBeta(sp_easi_PEG), mode='maxabs'))

sp_easi_PEG <- spiec.easi(df_phy_filtered_PEG, method='mb', lambda.min.ratio=1e-2,
                         nlambda=20, icov.select.params=list(rep.num=50))
ig2.mb_PEG <- adj2igraph(sp_easi_PEG$refit,  vertex.attr=list(name=taxa_names(df_phy_filtered_PEG)))

png(file = "./Figures_network/NETWORK-PEG-C20-A50.png", width = 8, height = 8, res = 500, units = "in")
plot_network_custom(ig2.mb_PEG, df_phy_filtered_PEG, type='taxa',
             line_weight = 2, hjust = +0.25, color = "Family", label_size = 3,
             label = "Genus")+
  scale_color_brewer(palette = "Paired")+
  geom_point(aes(size = vsize_PEG, color = Family_PEG))+
  geom_text(aes_string(label = "Genus"), size = 3, 
            hjust = +0.25, na.rm = TRUE)+
  guides(size = FALSE,
         color  = guide_legend(title = "Family", nrow = 3))+
  theme(legend.position="bottom", legend.text=element_text(size=12),
        text = element_text(size = 12),
        plot.margin = unit(c(1,1,1,1), "cm"))
dev.off()

################################################################################
### Run DEseq2 between Donors (default pipeline)
################################################################################

# Analysis between donors for TPT = 0h
df_phy_T0 <- prune_samples(sample_data(df_phy)$TPT == "0h", df_phy)
diagdds <- phyloseq_to_deseq2(df_phy_T0, ~ Donor)
diagdds <- DESeq2::DESeq(diagdds, test = "Wald", fitType = "parametric")

# Summarize results
res = DESeq2::results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(df_phy_T0)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Make default result plot
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

p_deseq_donor_T0 <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + 
  geom_point(shape = 21, size=6, alpha = 0.7, color = "black") + 
  theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust=1))+
  scale_color_brewer(palette = "Paired")+
  ylab("Log2 fold change")

png(file = "./Figures_DESEQ/DESEQ-DONOR-P0.01-T0.png", width = 8, height = 6, res = 500, units = "in")
print(p_deseq_donor_T0)
dev.off()


# Analysis between donors for TPT = 16h
df_phy_T16 <- prune_samples(sample_data(df_phy)$TPT == "16h", df_phy)
diagdds <- phyloseq_to_deseq2(df_phy_T16, ~ Donor)
diagdds <- DESeq2::DESeq(diagdds, test = "Wald", fitType = "parametric")

# Summarize results
res = DESeq2::results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(df_phy_T16)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Make default result plot
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

p_deseq_donor_T16 <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + 
  geom_point(shape = 21, size=6, alpha = 0.7, color = "black") + 
  theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust=1))+
  scale_color_brewer(palette = "Paired")+
  ylab("Log2 fold change")

png(file = "./Figures_DESEQ/DESEQ-DONOR-P0.01-T16.png", width = 8, height = 6, res = 500, units = "in")
print(p_deseq_donor_T16)
dev.off()

# Perform DESEQ analysis for each donor separately (nothing is significant)
for(i in unique(sample_data(df_phy)$Donor)){
  # Analysis between donors for TPT = 16h
  df_phy_tmp <- prune_samples(sample_data(df_phy)$TPT == "16h", df_phy)
  df_phy_tmp <- prune_samples(sample_data(df_phy)$Donor == i, df_phy)
  diagdds <- phyloseq_to_deseq2(df_phy_tmp, ~ TRT)
  diagdds <- DESeq2::DESeq(diagdds, test = "Wald", fitType = "parametric")
  
  # Summarize results
  res = DESeq2::results(diagdds, cooksCutoff = FALSE)
  alpha = 0.01
  sigtab = res[which(res$padj < alpha), ]
  if(nrow(sigtab) > 0){
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(df_phy_tmp)[rownames(sigtab), ], "matrix"))
    head(sigtab)
    
    # Make default result plot
    theme_set(theme_bw())
    scale_fill_discrete <- function(palname = "Set1", ...) {
      scale_fill_brewer(palette = palname, ...)
    }
    # Phylum order
    x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
    x = sort(x, TRUE)
    sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
    
    # Genus order
    x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
    x = sort(x, TRUE)
    sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
    
    p_deseq_donor_tmp <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, fill=Phylum)) + 
      geom_point(shape = 21, size=6, alpha = 0.7, color = "black") + 
      theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust=1))+
      scale_color_brewer(palette = "Paired")+
      ylab("Log2 fold change")
    
    png(file = paste0("./Figures_DESEQ/DESEQ-DONOR-P0.01-T16-",i,".png",sep=""), width = 8, height = 6, res = 500, units = "in")
    print(p_deseq_donor_tmp)
    dev.off()
  }
}
