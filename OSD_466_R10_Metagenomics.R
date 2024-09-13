df = read.csv("/Users/oluwamayowaakinsuyi/Downloads/metag1.csv")
View(df)
######### clean the data
# Pick the last two words from rows with '<' or '>'
df$Species <- sub(".*>", "", df$Species)

View(df)
dim(df)

##### remove viruses 
df2 <- df[!grepl("phage|virus", df$Species, ignore.case = TRUE), ]
View(df2)
dim(df2)


# Assuming your dataframe is named 'your_dataframe' remove every rank thats not genus 
df_new1 <- subset(df2, !grepl("G", taxRank, ignore.case = TRUE))
View(df_new1)
df_new1 <- subset(df_new1, !grepl("F", taxRank, ignore.case = TRUE))
View(df_new1)
df_new1 <- subset(df_new1, !grepl("U", taxRank, ignore.case = TRUE))
View(df_new1)
df_new1 <- subset(df_new1, !grepl("C", taxRank, ignore.case = TRUE))
View(df_new1)
df_new1 <- subset(df_new1, !grepl("P", taxRank, ignore.case = TRUE))
View(df_new1)
df_new1 <- subset(df_new1, !grepl("O", taxRank, ignore.case = TRUE))
View(df_new1)
df_new1 <- subset(df_new1, !grepl("D", taxRank, ignore.case = TRUE))
View(df_new1)

df_new2 <- subset(df_new1, taxRank != "-")
View(df_new2)
dim(df_new2)

###########
########### set new unique row names 
row_names <- paste0("sets", 1:nrow(df_new2))
rownames(df_new2) <- row_names
View(df_new2)

##### Make the phyloseq object
###extract tax table 
cs = c(1,2)
tax = df_new2[,cs]
View(tax)
dim(tax)
tax = as.matrix(tax)


#########################
cv = c(1,2)
otu = df_new2[,-cv]
View(otu)
otu = as.matrix(otu)
##otus = t(otu)
View(otu)
#change NAs to Zeros
otu[is.na(otu)] <-0
#View(otu)
#dim(otu)

library(phyloseq)
metadata = read.csv("/Users/oluwamayowaakinsuyi/Desktop/Nasa_paper/datanew.csv")
View(metadata)
######make phyloseq 
ps <- phyloseq(otu_table(otu, taxa_are_rows= T),
               tax_table(tax))
ps

View(metadata)
row.names(metadata) = metadata$X.Sample_id
sample_data(ps) = sample_data(metadata) 
ps





######
####View the reads found in each sample
library(data.table)
all.reads = data.table(as(sample_data(ps), "data.frame"),
                       TotalReads = sample_sums(ps), keep.rownames = TRUE)
View(all.reads)
####

#####Rarefy the dataset
set.seed(23)
ps_all_r= rarefy_even_depth(ps, sample.size = 1777084, replace = FALSE)
data.frame = as(sample_data(ps_all_r), "data.frame")
plyr::count(data.frame, "Group")
ps_all_r
#######
ps_all_r

ps.gen <- phyloseq::tax_glom(ps_all_r, "Species", NArm = TRUE)
View(tax_table(ps.gen))




library(ggplot2)
library(ggpubr)
######Alpha diversity
alpha_means1 <- "Shannon"
p <- plot_richness(ps_all_r, x = "Group", measures = alpha_means1, color = "Group") +
  geom_boxplot() +
  stat_compare_means() +
  xlab("") +
  ylab("Shannon Diversity") +  # Adding label for Shannon diversity on the y-axis
  theme(axis.text.x = element_blank(),
        strip.text = element_blank(),  # Remove facet labels
        strip.background = element_blank())  # Remove facet background
p




########Merge the alpha diversity  figure 
tiff("Metafig1.tiff", units = "in", width = 8, height = 7.5, res = 300)
Fig4=ggarrange(p,p1,
               labels = c("A", "B"),
               common.legend = T,legend = "bottom")

Fig4
dev.off() 


########wilcox
# Extract abundance data from phyloseq object
########
# Initialize an empty list to store results
wilcox_results_list <- list()

# Loop through taxa and perform Wilcoxon test
for (taxon_row in rownames(abundance_data)) {
  # Subset abundance data for the current taxon
  taxon_abundance <- abundance_data[taxon_row, ]
  
  # Extract abundance values for each group
  control_abundance <- taxon_abundance[groups == "control"]
  flight_abundance <- taxon_abundance[groups == "flight"]
  
  # Apply Wilcoxon rank sum test
  wilcox_result <- wilcox.test(control_abundance, flight_abundance)
  
  # Store the p-value in the list
  wilcox_results_list[[taxon_row]] <- wilcox_result$p.value
}

# View results
View(wilcox_results_list)

# Combine the results into a data frame
results_df <- data.frame(taxon = names(wilcox_results_list),
                         p_value = unlist(wilcox_results_list))

# View the results
class(results_df)
taxa = as.data.frame(taxa)
View(taxa)
###################
# Extract the species column from the tax table
species_column <- taxa$Species

# Add the species column to results_df
results_df$Species <- species_column

# View the updated results_df
View(results_df)



############## Deseq2
library(microbiome)
# filter sparse features, with > 90% zeros
ps.tax <- prune_taxa(rowSums(otu_table(ps.gen) == 0) < ncol(otu_table(ps.gen)) * 0.9, ps.gen)
ps_ds = phyloseq_to_deseq2(ps.tax, ~ Group)
library(phyloseq)
ds <- estimateSizeFactors(ps_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res_df = as.data.frame(res)
View(res_df)
dim(res_df)
tax = tax_table(ps.gen) 
View(tax)
res = res[order(res$padj, na.last=NA), ]
taxa_sig = rownames(res[1:20, ]) # select bottom 20 with lowest p.adj values
ps.taxa.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps.taxa.rel.sig <- prune_taxa(taxa_sig, ps.taxa.rel)



#########
########################3
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}
############
ps_all_r_clr = microbiome::transform(ps.gen, 'clr')
ps_alld1 = subset_taxa(ps_all_r_clr, Species == "Akkermansia glycaniphila")
ps_alld1 <- phyloseq::psmelt(ps_alld1) %>%
  ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
  facet_wrap(~ Species) +
  labs(x = "", y = "CLR Abundance") +
  geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
  theme(axis.text.x = element_blank(),
        #axis.title.y = element_blank(),  # Remove y-axis label
        strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
  stat_compare_means()
ps_alld1


ps_all_r_clr = microbiome::transform(ps.gen, 'clr')
ps_d3 = subset_taxa(ps_all_r_clr, Species == "Serratia liquefaciens")
ps_d3 <- phyloseq::psmelt(ps_d3) %>%
  ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
  facet_wrap(~ Species) +
  labs(x = "", y = "CLR Abundance") +
  geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),  # Remove y-axis label
        strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
  stat_compare_means()
ps_d3

ps_all_r_clr
#######
ps_d4 = subset_taxa(ps_all_r_clr, Species == "Citrobacter rodentium")
ps_d4 <- phyloseq::psmelt(ps_d4) %>%
  ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
  facet_wrap(~ Species) +
  labs(x = "", y = "CLR Abundance") +
  geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
  theme(axis.text.x = element_blank(),
        #axis.title.y = element_blank(),  # Remove y-axis label
        strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
  stat_compare_means()
ps_d4







library(ggpubr)
ps_d4 = phyloseq::psmelt(ps_d4) %>%
  ggplot(data = ., aes(x = Group, y = Abundance, colour = Group )) +
  geom_boxplot(outlier.colour  = "NA") +  
  facet_wrap(~ Species) + 
  labs(x = "", y = "CLR Abundance")  + 
  geom_jitter(aes(color = Group), height = 0, width = .2)  + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank()) +  # This line removes the Y-axis label
  stat_compare_means()
ps_d4



#############
tiff("Meta3f.tiff", units = "in", width = 9, height = 7.5, res = 300)
Fig4=ggarrange(ps_alld1,ps_all,ps_alla,
               labels = c("A","B","C"),
               common.legend = T,legend = "bottom",
               font.label = list(size = 10),
               nrow = 1)
Fig4
dev.off()


library(vegan)


tiff("Meta3f.tiff", units = "in", width = 9, height = 7.5, res = 300)
library(ggpubr)
Fig4 <- ggarrange(
  ps_d4, ps_d3, ps_all1, ps_all3, ps_all2, ps_all5,
  labels = c("A", "B", "C", "D", "E", "F", "G"),
  common.legend = TRUE,
  legend = "bottom",
  font.label = list(size = 10),
  ncol = 3, nrow = 3
)

Fig4


library(vegan)
library(ggpubr)
library(ggpubr)

tiff("Meta4f.tiff", units = "in", width = 9, height = 7.5, res = 300)
Fig4 <- ggarrange(
  ps_d4, ps_d3, ps_all1, ps_all3, ps_all2, ps_all5,
  labels = c("A", "B", "C", "D", "E", "F", "G"),
  common.legend = TRUE,
  legend = "bottom",
  font.label = list(size = 10),
  ncol = 3, nrow = 2,
  widths = c(3, 3, 3, 3),  # Adjust the widths of the columns
  heights = c(3, 3)  # Adjust the heights of the rows
)

Fig4

dev.off()

library(rstatix)
####################
# Convert the phyloseq object to a data frame
x <- as.data.frame(psmelt(ps_all_r_clr))
View(x)
# Initialize an empty list to store results
effect_sizes <- list()

# Loop through each unique OTU
for (otu in unique(x$OTU)) {
  # Subset the data for the current OTU
  x_subset <- subset(x, OTU == otu)
  
  # Calculate effect size using wilcox_effsize
  effect_size <- x_subset %>%
    wilcox_effsize(Abundance ~ Group)
  
  # Add the effect size results to the list
  effect_sizes[[otu]] <- effect_size
}

library(dplyr)
# Combine the effect size results into a single data frame
effect_sizes_summary <- bind_rows(effect_sizes, .id = "OTU")

# Print the summary data frame
View(effect_sizes_summary)


View(tax_table(ps_all_r_clr))



p_effect= merge (effect_sizes_summary,wilcox_results_df)

dim(wilcox_results_df)

