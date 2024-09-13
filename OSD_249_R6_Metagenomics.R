df = read.csv("/Users/oluwamayowaakinsuyi/Downloads/file2.csv")
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



###extract tax table 
cs = c(1,2)
tax1 = df_new2[,cs]
View(tax)
dim(tax1)
tax1 = as.matrix(tax1)


######################### Extract otu table 
cv = c(1,2)
otu1 = df_new2[,-cv]
View(tax)
otu1 = as.matrix(otu)
##otus = t(otu)
View(otu)
#change NAs to Zeros
otu1[is.na(otu1)] <-0
#View(otu)
#dim(otu)

library(phyloseq)
metadata = read.csv("/Users/oluwamayowaakinsuyi/Desktop/Nasa_paper/Datam.csv")
View(metadata)
######make phyloseq 
ps1 <- phyloseq(otu_table(otu1, taxa_are_rows= T),
               tax_table(tax1))
ps1

View(metadata)
row.names(metadata) = metadata$X.Sample_ID
sample_data(ps1) = sample_data(metadata) 
ps1

library(phyloseq)
 View(tax_table(ps_all_r))



######
####View the reads found in each sample
library(phyloseq) 
all.reads1 = data.table(as(sample_data(ps), "data.frame"),
                       TotalReads = sample_sums(ps), keep.rownames = TRUE)
View(all.reads)
####

#####Rarefy the dataset
set.seed(23)
ps_all_r1= rarefy_even_depth(ps, sample.size = 1604378, replace = FALSE)
data.frame = as(sample_data(ps_all_r1), "data.frame")
plyr::count(data.frame, "Group")
ps_all_r1
#######
ps_all_r1



library(ggplot2)
library(ggpubr)

######Alpha diversity
########

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
  scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
  scale_fill_brewer(palette = palname, ...)
}
alpha_means <- "Shannon"
b<- plot_richness(ps_all_r, x = "Group", measures = alpha_means, color = "Group") +
  geom_boxplot(width = 0.5, alpha = 0.5) +  # Adjust width and add transparency to the boxplot
  geom_point(position = position_jitterdodge(), alpha = 0.5) +  # Add jittering to the points for better visibility
  stat_compare_means() +
  xlab("") +
  ylab("Shannon Diversity") +  # Adding label for Shannon diversity on the y-axis
  theme(axis.text.x = element_blank(),
        axis.title.y = element_blank(),  # Remove y-axis label
        strip.text = element_blank(),  # Remove facet labels
        strip.background = element_blank())  # Remove facet background
b





ps.gen1 <- phyloseq::tax_glom(ps_all_r1, "Species", NArm = TRUE)
View(tax_table(ps.gen))

ps.gen1
############Deseq2
dds = phyloseq_to_deseq2(ps.gen1, ~Group)

#Perform geometric mean normalization
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dds), 1, gm_mean)
psfdds2 = estimateSizeFactors(dds, geoMeans = geoMeans)
psfdds2

dds2 = DESeq(dds, fitType="local")
dds2

##Investigate test results table
res = results(dds2, cooksCutoff = FALSE)
res
alpha = 0.9
sigtab = res[which(res$pvalue < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.gen1)[rownames(sigtab), ], "matrix"))
View(sigtab)

######### Agglomerate 
ps.gen1 <- phyloseq::tax_glom(ps_all_r1, "Species", NArm = TRUE)
View(tax_table(ps.gen1))



########
 # Extract abundance data from phyloseq object
 abundance_data1 <- as.data.frame(otu_table(ps.gen1))
dim(abundance_data1)
 taxa_data1 <- as.data.frame(tax_table(ps.gen1))
dim(taxa_data1)
 # Initialize an empty list to store results
 wilcox_results_list <- list()
 
 # Perform pairwise Wilcoxon test for each taxon
 wilcox_results_list <- apply(abundance_data1, 1, function(abundance_data1_row) {
   wilcox_result <- wilcox.test(abundance_data1_row ~ groups, Paired = True)
   return(wilcox_result$p.value)
 })
 wilcox_results_df1 <- data.frame(Taxon = rownames(abundance_data), P_Value = wilcox_results_list)
 View(wilcox_results_df)


 # Merge with taxa_data
 wilcox_results_df1 <- cbind(wilcox_results_df1, taxa_data)
 # View the results
 View(wilcox_results_df)
 



##############################################
ps_all_r_clr1 = microbiome::transform(ps.gen, 'clr')
 ########
 
 ps_alla = subset_taxa(ps_all_r_clr1, Species =="Lactococcus cremoris")
 ps_alla 
 
 
 ps_all1 = subset_taxa(ps_all_r_clr, Species =="Campylobacter hyointestinalis")
 ps_all1
 ps_all3 = subset_taxa(ps_all_r_clr, Species =="Enterobacter cloacae")
 ps_all5 = subset_taxa(ps_all_r_clr1, Species =="Klebsiella aerogenes")
 ps_all6 = subset_taxa(ps_all_r_clr1, Species =="Proteus hauseri")
 ps_all6
 
 
 ##########
 theme_set(theme_bw())
 pal = "Set1"
 scale_colour_discrete <- function(palname = pal, ...) {
   scale_colour_brewer(palette = palname, ...)
 }
 scale_fill_discrete <- function(palname = pal, ...) {
   scale_fill_brewer(palette = palname, ...)
 }

 library(microbiome)
 
library(phyloseq)
library(ggplot2)
library(ggpubr)

 
 

 

 library(ggpubr)
 #####
 ps_alla = subset_taxa(ps_all_r_clr, Species =="Gardnerella vaginalis")
 ps_alla 
 ps_alla <- phyloseq::psmelt(ps_alla) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),
         axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_alla
 
 

 ps_allb = subset_taxa(ps_all_r_clr, Species =="Bifidobacterium longum")
 ps_allb 
 ps_allb <- phyloseq::psmelt(ps_allb) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),
         axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_allb
 
 ##################################
 ps_allc = subset_taxa(ps_all_r_clr, Species =="Collinsella stercoris")
 ps_allc 
 ps_allc <- phyloseq::psmelt(ps_allc) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),
         axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_allc
 
 ##############
 ##################################
 ps_alld = subset_taxa(ps_all_r_clr, Species =="Collinsella sp. zg1085")
 ps_alld 
 ps_alld <- phyloseq::psmelt(ps_alld) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),
         axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_alld
 
 
 
 ps_alle = subset_taxa(ps_all_r_clr, Species =="Bacteroides luhongzhouii")
 ps_alle 
 ps_alle <- phyloseq::psmelt(ps_alle) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_alle
 
 
 
 ps_allf = subset_taxa(ps_all_r_clr, Species =="Roseburia sp. NSJ-69")
 ps_allf 
 ps_allf <- phyloseq::psmelt(ps_allf) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(), axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_allf
 
 
 
 tiff("Metaall.tiff", units = "in", width = 18, height = 7.5, res = 300)
 Fig4=ggarrange( ps_alle,ps_allb,ps_allc,ps_alld,ps_allf,ps_alla,
                labels = c("A","B","C","D","E","F"),
                common.legend = T,legend = "bottom",
                font.label = list(size = 10),
                nrow = 1)
 Fig4
 dev.off()
 View(sample_data(ps_all_r))
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 library(ggplot2)#######################
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
 
 # Combine the effect size results into a single data frame
 effect_sizes_summary <- bind_rows(effect_sizes, .id = "OTU")
 
 # Print the summary data frame
View(effect_sizes_summary)
 
 
p_effect= merge (effect_sizes_summary,wilcox_results_df)
 
dim(wilcox_results_df)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 View(taxa_data)
 
 
 
 
 
 
 ############
 ps_all1 = subset_taxa(ps_all_r_clr1, Species =="Mediterraneibacter")
 ps_all1 
 ps_all1 <- phyloseq::psmelt(ps_all1) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_all1
 
 
 ps_all_r_clr
 

 #####
 ps_all3 = subset_taxa(ps_all_r_clr1, Species =="Enterobacter cloacae")
 ps_all3 <- phyloseq::psmelt(ps_all3) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),   
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_all3 
 
 ###########
 ps_all2 = subset_taxa(ps_all_r_clr1, Species =="Klebsiella aerogenes")
 ps_all2 <- phyloseq::psmelt(ps_all2) %>%
   ggplot(data = ., aes(x = Group, y = Abundance, colour = Group)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   facet_wrap(~ Species) +
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_all2
 
 ps_all5 = subset_taxa(ps_all_r_clr1, Species =="Proteus hauseri")
 ps_all5 <- phyloseq::psmelt(ps_all5) %>%
   ggplot(data = ., aes(x = Groups, y = BA, colour = Groups)) +
   geom_boxplot(outlier.colour = NA, alpha = 0.5) +  # Adjust boxplot transparency
   labs(x = "", y = "CLR Abundance") +
   geom_jitter(aes(color = Group), height = 0, width = 0.2, alpha = 0.5) +  # Adjust point transparency
   theme(axis.text.x = element_blank(),axis.title.y = element_blank(),  # Remove y-axis label
         strip.text = element_text(size = 14, face = "italic")) +  # Italicize facet labels only
   stat_compare_means()
 
 ps_all5
 
 
 tiff("Meta2.tiff", units = "in", width = 11, height = 7.5, res = 300)
 Fig4=ggarrange(ps_all1,ps_all3,ps_all5,ps_all6,
                labels = c("A", "B","C","D"),
                common.legend = T,legend = "bottom",
                font.label = list(size = 10),
                nrow = 1)
 Fig4
 dev.off()
 View(sample_data(ps_all_r))
 
 
 
tiff("Meta1.tiff", units = "in", width = 9, height = 7.5, res = 300)
Fig4=ggarrange(p,ps_all,ps_alla,
                labels = c("A", "B","C"),
                common.legend = T,legend = "bottom")
 
Fig4
dev.off() 











#################
tiff("Meta1.tiff", units = "in", width = 9, height = 6, res = 300)
Figall=ggarrange(ps_d3,ps_d2,ps_d1,
                 labels = c("A", "B","C"),
                 common.legend = T,legend = "bottom", nrow = 1)

Figall
dev.off()

View(tax_table(ps_all_r_clr))
library(phyloseq)

######################effect size
library(rstatix)#######################

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


p_effect= merge (effect_sizes_summary,wilcox_results_df)

dim(wilcox_results_df)
 