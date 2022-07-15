#Arthropoda Biodiversity#

library(devtools)
library(decontam)
library(ranacapa)
library(dplyr)
library(plyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(plotly)
library(optparse)
library(RColorBrewer)
library(microbiome)

library(pairwiseAdonis)

color_blind1 <- c("#D55E00" , "#0072B2") 
color_blind2 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")


#import phyloseq objects
merged_18S_all4 <- readRDS("merged_18s_all4.rds")
merged_CO1_all4 <-readRDS("merged_CO1_all4.rds")

#merge
merged_18S_CO1_all4 <- merge_phyloseq(merged_18S_all4, merged_CO1_all4)

Phylum_keep <- c("Arthropoda")
Oasis_keep <- c("SP","TPO")

#filter to just Arthropoda
arthro <- subset_taxa(merged_18S_CO1_all4, Phylum %in% Phylum_keep)
arthro_oasis <- subset_samples(arthro, Oasis!="Tank")
arthro_oasis_sed <- subset_samples(arthro_oasis, Substrate=="Sediment")

arthro_prune <- prune_taxa(taxa_sums(arthro_oasis_sed) > 0, arthro_oasis_sed) 
arthro_oasis_prune <- prune_taxa(taxa_sums(arthro_oasis) > 0, arthro_oasis)

arthro_10tax_10 <- prune_samples(sample_sums(arthro_prune) > 10, arthro_prune)

arthro_10tax_200 <- prune_samples(sample_sums(arthro_10tax)>200, arthro_10tax)
arthro_10tax_500 <- prune_samples(sample_sums(arthro_10tax) > 500, arthro_10tax)


#beta diversity
#use only sediment samples for beta diversity, 10 minumum reads/taxa, 500 minimum reads/sample

#conduct test - are communities different between sample date?
#ordination- jaccard- uses presence absence only
ord_arthro <- ordinate(arthro_10tax_10, "NMDS", "jaccard", k=2, trymax=100)
ord_arthro <- ordinate(arthro_10tax_10, "NMDS", "jaccard", k=2, trymax=500, previous.best=ord_arthro)

plot_ord_arthro <- plot_ordination(arthro_10tax_10, ord_arthro, type="samples", 
                                   color="Oasis", axes= 1:2,label = NULL) 
plot_ord_arthro <- plot_ord_arthro + geom_point(size=3)+theme_bw() + labs(title="k=2, stress=0.232") + scale_color_manual(labels= c("SP - Invaded","TPO - Pristine"), values=c(color_blind1)) + 
  facet_wrap(~Timepoint, labeller=labeller(Timepoint = c("2019_06" = "Jun 2019 - Pre-Eradication",
                                                         "2019_12" = "Dec 2019 - 6 Mos Post-Eradication",
                                                         "2020_06" = "Jun 2020 - 12 Mos Post-Eradication",
                                                         "2020_12" = "Dec 2020 - 18 Mos Post-Eradication"))) + theme(
                                                           plot.title = element_text(color="black", size=14),
                                                           axis.title.x = element_text(color="black", size=14, face="bold"),
                                                           axis.title.y = element_text(color="black", size=14, face="bold", margin=margin(t=20)),
                                                           legend.title = element_text(size=12), legend.text = element_text(size=11),
                                                           axis.text = element_text(size=12)
                                                         )
plot_ord_arthro



#no repeats
dist_arthro <- distance(arthro_10tax_10, method="jaccard") #make distance matrix

metadata_arthro<- as(sample_data(arthro_10tax_10), "data.frame")

perm <- how(nperm=9999)
setBlocks(perm) <- with(metadata_arthro, Timepoint)
adonis_dist <- adonis2(dist_arthro~Oasis*Timepoint, data=metadata_arthro, permutations=perm) #run adonis (PERMANOVA)


adonis_dist


pairwiseAdonis <-pairwise.adonis2(dist_arthro~Oasis_Timepoint, metadata_arthro,strata="Timepoint", nperm=9999)

pairwiseAdonis

#betadispersion - if significant, adonis result could be due to differences in dispersion
beta <- betadisper(dist_arthro, metadata_arthro$Oasis_Timepoint)

permutest(beta)
TukeyHSD(beta)

#alpha diversity
alpha.diversity <- estimate_richness(arthro_10tax_10, measures =c("Observed", "Chao1", "Shannon"))
data <- cbind(sample_data(arthro_10tax_10), alpha.diversity)

aov <- aov(Chao1~ Oasis*Timepoint, data=data)

summary(aov)


Tukey <- TukeyHSD(aov)
Tukey


#violin plot with mean points

alpha.diversity <- estimate_richness(arthro_10tax_10, measures =c("Observed", "Chao1", "Shannon"))
data <- cbind(sample_data(arthro_10tax_10), alpha.diversity)

str(data)

data$Oasis <- as.factor(data$Oasis)


p_alpha <- ggplot(data, aes(x=Timepoint,  
                            y= Chao1, fill=Oasis)) + geom_violin() + stat_summary(fun.data = "mean_cl_boot",
                                                                                  geom = "pointrange",
                                                                                  color = "white", position=position_dodge(0.9)) + facet_grid(~Oasis) #+ geom_boxplot(width=0.1, color="white", alpha=0.2, position=position_dodge(0.9), outlier.shape=NA)


p_alpha <- p_alpha + scale_fill_manual(values=c(color_blind1), labels=c(
  "SP - Invaded", "TPO - Pristine")) + scale_x_discrete(labels=c("Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020")) 
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), axis.text.y=element_text(size=12)) + 
  labs(x="Sample Time Point",  y = "Chao1 Richness - Arthropoda") +  #+ scale_y_continuous(limits=c(0,25))
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(size=12), legend.text = element_text(size=11))
    #, labeller=labeller(Timepoint = c("2019_06" = "Jun 2019 - Pre-Eradication",
                                                            #"2019_12" = "Dec 2019 - 6 Mos Post-Eradication",
                                                            #"2020_06" = "Jun 2020 - 12 Mos Post-Eradication",
                                                            #"2020_12" = "Dec 2020 - 18 Mos Post-Eradication"))) + theme(
                                                            #  plot.title = element_text(color="black", size=14),
                                                             # axis.title.x = element_text(color="black", size=14, face="bold"),
                                                            #  axis.title.y = element_text(color="black", size=14, face="bold", margin=margin(t=20)),
                                                           #   legend.title = element_text(size=12), legend.text = element_text(size=11),
                                                          #    axis.text = element_text(size=12)
                                                          #  )
p_alpha


#barplot
library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)
library(hrbthemes)
library(gcookbook)

install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)

sample_data(arthro_10tax_10)$Oasis <- as.factor(sample_data(arthro_10tax_10)$Oasis)
sample_data(arthro_10tax_10)$Timepoint <- as.factor(sample_data(arthro_10tax_10)$Timepoint)
sample_data(arthro_10tax_10)$Oasis_Timepoint <- as.factor(sample_data(arthro_10tax_10)$Oasis_Timepoint)


arthro_10tax_10_class <- arthro_10tax_10 %>% aggregate_taxa(level="Class") %>% microbiome::transform(transform="compositional")
arthro_10tax_10_order <- arthro_10tax_10 %>% aggregate_taxa(level="Order") %>% microbiome::transform(transform="compositional")
arthro_10tax_10_family <- arthro_10tax_10 %>% aggregate_taxa(level="Family") %>% microbiome::transform(transform="compositional")
arthro_10tax_10_genus <- arthro_10tax_10 %>% aggregate_taxa(level="Genus") %>% microbiome::transform(transform="compositional")



barplot <- plot_composition(arthro_10tax_10_class, average_by="Oasis_Timepoint", transform="compositional") + 
  scale_fill_brewer("Class", palette="Paired") + 
  labs(x = "Sample Time Point", y = "Relative Abundance (%)",
       title = "Arthropod Relative Abundance") + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), 
        axis.text.y=element_text(size=12), panel.border=element_blank(),
          panel.grid = element_blank(), axis.line=element_line(colour="black")) + geom_hline(yintercept=0) +
   scale_x_discrete(labels=c("Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020", "Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020")) 
  
barplot 
layer_data(barplot, 1)

#look into differences of specific taxa between timepoints
install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)

mypalette<-brewer.pal(12, "Set3")

pn <- plot_taxa_boxplot(arthro_10tax_10,
                        taxonomic.level = "Class",
                        top.otu = 11, 
                        group = "Oasis_Timepoint",
                        add.violin= FALSE,
                        title = "Top Class", 
                        keep.other = FALSE,
                        #group.order = c("H","CRC","nonCRC"),
                        group.colors=mypalette,
                        dot.size = 1)
pn

plot_alpha <-plot_alpha_diversities(arthro_10tax_10,
                                    type = "dominance",
                                    index.val = "all",
                                    plot.type = "violin",
                                    variableA = "Oasis_Timepoint",
                                    palette = mypalette)

plot_alpha

p.m <- plot_diversity_stats(arthro_10tax_10, group = "Oasis_Timepoint", 
                            index = "diversity_absolute", 
                            #group.order = c("H", "CRC", "nonCRC"), 
                            group.colors = mypalette,
                            label.format="p.format",
                            stats = TRUE)
p.m + ylab("Alpha Diversity") + xlab("")

top_taxa <- top_taxa(arthro_10tax_10_class, 4)
top_taxa_order <- top_taxa(arthro_10tax_10_order, 4)
top_taxa_family <- top_taxa(arthro_10tax_10_family, 4)
top_taxa_genus <- top_taxa(arthro_10tax_10_genus, 4)

taxa_AC <- c("Arachnida", "Branchiopoda", "Chilopoda", "Collembola")
taxa_DM <- c("Diplopoda", "Hexanauplia", "Insecta", "Malacostraca")
taxa_OU <- c('Ostracoda', "Protura", "unknown")
arthro_oasis <- subset_samples(arthro, Oasis!="Tank")
arthro_10tax_10_class_SP <- subset_samples(arthro_10tax_10_class, Oasis=="SP")
arthro_10tax_10_class_TPO <- subset_samples(arthro_10tax_10_class, Oasis=="TPO")
mytaxa1 <- c("Arachnida", "Hexanauplia") 
mytaxa2 <- c("Insecta", "Ostracoda")

p_taxa <- plot_listed_taxa(arthro_10tax_10_genus, top_taxa_genus, 
                      group= "Oasis_Timepoint",
                      #group.order = c("H","CRC","nonCRC"),
                      group.colors = mypalette,
                      add.violin = FALSE,
                      violin.opacity = 0.3,
                      dot.opacity = 0.25,
                      box.opacity = 0.25,
                      panel.arrange = "wrap")
p_taxa <- p_taxa + theme(axis.text.x = element_text(angle = 45, size=10, hjust=1, vjust=1)) +
             ylab("Relative Abundance per Sample") + scale_y_continuous(labels=scales::percent)
                    

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
p_taxa
library(ggpubr)
comps <- make_pairs(sample_data(arthro_10tax_10_genus)$Oasis_Timepoint)
comps_index <- comps[c(1:3, 23:25)]
p_taxa <- p_taxa + stat_compare_means(
  comparisons = comps_index,
  label = "p.format",
  tip.length = 0.05,
  vjust = 1.3,
  method = "wilcox.test",
  symnum.args = symnum.args)
p_taxa
warnings()
