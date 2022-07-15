#CO1 biodiversity

library(devtools)

library(dplyr)
library(plyr)
library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(vegan)


library(microbiomeSeq) #R 4.0
library(pairwiseAdonis)
library(microbiome)
library(fantaxtic)

mixed <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
           "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
           "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
           "#8A7C64", "#599861")

oasis_color <- c("dodgerblue4", "goldenrod2")
oasis_blue <-c("yellowgreen", "dodgerblue4")
color_blind1 <- c("#D55E00" , "#0072B2") 

#####phyloseq obj creation and manipulation#####
#read in phyloseq object
merged_CO1_all4 <- readRDS("merged_CO1_all4.rds")


#look for all unknown taxa to remove 
merged_CO1_all4
#1703 taxa
write.csv(merged_CO1_all4@tax_table, "CO1_tax.csv")  #in Excel sort by Domain, Phylum, class etc and scroll to end to find unknown, unknown, unknown (sometimes in  the middle without sorting)
#3 unknown, 2 Domain=="unknown", 1 Domain==""

Domain_keep <- c("Eukaryota", "Archaea", "Bacteria")
merged_CO1_all4_nounknown <- subset_taxa(merged_CO1_all4, Domain %in% Domain_keep) #remove fully unknown assigned taxa
merged_CO1_all4_nounknown #check that number of taxa is what you expect after NA removal
#1700 taxa :)

#make pruned phyloseq obj
phyCO1 <- prune_taxa(taxa_sums(merged_CO1_all4_nounknown) > 10, merged_CO1_all4_nounknown) 
phyCO1_500 <- prune_samples(sample_sums(phyCO1) > 500, phyCO1)

#sediments only
phyCO1_500_sediments <- subset_samples(phyCO1_500, Substrate == "Sediment")
phyCO1_500_sediments <- prune_taxa(taxa_sums(phyCO1_500_sediments) > 10, phyCO1_500_sediments) 

#ordination
ord_CO1 <- ordinate(phyCO1_500_sediments, "NMDS", "jaccard", k=2, trymax=100)
#ord_FITS <- ordinate(phy16S_500_oases_noK0601_A1, "NMDS", "jaccard", k=2, trymax=200,  previous.best=ord_16S)

#k=2 500 sediments only stress=0.106

#plot
plot_ord_CO1_12 <- plot_ordination(phyCO1_500_sediments, ord_CO1, type="samples", color="Oasis", shape="Timepoint", axes= 1:2,label = NULL)
plot_ord_CO1_12 <- plot_ord_CO1_12+ scale_color_manual(values=c(oasis_color)) + labs(title="CO1, stress=0.10, k=2")#+ facet_wrap(~Processing_Batch)
plot_ord_CO1_12


#faceted
plot_ord_CO1_facet <- plot_ordination(phyCO1_500_sediments, ord_CO1, type="samples", 
                                       color="Oasis", axes= 1:2,label = NULL) 
plot_ord_CO1_facet  <- plot_ord_CO1_facet  + geom_point(size=3)+theme_bw()+ theme(aspect.ratio = 1) + scale_color_manual(labels= c("SP - Invaded","TPO - Pristine"), values=c(color_blind1)) + 
  facet_wrap(~Timepoint, labeller=labeller(Timepoint = c("2019_06" = "Jun 2019 \nPre-Eradication",
                                                         "2019_12" = "Dec 2019 \n5 Mos Post-Eradication",
                                                         "2020_06" = "Jun 2020 \n11 Mos Post-Eradication",
                                                         "2020_12" = "Dec 2020 \n17 Mos Post-Eradication"))) + theme(
                                                           plot.title = element_text(color="black", size=14),
                                                           axis.title.x = element_text(color="black", size=14, face="bold"),
                                                           axis.title.y = element_text(color="black", size=14, face="bold", margin=margin(t=20)),
                                                           #legend.title = element_text(size=12), legend.text = element_text(size=11),
                                                           legend.position = "none",
                                                           axis.text = element_text(size=12)
                                                           #+ labs(title="CO1 - k=2, stress=0.10") 
                                                         )

plot_ord_CO1_facet

#2019_06 for both SP and TPO clusters together, well apart from other samples 

dist_CO1 <- distance(phyCO1_500_sediments, method="jaccard") #make distance matrix

metadata_CO1<- as(sample_data(phyCO1_500_sediments), "data.frame")

perm <- how(nperm=9999)
setBlocks(perm) <- with(metadata_CO1, Timepoint)
adonis_dist <- adonis2(dist_CO1~Oasis*Timepoint, data=metadata_CO1, permutations=perm) #run adonis (PERMANOVA)
adonis_dist

pairwiseAdonis <-pairwise.adonis2(dist_CO1~Oasis_Timepoint,data=metadata_CO1, perm=9999, strata="Timepoint")#, strata="Oasis") #strata="Timepoint", perm=9999)
pairwiseAdonis

#TPO and SP still significantly different in 2019_06
#All oasis comparisions are significantly different with Timepoint=strata

#betadispersion - if significant, adonis result could be due to differences in dispersion
beta <- betadisper(dist_CO1, metadata_CO1$Timepoint)
permutest(beta)
TukeyHSD(beta)

#betatest is significant for some oasis_timepoints: SP_2019_06 vs TPO 2019_06, vs SP2020_06, vs SP 2020_12, vs TPO 2020_06
#significant for timepoint (most significant)
#not significant for just oasis

#alpha diversity
alpha.diversity <- estimate_richness(phyCO1_500_sediments, measures =c("Observed", "Chao1", "Shannon"))
data <- cbind(sample_data(phyCO1_500_sediments), alpha.diversity)
aov <- aov(Chao1~ Oasis*Timepoint, data=data)

summary(aov)

Tukey <- TukeyHSD(aov)
Tukey

p_alpha <- ggplot(data, x="Timepoint",  
                  y= "Chao1",
                  color="Oasis")


p_alpha <- p_alpha + geom_boxplot(aes(x=Timepoint,
                                      y=Chao1, fill=Oasis), outlier.shape=NA) #+ scale_color_manual(values=c("red", "blue", "green", "purple")) +theme(axis.text.x = element_text(angle = 45, hjust=0.75))
p_alpha <- p_alpha + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  ggtitle("CO1 Biodiversity") + scale_fill_manual(values=c(oasis_color))  #+ scale_y_continuous(limits=c(0,25))
p_alpha

#low biodiversity in 2019_06, swamped with proteobacteria
#biodiversity in CO1 seems to be declining over time? but maybe swamped with other things?
#2020_12 and 2020_06 similar in alpha diversity (not signficantly different)

#violin plot with mean points

p_alpha <- ggplot(data, aes(x=Timepoint,  
                            y= Chao1, fill=Oasis)) + geom_violin() + stat_summary(fun.data = "mean_cl_boot",
                                                                                  geom = "pointrange",
                                                                                  color = "white", position=position_dodge(0.9))#+ geom_boxplot(width=0.1, color="white", alpha=0.2, position=position_dodge(0.9), outlier.shape=NA)


p_alpha <- p_alpha + scale_fill_manual(values=c(color_blind1), labels=c(
  "SP - Invaded", "TPO - Pristine")) + scale_x_discrete(labels=c("Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020")) 
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), axis.text.y=element_text(size=12)) + 
  labs(x="Sample Time Point",  y = "Chao1 Richness - FITS") +  #+ scale_y_continuous(limits=c(0,25))
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(size=12), legend.text = element_text(size=11)
  )
p_alpha