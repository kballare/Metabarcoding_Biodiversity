#Taxon-Specific Biodiversity - PITS
library(devtools)
library(decontam)
library(ranacapa)
library(plyr)
library(dplyr)

library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(plotly)
library(optparse)
library(microbiomeSeq) #R 4.0
library(stringr)

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

install_github("zdk123/SpiecEasi")

install.packages("remotes")
remotes::install_github("umerijaz/microbiomeSeq")


library(SpiecEasi)
oasis_color <- c("dodgerblue4", "goldenrod2")
color_blind1 <- c("#D55E00" , "#0072B2") 

#import phyloseq object
merged_PITS_all4 <- readRDS("Physeq_PITSmerged4_converted.rds")



#look for all unknown taxa to remove 
write.csv(merged_PITS_all4@tax_table, "PITS_tax.csv")  #in Excel sort by Domain, Phylum, class etc and scroll to end to find unknown, unknown, unknown (sometimes in  the middle without sorting)
#310 taxa, 1 unknown,  Domain=="unknown"

PITS_all4_nounknown <- subset_taxa(merged_PITS_all4, Domain == "Eukaryota") #remove NA assigned taxa
PITS_all4_nounknown #check that number of taxa is what you expect after NA removal
#309 taxa :) 

PITS_all4_nounknown <- prune_taxa(taxa_sums(PITS_all4_nounknown) > 10, PITS_all4_nounknown) #prune to 10 reads minumum


PITS_all4_nounknown_500 <- prune_samples(sample_sums(PITS_all4_nounknown) > 500, PITS_all4_nounknown)
PITS_all4_nounknown_500 #165 samples

PITS_all4_nounknown_500_sediments <- subset_samples(PITS_all4_nounknown_500, Substrate == "Sediment")
PITS_all4_nounknown_500_sediments <- prune_taxa(taxa_sums(PITS_all4_nounknown_500_sediments) > 10, PITS_all4_nounknown_500_sediments) 

PITS_all4_nounknown_500_sediments_noK0726_L7 <- subset_samples(PITS_all4_nounknown_500_sediments, sum.taxonomy != "K0726_L7")



#beta diversity
#use only sediment samples for beta diversity, 10 minumum reads/taxa, 500 minimum reads/sample

#conduct test - are communities different between sample date?
#ordination- jaccard- uses presence absence only
ord_PITS <- ordinate(PITS_all4_nounknown_500_sediments_noK0726_L7, "NMDS", "jaccard", k=3, trymax=100)
ord_PITS <- ordinate(PITS_all4_nounknown_500_sediments_noK0726_L7, "NMDS", "jaccard", k=3, trymax=500, previous.best=ord_PITS)

#all samples
# k=2 stress = 0.223
#k=3 stress = 0.129

#outlier removed
#k=2 no covergence
#k=3 stress=0.176


plot_ord_PITS_12 <- plot_ordination(PITS_all4_nounknown_500_sediments_noK0726_L7, ord_PITS, type="samples", color="Oasis", shape="Timepoint", axes= 1:2,label = NULL)
plot_ord_PITS_12 <- plot_ord_PITS_12+ scale_color_manual(values=c(oasis_color)) + labs(title="PITS, stress=0.176, k=3")#+ facet_wrap(~Processing_Batch)
plot_ord_PITS_12

plot_ord_PITS_23 <- plot_ordination(PITS_all4_nounknown_500_sediments_noK0726_L7, ord_PITS, type="samples", color="Oasis", shape="Timepoint", axes= 2:3,label = NULL)
plot_ord_PITS_23 <- plot_ord_PITS_13+ scale_color_manual(values=c(oasis_color)) + labs(title="PITS, stress=0.176, k=3")#+ facet_wrap(~Processing_Batch)
plot_ord_PITS_23


#faceted
plot_ord_PITS_facet <- plot_ordination(PITS_all4_nounknown_500_sediments_noK0726_L7, ord_PITS, type="samples", 
                                       color="Oasis", axes= 1:2,label = NULL) 
plot_ord_PITS_facet  <- plot_ord_PITS_facet  + geom_point(size=3)+theme_bw() + theme(aspect.ratio=1)  + scale_color_manual(labels= c("SP - Invaded","TPO - Pristine"), values=c(color_blind1)) + 
  facet_wrap(~Timepoint, labeller=labeller(Timepoint = c("2019_06" = "Jun 2019 \nPre-Eradication",
                                                         "2019_12" = "Dec 2019 \n5 Mos Post-Eradication",
                                                         "2020_06" = "Jun 2020 \n11 Mos Post-Eradication",
                                                         "2020_12" = "Dec 2020 \n17 Mos Post-Eradication"))) + theme(
                                                           plot.title = element_text(color="black", size=14),
                                                           axis.title.x = element_text(color="black", size=14, face="bold"),
                                                           axis.title.y = element_text(color="black", size=14, face="bold", margin=margin(t=20)),
                                                           #legend.title = element_text(size=12), legend.text = element_text(size=11),
                                                           legend.position="none",
                                                           axis.text = element_text(size=12))
                                                           #axis.text.x = element_text(angle=45,hjust=1))

                                                           #+ labs(title="PITS - k=3, stress=0.176")
                                                         

plot_ord_PITS_facet

#beta stats
dist_PITS <- distance(PITS_all4_nounknown_500_sediments_noK0726_L7, method="jaccard") #make distance matrix

metadata_PITS<- as(sample_data(PITS_all4_nounknown_500_sediments_noK0726_L7), "data.frame")

perm <- how(nperm=9999)
setBlocks(perm) <- with(metadata_FITS, Timepoint)
adonis_dist <- adonis2(dist_PITS~Oasis*Timepoint, data=metadata_PITS, permutations=perm) #run adonis (PERMANOVA)
adonis_dist

pairwiseAdonis <-pairwise.adonis2(dist_PITS~Oasis_Timepoint,data=metadata_PITS, strata="Timepoint", perm=9999)
pairwiseAdonis


#betadispersion - if significant, adonis result could be due to differences in dispersion
beta <- betadisper(dist_PITS, metadata_PITS$Timepoint)
permutest(beta)

TukeyHSD(beta) 


#alpha diversity
alpha.diversity <- estimate_richness(PITS_all4_nounknown_500_sediments, measures =c("Observed", "Chao1", "Shannon"))
data <- cbind(sample_data(PITS_all4_nounknown_500_sediments), alpha.diversity)
aov <- aov(Chao1~ Oasis*Timepoint, data=data)

summary(aov)

Tukey <- TukeyHSD(aov)
Tukey


p_alpha <- ggplot(data, x="Timepoint",  
                  y= "Chao1",
                  color="Oasis") 


p_alpha <- p_alpha + geom_boxplot(aes(x=Timepoint,
                                      y=Chao1, fill=Oasis), outlier.shape=1) + scale_fill_manual(values=c(oasis_color)) #+theme(axis.text.x = element_text(angle = 45, hjust=0.75))
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  ggtitle("PITS Biodiversity") #+ scale_y_continuous(limits=c(0,25))
p_alpha


#violin plot with mean points

p_alpha <- ggplot(data, aes(x=Timepoint,  
                            y= Chao1, fill=Oasis)) + geom_violin() + stat_summary(fun.data = "mean_cl_boot",
                                                                                  geom = "pointrange",
                                                                                  color = "white",  position=position_dodge(0.9))#+ geom_boxplot(width=0.1, color="white", alpha=0.2, position=position_dodge(0.9), outlier.shape=NA)


p_alpha <- p_alpha + scale_fill_manual(values=c(color_blind1), labels=c(
  "SP - Invaded", "TPO - Pristine")) + scale_x_discrete(labels=c("Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020")) 
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), axis.text.y=element_text(size=12)) + 
  labs(x="Sample Time Point",  y = "Chao1 Richness - PITS") +  #+ scale_y_continuous(limits=c(0,25))
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(size=12), legend.text = element_text(size=11)
  )
p_alpha
