#FITS biodiversity

######load packages######
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
merged_FITS_all4 <- readRDS("merged_FITS_all4.rds")

#look for all unknown taxa to remove 
write.csv(merged_FITS_all4@tax_table, "FITS_tax.csv")  #in Excel sort by Domain, Phylum, class etc and scroll to end to find unknown, unknown, unknown (sometimes in  the middle without sorting)
#1520 taxa, 3 unknown, all has Domain=="unknown"

#filter unknowns
merged_FITS_all4_nounknowns <- subset_taxa(merged_FITS_all4, Domain != "unknown") #remove fully unknown assigned taxa
merged_FITS_all4_nounknowns #check that number of taxa is what you expect after unknown/NA removal
#1517 taxa :)

#make pruned phyloseq obj
phyFITS <- prune_taxa(taxa_sums(merged_FITS_all4_nounknowns) > 10, merged_FITS_all4_nounknowns) 
phyFITS_500 <- prune_samples(sample_sums(phyFITS) > 500, phyFITS)

#sediments only
phyFITS_500_sediments <- subset_samples(phyFITS_500, Substrate == "Sediment")
phyFITS_500_sediments <- prune_taxa(taxa_sums(phyFITS_500_sediments) > 10, phyFITS_500_sediments) 

#ordination
ord_FITS <- ordinate(phyFITS_500_sediments, "NMDS", "jaccard", k=2, trymax=100)
#ord_FITS <- ordinate(phy16S_500_oases_noK0601_A1, "NMDS", "jaccard", k=2, trymax=200,  previous.best=ord_16S)

#k=2 500 sediments only stress=0.19

#plot
plot_ord_FITS_12 <- plot_ordination(phyFITS_500_sediments, ord_FITS, type="samples", color="Oasis", shape="Timepoint", axes= 1:2,label = NULL)
plot_ord_FITS_12 <- plot_ord_FITS_12+ scale_color_manual(values=c(oasis_color)) + labs(title="FITS, stress=0.19, k=2")#+ facet_wrap(~Processing_Batch)
plot_ord_FITS_12


#faceted
plot_ord_FITS_facet <- plot_ordination(phyFITS_500_sediments, ord_FITS, type="samples", 
                                   color="Oasis", axes= 1:2,label = NULL) 
plot_ord_FITS_facet  <- plot_ord_FITS_facet  + geom_point(size=3)+theme_bw()  +  theme(aspect.ratio=1)+ scale_color_manual(labels= c("SP - Invaded","TPO - Pristine"), values=c(color_blind1)) + 
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
                                                           # labs(title="FITS - k=2, stress=0.19")
                                                         

plot_ord_FITS_facet

#2019_12 separating on x axis from other samples

#beta stats
dist_FITS <- distance(phyFITS_500_sediments, method="jaccard") #make distance matrix

metadata_FITS<- as(sample_data(phyFITS_500_sediments), "data.frame")

perm <- how(nperm=9999)
setBlocks(perm) <- with(metadata_FITS, Timepoint)
adonis_dist <- adonis2(dist_FITS~Oasis*Timepoint, data=metadata_FITS, permutations=perm) #run adonis (PERMANOVA)
adonis_dist

pairwiseAdonis <-pairwise.adonis2(dist_FITS~Oasis_Timepoint,data=metadata_FITS, strata="Timepoint", perm=9999)
pairwiseAdonis

#all oasis to oasis comparisons in same year are significantly different with strata=Timepoint, all other comparisons pr(>F)=1 


#betadispersion - if significant, adonis result could be due to differences in dispersion
beta <- betadisper(dist_FITS, metadata_FITS$Oasis_Timepoint)
permutest(beta)

TukeyHSD(beta) 
#2019_12 is significantly different in dispersion from other groups - so difference is likely due to that?


#alpha diversity
alpha.diversity <- estimate_richness(phyFITS_500_sediments, measures =c("Observed", "Chao1", "Shannon"))
data <- cbind(sample_data(phyFITS_500_sediments), alpha.diversity)
aov <- aov(Chao1~ Oasis*Timepoint, data=data)

summary(aov)

Tukey <- TukeyHSD(aov)
Tukey

p_alpha <- ggplot(data, x="Timepoint",  
                  y= "Chao1",
                  color="Oasis") 


p_alpha <- p_alpha + geom_boxplot(aes(x=Timepoint,
                                      y=Chao1, fill=Oasis), outlier.shape=NA) + scale_fill_manual(values=c(oasis_color)) #+theme(axis.text.x = element_text(angle = 45, hjust=0.75))
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  ggtitle("FITS Biodiversity") #+ scale_y_continuous(limits=c(0,25))
p_alpha







p_alpha <- plot_richness(phyFITS_500_sediments, x="Timepoint",
                         measures=c("Observed", "Chao1", "Shannon"),
                         color="Oasis", title = "FITS")# scale_color_manual(values = c("red", "blue", "green"))

p_alpha

p_alpha <- p_alpha + geom_boxplot(data=p_alpha$data, aes(x=Timepoint, y=value, fill=Oasis), alpha=0.1) + 
  scale_color_manual(values= oasis_blue)  #+theme_bw() 

p_alpha

#2020_06 highest biodiversity overall fpr both oases

#violin plot with mean points

p_alpha <- ggplot(data, aes(x=Timepoint,  
                            y= Chao1, fill=Oasis)) + geom_violin() + stat_summary(fun.data = "mean_cl_boot",
                                                                                  geom = "pointrange",
                                                                                  color = "white", size = 0.25, position=position_dodge(0.9))#+ geom_boxplot(width=0.1, color="white", alpha=0.2, position=position_dodge(0.9), outlier.shape=NA)


p_alpha <- p_alpha + scale_fill_manual(values=c(color_blind1), labels=c(
  "SP - Invaded", "TPO - Pristine")) + scale_x_discrete(labels=c("Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020")) 
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), axis.text.y=element_text(size=12)) + 
  labs(x="Sample Time Point",  y = "Chao1 Richness - CO1") +  #+ scale_y_continuous(limits=c(0,25))
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(size=12), legend.text = element_text(size=11)
  )
p_alpha
