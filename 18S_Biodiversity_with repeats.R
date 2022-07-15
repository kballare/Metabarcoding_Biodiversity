#18S Biodiversity
library(devtools)
library(plyr)
library(dplyr)

library(tibble)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(vegan)


library(microbiomeSeq) #R 4.0

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


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
color_blind2 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

#read in phyloseq object
merged_18S_all4 <- readRDS("merged_18S_all4.rds")
merged_18S_all4_withrepeats <- readRDS("merged_18S_all4_with_repeats1.rds")

#add new variables to metadata- add in excel, import and add individually as vbectors (make sure samples are in exact same order)
#merge phyloseq function should work too but can't figure it out
meta_18S_repeats <- sample_data(merged_18S_all4_withrepeats)
meta_18S_repeats <- as.matrix(meta_18S_repeats)
write.csv(meta_18S_repeats, "meta_18S_repeats.csv")

#make new variables in excel, in exact same order as original , reimport new variables
meta_18S_repeats_newvar <- as.data.frame(meta_18S_repeats_newvar)
#row.names(meta_18S_repeats_newvar) <- meta_18S_repeats_newvar$sum.taxonomy1
#meta_18S_repeats_newvar <- meta_18S_repeats_newvar[,-1]
meta_18S_repeats <- as.data.frame(meta_18S_repeats)

#create vectors
Sample_Name <- meta_18S_repeats_newvar$Sample_Name
Original_or_Repeat <- meta_18S_repeats_newvar$Original_or_Repeat
Processing_Batch <- meta_18S_repeats$Processing_Batch

#add to phyloseq object and check
sample_data(merged_18S_all4_withrepeats)$Sample_Name <- Sample_Name
sample_data(merged_18S_all4_withrepeats)$Original_or_Repeat <- Original_or_Repeat
sample_data(merged_18S_all4_withrepeats)$Processing_Batch <- Processing_Batch

#head(sample_data(merged_18S_all4_withrepeats))$Oasis_Timepoint
#sample_data(merged_18S_all4_withrepeats)$Oasis_Timepoint


#look for all unknown taxa to remove
write.csv(merged_18S_all4_withrepeats@tax_table, "18S_tax_repeats.csv")
tax_18S_repeats <- read.csv("18S_tax_repeats.csv")  #sort by Domain, Phylum, class etc and scroll to end to find unknown, unknown, unknown (sometimes in  the middle without sorting)
#1 taxon is all unknowns out of 2392 taxa

merged_18S_all4_withrepeats_nounknown <- subset_taxa(merged_18S_all4_withrepeats, Domain != "unknown") #remove fully unknown assigned taxa
merged_18S_all4_withrepeats_nounknown #check that number of taxa is what you expect after NA removal
#2391 taxa :)

merged_18S_all4_norepeats_nounknown <- subset_samples(merged_18S_all4_withrepeats_nounknown, Original_or_Repeat == "Original")
write.csv(merged_18S_all4_norepeats_nounknown@tax_table, "18_tax_norepeats.csv")

#make pruned phyloseq obj
phy18S <- prune_taxa(taxa_sums(merged_18S_all4_withrepeats_nounknown) > 10, merged_18S_all4_withrepeats_nounknown) #p
phy18S_500 <- prune_samples(sample_sums(phy18S) > 500, phy18S)

#oases only
phy18S_500_oases <- subset_samples(phy18S_500, Oasis != "Tank")
phy18S_500_oases <- prune_taxa(taxa_sums(phy18S_500_oases) > 10, phy18S_500_oases) 



#sediments only
phy18S_500_sediments <- subset_samples(phy18S_500, Substrate == "Sediment")
phy18S_500_sediments <- prune_taxa(taxa_sums(phy18S_500_sediments) > 10, phy18S_500_sediments) 
phy18S_500_sediments_norepeat <- subset_samples(phy18S_500_sediments, Original_or_Repeat != "Repeat")

#remove outliers
phy18S_500_sediments_noK0601_A1 <- subset_samples(phy18S_500_sediments, sample_names(phy18S_500_sediments)!="K0601_A1")
#phy18S_all4_500_rep_sediments_noK0601_A1 <- subset_samples(phy18S_all4_rep_500_sediments, sample_names(phy18S_all4_rep_500_sediments)!="K0601_A1")
#r <- rarefy_even_depth(phy18S_all4_500_sediments_noK0601_A1, sample.size = 5000, verbose = FALSE, replace = TRUE)
phy18S_500_oases_noK0601_A1 <- subset_samples(phy18S_500_oases, sample_names(phy18S_500_oases)!="K0601_A1")

phy18S_500_sediments_norepeat_noK0601_A1 <- subset_samples(phy18S_500_sediments_norepeat, sample_names(phy18S_500_sediments_norepeat)!="K0601_A1")


#make separate physeq objects for each timepoint
phy18S_2019_06 <- subset_samples(phy18S_500_sediments_noK0601_A1, Timepoint=="2019_06")
phy18S_2019_12 <- subset_samples(phy18S_500_sediments_noK0601_A1, Timepoint=="2019_12")
phy18S_2020_06 <- subset_samples(phy18S_500_sediments_noK0601_A1, Timepoint=="2020_06")
phy18S_2020_12 <- subset_samples(phy18S_500_sediments_noK0601_A1, Timepoint=="2020_12")

#make separate physeq obj for each oasis
phy18S_500_sediments_norepeat_SP <- subset_samples(phy18S_500_sediments_norepeat_noK0601_A1, Oasis=="SP")
phy18S_500_sediments_norepeat_TPO <- subset_samples(phy18S_500_sediments_norepeat_noK0601_A1, Oasis=="TPO")



#make separate physeq objs for each timepoint and oasis to test for betadispersion between repeats
phy18S_2019_06_SP <- subset_samples(phy18S_2019_06, Oasis=="SP")
phy18S_2019_06_TPO <- subset_samples(phy18S_2019_06, Oasis=="TPO")

#beta diversity
#use only sediment samples for beta diversity, 10 minumum reads/taxa, 500 minimum reads/sample

#conduct test - are communities different between sample date?
#ordination- jaccard- uses presence absence only
ord_18S <- ordinate(phy18S_500_sediments_norepeat_noK0601_A1, "NMDS", "jaccard", k=3, trymax=100)
ord_18S <- ordinate(phy18S_500_sediments_norepeat_noK0601_A1, "NMDS", "jaccard", k=2, trymax=500, previous.best=ord_18S)

#repeats
#oases only, 500 reads, outliers retained
#no convergence k=2 -4
#remove K0601_A1, oases only, no covergence k=2
#k=3, stress=0.156

#norepeats, sediments only
#k=3, stress = 0.167

#no faceting #oasis color
 
plot_ord_18S_12 <- plot_ordination(phy18S_500_sediments_norepeat_noK0601_A1, ord_18S, type="samples", 
                                   color="Oasis", shape = "Timepoint")
plot_ord_18S_12 <- plot_ord_18S_12 +theme_bw() + labs(title="k=3 \n stress=0.167") + geom_point(size=3) + theme(plot.title=element_text(hjust=0.05, margin=margin(t=20,b=-30))) + 
  scale_color_manual(values=c(color_blind1)) + scale_shape_manual(values=c(16, 17, 1, 2)) #+ scale_size_manual(values=c(10,10,10,10))
#+ facet_wrap(~Timepoint, labeller=labeller(Timepoint = c("2019_06" = "June 2019, pre-eradication",
plot_ord_18S_12                                                                                                     #"2019_12" = "Dec 2019",
                                                                                                    # "2020_06" = "June 2020",
      
#no faceting #oasis shape
plot_ord_18S_12 <- plot_ordination(phy18S_500_sediments_norepeat_noK0601_A1, ord_18S, type="samples", 
                                   color="Timepoint", shape = "Oasis")
plot_ord_18S_12 <- plot_ord_18S_12 +theme_bw() + labs(title="k=3 \n stress=0.167") + geom_point(size=3) + theme(plot.title=element_text(hjust=0.05, margin=margin(t=20,b=-30))) + 
  scale_color_manual(values=c(color_blind2)) + scale_shape_manual(values=c(16, 17)) #+ scale_size_manual(values=c(10,10,10,10))
#+ facet_wrap(~Timepoint, labeller=labeller(Timepoint = c("2019_06" = "June 2019, pre-eradication",
plot_ord_18S_12                                                                                                     # "2020_12" = "Dec 2020")))plot_ord_18S_12

#faceted
plot_ord_18S_12 <- plot_ordination(phy18S_500_sediments_norepeat_noK0601_A1, ord_18S, type="samples", 
                                   color="Oasis", axes= 1:2,label = NULL) 
plot_ord_18S_12 <- plot_ord_18S_12 + geom_point(size=3)+theme_bw() + labs(title="k=3, stress=0.167") + scale_color_manual(labels= c("SP - Invaded","TPO - Pristine"), values=c(color_blind1)) + 
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

plot_ord_18S_12


plot_ordination(phy16S_500_oases, ord_16S, type="samples", color="Oasis", shape= "Original_or_Repeat", axes= 1:2,label = NULL)
plot_ord_16S_12 <- plot_ord_16S_12  + scale_shape_manual(values=c(16,4)) + scale_color_manual(values=c(oasis_color, "red")) + facet_wrap(~Processing_Batch)
plot_ord_16S_12
 
plot_ord_18S_rep <- plot_ordination(phy18S_500_sediments_noK0601_A1, ord_18S, type="samples", color="Oasis", shape="Original_or_Repeat",axes= 1:2,label = NULL)
plot_ord_18S_rep <- plot_ord_18S_rep + facet_grid(~Timepoint) + scale_shape_manual(values=c(16,4)) + scale_color_manual(values=c(oasis_color))
plot_ord_18S_rep


plot_ord_18S_rep_13 <- plot_ordination(phy18S_500_sediments_noK0601_A1, ord_18S, type="samples", color="Oasis", shape="Original_or_Repeat",axes= 1:3,label = NULL)
plot_ord_18S_rep_13 <- plot_ord_18S_rep_13 + facet_grid(~Timepoint) + scale_shape_manual(values=c(16,4)) + scale_color_manual(values=c(oasis_color))
plot_ord_18S_rep_13

plot_ord_18S_rep_23 <- plot_ordination(phy18S_500_sediments_noK0601_A1, ord_18S, type="samples", color="Oasis", shape="Original_or_Repeat",axes= 2:3,label = NULL)
plot_ord_18S_rep_23 <- plot_ord_18S_rep_23 + facet_grid(~Timepoint) + scale_shape_manual(values=c(16,4)) + scale_color_manual(values=c(oasis_color))
plot_ord_18S_rep_23

#all repeats look pretty good, clustering mostly with originals on NMDS axis 1 and 2

MDS_18S <- ordinate(phy18S_500_sediments_noK0601_A1, "MDS", "jaccard")
plot_MDS_18S <- plot_ordination(phy18S_500_sediments_noK0601_A1, MDS_18S, type="samples", color="Oasis", shape=
                                "Original_or_Repeat", label = NULL)
plot_MDS_18S <- plot_MDS_18S + facet_grid(~Timepoint) + scale_shape_manual(values=c(8, 16))

plot_MDS_18S


#adonis with repeats
#test with TPO and SP separate? too many comparisions?
dist_18S <- distance(phy18S_500_oases_noK0601_A1, method="jaccard") #make distance matrix

metadata_18S<- as(sample_data(phy18S_500_oases_noK0601_A1), "data.frame")

adonis_dist <- adonis2(dist_18S~ Processing_Batch, data=metadata_18S, permutations=9999) #run adonis (PERMANOVA)

adonis_dist

pairwiseAdonis <-pairwise.adonis(dist_18S, metadata_18S$Processing_Batch, perm=9999)

pairwiseAdonis


#no repeats
dist_18S <- distance(phy18S_500_sediments_norepeat_noK0601_A1, method="jaccard") #make distance matrix

metadata_18S<- as(sample_data(phy18S_500_sediments_norepeat_noK0601_A1), "data.frame")

perm <- how(nperm=9999)
setBlocks(perm) <- with(metadata_18S, Timepoint)
adonis_dist <- adonis2(dist_18S~ Oasis_Timepoint, data=metadata_18S, permutations=perm) #run adonis (PERMANOVA)


adonis_dist


pairwiseAdonis <-pairwise.adonis2(dist_18S~Oasis_Timepoint, metadata_18S,strata="Timepoint", nperm=9999)

pairwiseAdonis

#betadispersion - if significant, adonis result could be due to differences in dispersion
beta <- betadisper(dist_18S, metadata_18S$Timepoint)
permutest(beta)
TukeyHSD(beta)

#beta stats tests by timepoint and pond
  dist_18S_2019_06_SP <- distance(phy18S_2019_06_SP, method="jaccard") #make distance matrix

metadata_18S_2019_06_SP<- as(sample_data(phy18S_2019_06_SP), "data.frame")

adonis_dist <- adonis2(dist_18S_2019_06_SP~ Oasis_Timepoint, data=metadata_18S_2019_06_SP, perm=9999) #run adonis (PERMANOVA)

adonis_dist

pairwiseAdonis <-pairwise.adonis(dist_18S_2019_06_SP, metadata_18S_2019_06_SP$Oasis_Timepoint, perm=9999)

pairwiseAdonis


        #pairs                    Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#1     TPO_2019_06 vs SP_2019_06  1 1.4553831 3.637226 0.07327617  0.0001     0.0006  **
#  2   TPO_2019_06 vs SP_2019_06_r  1 1.3682479 3.440277 0.09989110  0.0001     0.0006  **
#  3  TPO_2019_06 vs TPO_2019_06_r  1 0.6347645 1.543781 0.05408468  0.0166     0.0996    
#4    SP_2019_06 vs SP_2019_06_r  1 0.8370626 2.330071 0.08525668  0.0021     0.0126   .
#5   SP_2019_06 vs TPO_2019_06_r  1 0.7573163 2.051113 0.08898108  0.0041     0.0246   .
#6 SP_2019_06_r vs TPO_2019_06_r  1 0.7558314 2.704904 0.31073339  0.0341     0.2046   

#betadispersion - if significant, adonis result could be due to differences in dispersion/can't trust PERMANOVA results especially for unbalanced design
dist_18S_2019_06 <- distance(phy18S_2019_06, method="jaccard") #make distance matrix
metadata_18S_2019_06<- as(sample_data(phy18S_2019_06), "data.frame")
beta <- betadisper(dist_18S_2019_06, metadata_18S_2019_06$Oasis_Timepoint)
permutest(beta)

#ermutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
#Groups     3 0.27812 0.092708 13.163    999  0.001 ***
#  Residuals 52 0.36624 0.007043                                              

anova(beta)

#Response: Distances
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     3 0.27812 0.092708  13.163 1.639e-06 ***
 # Residuals 52 0.36624 0.007043   

TukeyHSD(beta)

#Fit: aov(formula = distances ~ group, data = df)

#$group
#diff         lwr         upr     p adj
#SP_2019_06_r-SP_2019_06    -0.13162878 -0.23473805 -0.02851951 0.0071558
#TPO_2019_06-SP_2019_06      0.03857417 -0.02623406  0.10338240 0.3989004
#TPO_2019_06_r-SP_2019_06   -0.26007662 -0.42490818 -0.09524506 0.0006157
#TPO_2019_06-SP_2019_06_r    0.17020295  0.06967182  0.27073408 0.0002251
#TPO_2019_06_r-SP_2019_06_r -0.12844784 -0.31031549  0.05341982 0.2514756
#TPO_2019_06_r-TPO_2019_06  -0.29865079 -0.46188201 -0.13541957 0.0000658

dist_18S_2019_12 <- distance(phy18S_2019_12, method="jaccard") #make distance matrix

metadata_18S_2019_12<- as(sample_data(phy18S_2019_12), "data.frame")

pairwiseAdonis <-pairwise.adonis(dist_18S_2019_12, metadata_18S_2019_12$Oasis_Timepoint, perm=9999)

pairwiseAdonis

#pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
#1     TPO_2019_12 vs SP_2019_12  1 0.7391943 1.711680 0.04204396  0.0008     0.0048   *
#  2   TPO_2019_12 vs SP_2019_12_r  1 0.6093765 1.444465 0.09352638  0.0038     0.0228   .
#3  TPO_2019_12 vs TPO_2019_12_r  1 0.4061759 0.953388 0.05976085  0.5556     1.0000    
#4    SP_2019_12 vs SP_2019_12_r  1 0.5340319 1.235695 0.03833314  0.0562     0.3372    
#5   SP_2019_12 vs TPO_2019_12_r  1 0.5909989 1.362373 0.04083561  0.0199     0.1194    
#6 SP_2019_12_r vs TPO_2019_12_r  1 0.4252241 1.010531 0.12615037  0.4365     1.0000    

#betadispersion - if significant, adonis result could be due to differences in dispersion/can't trust PERMANOVA results especially for unbalanced design
dist_18S_2019_12 <- distance(phy18S_2019_12, method="jaccard") #make distance matrix
metadata_18S_2019_12<- as(sample_data(phy18S_2019_12), "data.frame")
beta <- betadisper(dist_18S_2019_12, metadata_18S_2019_12$Oasis_Timepoint)
permutest(beta)

#Number of permutations: 999

#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
#Groups     3 0.042251 0.0140835 9.4147    999  0.001 ***
 # Residuals 46 0.068812 0.0014959   

anova(beta)

#Analysis of Variance Table

#Response: Distances
#Df   Sum Sq   Mean Sq F value   Pr(>F)    
#Groups     3 0.042251 0.0140835  9.4147 5.79e-05 ***
#  Residuals 46 0.068812 0.0014959     

TukeyHSD(beta)

#Fit: aov(formula = distances ~ group, data = df)

#$group
#diff          lwr         upr     p adj
#SP_2019_12_r-SP_2019_12    -0.09240895 -0.147395801 -0.03742211 0.0002796
#TPO_2019_12-SP_2019_12     -0.02395838 -0.059344516  0.01142776 0.2842750
#TPO_2019_12_r-SP_2019_12   -0.06359536 -0.113516712 -0.01367402 0.0074925
#TPO_2019_12-SP_2019_12_r    0.06845057  0.008929576  0.12797157 0.0183958
#TPO_2019_12_r-TPO_2019_12  -0.03963698 -0.094512634  0.01523866 0.2318239



dist_18S_2020_06 <- distance(phy18S_2020_06, method="jaccard") #make distance matrix

metadata_18S_2020_06<- as(sample_data(phy18S_2020_06), "data.frame")

pairwiseAdonis <-pairwise.adonis(dist_18S_2020_06, metadata_18S_2020_06$Oasis_Timepoint, perm=9999)

pairwiseAdonis

#pairs Df SumsOfSqs  F.Model         R2   p.value p.adjusted sig
#1     TPO_2020_06 vs SP_2020_06  1 1.0525037 2.630267 0.05640686 0.0001000     0.0006  **
 # 2  TPO_2020_06 vs TPO_2020_06_r  1 0.3722993 0.889535 0.05266782 0.7385000     1.0000    
#3   TPO_2020_06 vs SP_2020_06_r  1 0.5837648 1.412454 0.07671188 0.0145000     0.0870    
#4   SP_2020_06 vs TPO_2020_06_r  1 0.5720323 1.447478 0.04327615 0.0263000     0.1578    
#5    SP_2020_06 vs SP_2020_06_r  1 0.4688400 1.192361 0.03487213 0.1448000     0.8688    
#6 TPO_2020_06_r vs SP_2020_06_r  1 0.4386136 1.061718 0.17515130 0.2857143     1.0000    

dist_18S_2020_12 <- distance(phy18S_2020_12, method="jaccard") #make distance matrix

metadata_18S_2020_12<- as(sample_data(phy18S_2020_12), "data.frame")

pairwiseAdonis <-pairwise.adonis(dist_18S_2020_12, metadata_18S_2020_12$Oasis_Timepoint, perm=9999)

pairwiseAdonis

#pairs                            Df SumsOfSqs   F.Model         R2   p.value p.adjusted sig
#1     SP_2020_12 vs TPO_2020_12  1 0.9443867 2.3226857 0.06575620 0.0002000     0.0012   *
#  2    SP_2020_12 vs SP_2020_12_r  1 0.6612188 1.6202460 0.05124078 0.0070000     0.0420   .
#3   SP_2020_12 vs TPO_2020_12_r  1 0.5960038 1.4723444 0.05171139 0.0228000     0.1368    
#4   TPO_2020_12 vs SP_2020_12_r  1 0.6664880 1.5693764 0.12485714 0.0032000     0.0192   .
#5  TPO_2020_12 vs TPO_2020_12_r  1 0.3807539 0.9070537 0.10183544 0.6474000     1.0000    
#6 SP_2020_12_r vs TPO_2020_12_r  1 0.4278157 0.9796379 0.16382897 0.5714286     1.0000 


#### all repeats were not significant EXCEPT for 2019_06.  BetaDispersion tests was significant- not sure how to interpret
### no consistent evidence of batch effects?  add sample date stratum anyway to adonis with no repeats? 
### questions:  do we see differences in  oases over time? are oases more different in summer or winter? 
#what samples are they the most similar?  what  samples are they the most different?


#plot top 20 taxa
ps_tmp <- get_top_taxa(physeq_obj = phy18S_500_sediments_noK0601_A1, n = 21, relative = TRUE,
                       discard_other = TRUE, other_label = "Other")
ps_tmp <- name_taxa(ps_tmp, label = "unknown", species = T, other_label = "Other")

fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Species", other_label = "Other") 
?plot_bar

plot_bar(ps_tmp, x="Oasis_Timepoint", fill="Phylum")

ps_tmp_norm <- normalise_data(ps_tmp, norm.method="proportion")

p <- plot_taxa(ps_tmp_norm, grouping_column="Oasis_Timepoint", number.taxa=20, method="hellinger", filename="ps_tmp_norm_LCBD.csv")
p <- p + theme(legend.key.size = unit(0.10, 'cm'), axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p

library(adespatial)
plot_taxa <- function(physeq,grouping_column,method="hellinger",number.taxa=21,filename=NULL){
  #==extract components of the phyloseq object
  abund_table <- otu_table(physeq)
  meta_table <- data.frame(sample_data(physeq))
  #Enforce orientation of the phyloseq object
  if(taxa_are_rows(physeq) ){
    abund_table <- t(abund_table)
  }
  #===Calculate beta diversity and extract measure for local contribution to beta diversity
  beta_div<-beta.div(abund_table,method=method,sqrt.D=F,samp=T,nperm=999)
  df_LCBD<-data.frame(Sample=names(beta_div$LCBD),LCBD=beta_div$LCBD,p.LCBD=beta_div$p.LCBD)
  #=== add grouping information to the LCBD results
  df_LCBD<-data.frame(df_LCBD,Groups=meta_table[rownames(df_LCBD),grouping_column])
  if(!is.null(filename)){
    write.csv(df_LCBD,paste(filename,"_LCBD",".csv",sep=""))
  }
  select.top.taxa <- top.taxa(abund_table, number.taxa)
  new_x <- select.top.taxa$abund_table
  number.taxa <- select.top.taxa$number.taxa
  #arrange data for plotting in a format compatible to ggplot
  df<-NULL
  for (i in 1:dim(new_x)[2]){
    tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Groups=meta_table[,grouping_column])
    if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
  }
  df<-data.frame(df,df_LCBD[as.character(df$Sample),c("LCBD","p.LCBD","Groups")])
  #==plot the data
  colours <- microbiomeseq_cols()
  p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Groups, drop=TRUE,scale="free",space="free_x")
  p<- p+ guides(fill=guide_legend(ncol=1))+scale_fill_manual(values=colours[1:(number.taxa+1)])+theme_bw()+xlab("Samples")
  p<-p+ scale_y_continuous(expand = c(0.02,0))+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  p<-p+geom_point(aes(Sample,-0.02,size=LCBD))+theme(strip.background = element_rect(fill = "white"))
  return(p)
}
top.taxa <- function(abund_table, number.taxa){
  #==== sort the abundance table by total abundance of each taxa  in decreasing order
  abund_table<-abund_table[,order(colSums(abund_table),decreasing=TRUE)]
  #Extract list of top number.taxa Taxa
  taxa_list<-colnames(abund_table)[1:number.taxa]
  #remove "__Unknown__" and add it to others
  taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
  number.taxa<-length(taxa_list)
  #Generate a new table with everything added to Others
  new_x<-data.frame(abund_table[,colnames(abund_table) %in% taxa_list],Others=rowSums(abund_table[,!colnames(abund_table) %in% taxa_list]))
  out <- list("abund_table"=new_x, "number.taxa"=number.taxa)
  return(out)
}
microbiomeseq_cols <- function(){
  colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00",
               "#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000",
               "#FFFF00",grey.colors(1000));
  return(colours)
}

#####alpha diversity#####

alpha.diversity <- estimate_richness(phy18S_500_sediments_norepeat_noK0601_A1, measures =c("Observed", "Chao1", "Shannon"))
data <- cbind(sample_data(phy18S_500_sediments_norepeat_noK0601_A1), alpha.diversity)
aov <- aov(Chao1~ Oasis*Timepoint, data=data)

summary(aov)

Tukey <- TukeyHSD(aov)
Tukey

p_alpha <- plot_richness(phy18S_500_sediments_norepeat_noK0601_A1, x="Timepoint",
                         measures=c("Observed", "Chao1", "Shannon"),
                         color="Oasis", title = "18S")# scale_color_manual(values = c("red", "blue", "green"))


p_alpha <- p_alpha + geom_boxplot(data=p_alpha$data, aes(x=Timepoint, y=value, color=NULL), alpha=0.1) + 
  scale_color_manual(values= oasis_blue) #+theme_bw()

p_alpha
ggplot_build(p_alpha)$data


#violin plot
p_alpha <- ggplot(data, aes(x="Timepoint",  
                  y= "Chao1")) 


p_alpha <- p_alpha + geom_violin(aes(x=Timepoint,
                      y=Chao1, fill=Oasis)) + scale_fill_manual(values=c(color_blind))
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  ggtitle("18S Biodiversity") #+ scale_y_continuous(limits=c(0,25))
p_alpha


#violin plot with mean points

p_alpha <- ggplot(data, aes(x=Timepoint,  
                            y= Chao1, fill=Oasis)) + geom_violin() + stat_summary(fun.data = "mean_cl_boot",
                                                                                  geom = "pointrange",
                                                                                  color = "white", position=position_dodge(0.9))#+ geom_boxplot(width=0.1, color="white", alpha=0.2, position=position_dodge(0.9), outlier.shape=NA)
 

p_alpha <- p_alpha + scale_fill_manual(values=c(color_blind), labels=c(
  "SP - Invaded", "TPO - Pristine")) + scale_x_discrete(labels=c("Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020")) 
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), axis.text.y=element_text(size=12)) + 
  labs(x="Sample Time Point",  y = "Chao1 Richness - 18S") +  #+ scale_y_continuous(limits=c(0,25))
theme(
  plot.title = element_text(color="black", size=14, face="bold.italic"),
  axis.title.x = element_text(color="black", size=14, face="bold"),
  axis.title.y = element_text(color="black", size=14, face="bold"),
  legend.title = element_text(size=12), legend.text = element_text(size=11)
)



