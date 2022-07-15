######16S Biodiversity######
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

color_blind1 <- c("#D55E00" , "#0072B2") 

mixed <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
           "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
           "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
           "#8A7C64", "#599861")

oasis_color <- c("dodgerblue4", "goldenrod2")
oasis_blue <-c("yellowgreen", "dodgerblue4")

#####phyloseq obj creation and manipulation#####
#read in phyloseq object
merged_16S_all4_rep <- readRDS("merged_16S_all4_with_repeats1.rds")

#add new variables to metadata- add in excel, import and add individually as vbectors (make sure samples are in exact same order)
#merge phyloseq function should work too but can't figure it out
meta_16S_repeats <- sample_data(merged_16S_all4_rep)
meta_16S_repeats <- as.matrix(meta_16S_repeats)
write.csv(meta_16S_repeats, "meta_16S_repeats.csv")

#create new variables in a dataframe making sure in correct order (can create in excel and import csv) 
#create individual vectors of each variable
meta_16S_repeats <- read.csv("meta_16S_repeats.csv")
meta_16S_repeats <- as.data.frame(meta_16S_repeats)
Original_or_Repeat <- meta_16S_repeats$Original_or_Repeat
Processing_Batch <- meta_16S_repeats$Processing_Batch

#add to phyloseq object and check
sample_data(merged_16S_all4_rep)$Original_or_Repeat <- Original_or_Repeat
sample_data(merged_16S_all4_rep)$Processing_Batch <- Processing_Batch

head(sample_data(merged_16S_all4_rep))

#look for all unknown taxa to remove 
write.csv(merged_16S_all4_rep@tax_table, "16S_tax_repeats.csv")  #in Excel sort by Domain, Phylum, class etc and scroll to end to find unknown, unknown, unknown (sometimes in  the middle without sorting)
# 6055 taxa
#2 unknown taxa- one unknown, one blank domain

#filter out all unknowns
Domain_keep <- c("Eukaryota", "Archaea", "Bacteria")

merged_16S_all4_withrepeats_nounknown <- subset_taxa(merged_16S_all4_rep, Domain %in% Domain_keep) #remove fully unknown assigned taxa
merged_16S_all4_withrepeats_nounknown #check that number of taxa is what you expect after NA removal
#6053 taxa :)

#make pruned phyloseq obj
phy16S <- prune_taxa(taxa_sums(merged_16S_all4_withrepeats_nounknown) > 10, merged_16S_all4_withrepeats_nounknown) #p
phy16S_500 <- prune_samples(sample_sums(phy16S) > 500, phy16S)
phy16S_3000 <- prune_samples(sample_sums(phy16S) > 3000, phy16S)

#sediments only
phy16S_500_sediments <- subset_samples(phy16S_500, Substrate == "Sediment")
phy16S_500_sediments <- prune_taxa(taxa_sums(phy16S_500_sediments) > 10, phy16S_500_sediments) 

#remove outliers
phy16S_500_sediments_nooutliers <- subset_samples(phy16S_500_sediments, sample_names(phy16S_500_sediments)!="K0601_A1")
phy16S_500_sediments_nooutliers <- subset_samples(phy16S_500_sediments_nooutliers, sample_names(phy16S_500_sediments_nooutliers)!="K0738_B2")

phy16S_500_sed_norepeats_noK0601_A1 <- subset_samples(phy16S_500_sediments, sample_names(phy16S_500_sediments)!="K0601_A1")
phy16S_500_sed_norepeats_noK0601_A1 <- subset_samples(phy16S_500_sed_norepeats_noK0601_A1, Original_or_Repeat=="Original")
phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2 <- subset_samples(phy16S_500_sed_norepeats_noK0601_A1, sample_names(phy16S_500_sed_norepeats_noK0601_A1)!="K0738_B2")

#Oasis only
phy16S_500_oases <- subset_samples(phy16S_500, Oasis != "Tank")
phy16S_500_oases <- prune_taxa(taxa_sums(phy16S_500_oases) > 10, phy16S_500_oases) 

phy16S_500_oases_nooutliers <- subset_samples(phy16S_500_oases, sample_names(phy16S_500_oases)!="K0601_A1")
phy16S_500_oases_nooutliers <- subset_samples(phy16S_500_oases_nooutliers, sample_names(phy16S_500_oases_nooutliers)!="K0738_B2")
phy16S_500_oases_nooutliers <- subset_samples(phy16S_500_oases_nooutliers, sample_names(phy16S_500_oases_nooutliers)!="K0733_G5")


phy16S_500_oases_noK0601_A1 <- subset_samples(phy16S_500_oases, sample_names(phy16S_500_oases)!="K0601_A1")

phy16S_3000_oases <- subset_samples(phy16S_3000, Oasis != "Tank")
phy16S_3000_oases <- subset_samples(phy16S_3000_oases, Substrate != "Water")
phy16S_3000_oases <- prune_taxa(taxa_sums(phy16S_3000_oases) > 10, phy16S_3000_oases) 


#####ordinations#####
ord_16S <- ordinate(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2, "NMDS", "jaccard", k=2, trymax=100)
ord_16S <- ordinate(phy16S_500_oases_noK0601_A1, "NMDS", "jaccard", k=2, trymax=200,  previous.best=ord_16S)
#k=2 stress 0.18 (sediment only, no outliers 500)

#k=2 stress= 0.14 (all oasis samples included 500 ) **** going with this one for now

#k=2 stress 0.161 (all oases, outliers removed 500)

#k=2 stress 0.159 (all oases, 500 worst outlier removed (K0601_A1))

#k=2 stress 0.139 (all oasis samples, 3000)

#k=2 stress 0.123 (water samples removed, 3000)

#k=2 stress 0.168 no covergence (all samples, 3000)

#k=2 stress 0.17 (sediment only, no K0601_A1, 500, no repeats)

#k=2, stress = 0.185, (sediment only, no outliers, 500, no repeats)

plot_ord_16S_12 <- plot_ordination(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2, ord_16S, type="samples", color="Oasis", shape= "Original_or_Repeat", axes= 1:2,label = "sum.taxonomy")

plot_ord_16S_12 <- plot_ord_16S_12  + scale_shape_manual(values=c(16,4)) + scale_color_manual(values=c(oasis_color, "red")) + facet_wrap(~Processing_Batch)
plot_ord_16S_12

##saved plot with K0601 removed 500 minumum reads

#no faceting
plot_ord_16S_12 <- plot_ordination(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2, ord_16S, type="samples", 
                                   color="Oasis", shape= "Timepoint", axes= 1:2,label = "sum.taxonomy")
plot_ord_16S_12 <- plot_ord_16S_12 + theme_bw() + labs(title="16S, stress=0.185, k=2") +
  scale_color_manual(values=c(oasis_color)) #+ facet_wrap(~Timepoint, labeller=labeller(Timepoint = c("2019_06" = "June 2019, pre-eradication",

plot_ord_16S_12

#faceted NMDS
plot_ord_16S_facet <- plot_ordination(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2, ord_16S, type="samples", 
                                   color="Oasis", axes= 1:2,label = NULL) 
plot_ord_16S_facet <- plot_ord_16S_facet + geom_point(size=3)+theme_bw() +theme(aspect.ratio = 1) + scale_color_manual(labels= c("SP - Invaded","TPO - Pristine"), values=c(color_blind1)) + 
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
                                                           # + labs(title="16S - k=2, stress=0.185"))
                                                         

plot_ord_16S_facet

######plot top 20 taxa#####

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

ps_tmp <- get_top_taxa(physeq_obj = phy16S_500_oases_nooutliers, n = 21, relative = TRUE,
                       discard_other = TRUE, other_label = "Other")
ps_tmp <- name_taxa(ps_tmp, label = "unknown", species = T, other_label = "Other")

#fantaxtic_bar(ps_tmp, color_by = "Phylum", label_by = "Species", other_label = "Other") 


#plot_bar(ps_tmp, x="Oasis_Timepoint", fill="Phylum")

ps_tmp_norm <- normalise_data(ps_tmp, norm.method="proportion")

p <- plot_taxa(ps_tmp_norm, grouping_column="Oasis_Timepoint", number.taxa=20, method="hellinger", filename="ps_tmp_norm_LCBD.csv")
p <- p + theme(legend.key.size = unit(0.10, 'cm'), axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p

#####permanova and beta dispersion#####
  dist_16S <- distance(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2, method="jaccard") #make distance matrix

metadata_16S<- as(sample_data(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2), "data.frame")

perm <- how(nperm=9999)
setBlocks(perm) <- with(metadata_16S, Timepoint)
adonis_dist <- adonis2(dist_16S~Oasis*Timepoint, data=metadata_16S, permutations=perm) #run adonis (PERMANOVA)

pairwiseAdonis <-pairwise.adonis2(dist_16S~Oasis_Timepoint,data=metadata_16S, strata="Timepoint", perm=9999)

pairwiseAdonis
adonis_dist


#betadispersion - if significant, adonis result could be due to differences in dispersion
beta <- betadisper(dist_16S, metadata_16S$Oasis_Timepoint)
permutest(beta)

TukeyHSD(beta) 

#alpha diversity
alpha.diversity <- estimate_richness(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2, measures =c("Observed", "Chao1", "Shannon"))
data <- cbind(sample_data(phy16S_500_sed_norepeats_noK0601_A1_noK0738_B2), alpha.diversity)
aov <- aov(Observed~ Oasis*Timepoint, data=data)

summary(aov)

Tukey <- TukeyHSD(aov)
Tukey


#chao1 only
p_alpha <- ggplot(data, x="Timepoint",  y= "Chao1", color="Oasis") 
p_alpha <- p_alpha + geom_boxplot(aes(x=Timepoint,
             y=Chao1, fill=Oasis), outlier.shape=NA) + scale_fill_manual(values=c(oasis_color)) #+theme(axis.text.x = element_text(angle = 45, hjust=0.75))
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  ggtitle("16S Biodiversity") #+ scale_y_continuous(limits=c(0,25))
p_alpha

#facetted by alpha diversity measure
p_alpha <- plot_richness(phy16S_500_sed_norepeats_noK0601_A1, x="Timepoint",
                         measures=c("Observed", "Chao1", "Shannon"),
                         color="Oasis", title = "16S")# scale_color_manual(values = c("red", "blue", "green"))

p_alpha

p_alpha <- p_alpha + geom_boxplot(data=p_alpha$data, aes(x=Timepoint, y=value, color=NULL), alpha=0.1) + 
  scale_color_manual(values= oasis_blue) #+theme_bw()

p_alpha


#violin plot with mean points

p_alpha <- ggplot(data, aes(x=Timepoint,  
                            y= Chao1, fill=Oasis)) + geom_violin() + stat_summary(fun.data = "mean_cl_boot",
                                                                                  geom = "pointrange",
                                                                                  color = "white", position=position_dodge(0.9))#+ geom_boxplot(width=0.1, color="white", alpha=0.2, position=position_dodge(0.9), outlier.shape=NA)


p_alpha <- p_alpha + scale_fill_manual(values=c(color_blind1), labels=c(
  "SP - Invaded", "TPO - Pristine")) + scale_x_discrete(labels=c("Jun 2019", "Dec 2019", "Jun 2020", "Dec 2020")) 
p_alpha <- p_alpha + theme_bw()+ theme(axis.text.x = element_text(angle = 45, size=12, hjust=1, vjust=1), axis.text.y=element_text(size=12)) + 
  labs(x="Sample Time Point",  y = "Chao1 Richness - 16S") +  #+ scale_y_continuous(limits=c(0,25))
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(size=12), legend.text = element_text(size=11)
  )
p_alpha

#faceted NMDS
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
