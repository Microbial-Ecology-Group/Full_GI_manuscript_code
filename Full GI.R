# Library -----------------------------------------------------------------
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(vegan)
library(data.table)
library(metagenomeSeq)
library(stringr)
library(ggdendro)
library(cowplot)
library(pairwiseAdonis)
library(scales)
library(metagMisc)
library(GUniFrac)
library(viridis)
library(randomcoloR)
library(multcompView)
library(rcompanion)
library(readxl)
library(microbiome)



source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/changeSILVATaxaNames.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/g_unifrac.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/MergeLowAbund.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/uw_unifrac.R")
source("C:/Users/danie/OneDrive - West Texas A and M University/Bioninformatics/Source Scripts/w_unifrac.R")

# Example of a function that we can use (and re-use) to remove unwanted taxa
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

setwd("C:/Users/danie/OneDrive - West Texas A and M University/Research/Full GI/16S New Reads/")


# Color Palettes ----------------------------------------------------------
fgi_palette3 <- c("#D797F1",	"#89DC82",	"#AED86A",	"#4D53E0",	"#9AF798",	"#68DABE",	"#DEE77E",	"#8FAD48",	"#E0DDDA",	"#54D5DC",	"#E3762F",	"#96DDCA",	"#8B62AD",	"#4C8ED7",	"#84CBF1",	"#83E5B5",	"#CB7E79",	"#81AAAA",	"#45B48B",	"#DC647B",	"#EA83E6",	"#F4A6CC",	"#9848A1",	"#C5F46B",	"#D25793",	"#4BBFD1",	"#D4C3F6",	"#935758",	"#F7D649",	"#607E84",	"#58FAE2",	"#C6E089",	"#D5A3CE",	"#A29DD0",	"#479296",	"#EBEEF6",	"#6E6A26",	"#50F0F8",	"#3C7B5A",	"#B3B36C",	"#38F7EA",	"#EF965E",	"#857B99",	"#8FF74D",	"#F3AF31",	"#C77D95",	"#D4ACED",	"#D392D5",	"#BE9799",	"#A0DA41",	"#EFD1C6",	"#AB6CA6",	"#EDF270",	"#3DDEDB",	"#F15A6A",	"#D5F0B2",	"#728E63",	"#F39592",	"#4850A2",	"#EDF538",	"grey85")
oral_family_palette <- c("#D797F1",	"#AED86A",	"#F1F1C0",	"#9AF798",	"#68DABE",	"#DEE77E",	"#BEC54C",	"#A09FA7",	"#C7A644",	"#96DDCA",	"#4C8ED7",	"#C5E3EB",	"#84CBF1",	"#A851E4",	"#C72AF6",	"#81AAAA",	"#EA83E6",	"#F3BEED",	"#D25793",	"#A6C9A7",	"#BEE9D2",	"#D4C3F6",	"#935758",	"#8BADED",	"#F7D649",	"#C0F7A0",	"#607E84",	"#C6E089",	"#479296",	"#B679D5",	"#EBEEF6",	"#EFB19F",	"#3C7B5A",	"#38F7EA",	"#E7684F",	"#BBA7BD",	"#857B99",	"#8FF74D",	"#E9CC1D",	"#AF8FF6",	"#C77D95",	"#D392D5",	"#AB6CA6",	"#EDF270",	"#85C6D0",	"#F39592",	"#EDF538",	"#85E8F1",	"grey85")
rumen_family_palette <- c("#D797F1",	"#89DC82",	"#9AF798",	"#68DABE",	"#8FAD48",	"#E0DDDA",	"#427F23",	"#E3762F",	"#96DDCA",	"#8B62AD",	"#84CBF1",	"#83E5B5",	"#CB7E79",	"#81AAAA",	"#45B48B",	"#EA83E6",	"#F4A6CC",	"#9848A1",	"#C5F46B",	"#D943EB",	"#E1F2E2",	"#935758",	"#F7D649",	"#D5A3CE",	"#A29DD0",	"#EBEEF6",	"#CCCDA7",	"#38F7EA",	"#EF965E",	"#8FF74D",	"#F3AF31",	"#C77D95",	"#D4ACED",	"#D392D5",	"#BE9799",	"#A0DA41",	"#EFD1C6",	"#3DDEDB",	"#58F086",	"#D5F0B2",	"#CD5F5B",	"#728E63",	"#F39592",	"#4850A2",	"grey85")
abo_family_palette <- c("#D797F1",	"#EDCDA2",	"#89DC82",	"#AED86A",	"#9AF798",	"#68DABE",	"#88A5C4",	"#8FAD48",	"#E0DDDA",	"#E3762F",	"#96DDCA",	"#4C8ED7",	"#84CBF1",	"#83E5B5",	"#CB7E79",	"#81AAAA",	"#45B48B",	"#DC647B",	"#EA83E6",	"#F4A6CC",	"#9848A1",	"#C5F46B",	"#D943EB",	"#D4C3F6",	"#935758",	"#F7D649",	"#D5A3CE",	"#A29DD0",	"#EBEEF6",	"#38F7EA",	"#EF965E",	"#8FF74D",	"#F3AF31",	"#C77D95",	"#D4ACED",	"#D392D5",	"#BE9799",	"#A0DA41",	"#EFD1C6",	"#3DDEDB",	"#CD5F5B",	"#728E63",	"#F39592",	"#4850A2",	"grey85")
fecal_family_palette <- c("#D797F1",	"#89DC82",	"#4D53E0",	"#9AF798",	"#68DABE",	"#88A5C4",	"#DEE77E",	"#8FAD48",	"#71C6BD",	"#E0DDDA",	"#54D5DC",	"#E3762F",	"#96DDCA",	"#8B62AD",	"#4C8ED7",	"#84CBF1",	"#CB7E79",	"#81AAAA",	"#45B48B",	"#EA83E6",	"#F4A6CC",	"#9848A1",	"#C5F46B",	"#D943EB",	"#D4C3F6",	"#935758",	"#F7D649",	"#58FAE2",	"#D5A3CE",	"#EBEEF6",	"#6E6A26",	"#50F0F8",	"#B3B36C",	"#38F7EA",	"#EF965E",	"#8FF74D",	"#F3AF31",	"#C77D95",	"#D4ACED",	"#D392D5",	"#BE9799",	"#A0DA41",	"#EFD1C6",	"#3DDEDB",	"#DBF558",	"#D5F0B2",	"#728E63",	"#F39592")
la_family_palette <- c("#AED86A",	"#68DABE",	"#DEE77E"	,"#D25793"	,"#EBB85D"	,"#AB6CA6",	"#F3E6D2")

bodysite_palette <- c("red4", "royalblue1", "navy", "gold", "springgreen", "goldenrod4", 'goldenrod', "darkviolet", "darkorange1", "orangered2", "royalblue3")
phylum_palette <- c("#86E958",	"navy",	"#8777D4",	"#DFE36F",	"#D3B4E2",	"red2",		"#76B0D3",	"#A1AF7A",	"#DFAB57",	"#DC65B8",	"#D8E3D6",	"yellow3",	"Grey85")
firm_palette <- c("#c5ef99",	"maroon1",	"#88ef21",	"#ea877e",	"#c6aefc",	"#fff3af",	"#f9e189",	"#d5f293",	"#e29f7a",	"#e5799b",	"#04e2ea",	"#ccf79b",	"tomato",	"#df4eed",	"#16a37d",	"#f73722",	"#dc28fc",	"#3bdba3",	"#f7d296",	"#86ed42",	"#f4ae9f",	"#5b7705",	"deeppink",	"firebrick1",	"#59e5ce",	"#68a4cc",	"#fcb0de",	"#f9e957",	"#632393",	"#7ceac7",	"#f45a62",	"#a3f492",	"#32f26c",	"#db6980",	"#a4fcd1",	"#efe386",	"#ac64d6",	"#ed7f1e",	"#3570dd",	"#8a74ed",	"#b23140",	"#2261b5",	"#55fc71",	"#bef9a4",	"#697203",	"#db0697",	"#dd4840",	"#2a1570",	"#38589e",	"#33e854",	"#ea8c93",	"#92ef93",	"#f70978",	"#8752c4",	"red4",	"#8fe060",	"#776dd6",	"#9167e5",	"#3fba1a",	"#8860c9",	"#94efda",	"#698c10",	"#4587bc",	"#a8bf37",	"#14a5ff",	"#a0e4ff",	"red",	"#2c4caa",	"#f96240",	"#493ca8",	"#c1a801",	"#f2e44b",	"#4ff29e",	"indianred1",	"#6ed877",	"#e27c9e",	"#ff47e9",	"#d8634b",	"#8813bf",	"red",	"#015c89",	"#8884ed",	"#36db91",	"#bff7a0",	"#017f4d",	"#e28cd8",	"#2bc65a",	"orangered",	"#b73c3a",	"#e84ecb",	"#f9db40",	"#8bedf9",	"#fc0c90",	"#e8ed93",	"#c1f74c",	"#d8d52d",	"#63edbd",	"#c366e2",	"#55c8e8",	"#9eed50",	"#a1e2f4",	"#6772d3",	"#4de2cc",	"#d1894f",	"#ed9387",	"#16ff45",	"#4e63ba",	"#a0ea7e",	"#edbf5c",	"#8e1d04",	"#7fa2d1",	"#882fb7",	"#4adb15",	"#afed6d",	"#c4602b",	"#ffccd1",	"#ef92ce",	"deeppink4")
bac_palette <- c("#0ffcbd",	"mediumblue",	"#4080bc",	"#268ab5",	"#fcc9ae",	"#89edbd",	"#c7fca6",	"cyan",	"#a4f9c4",	"#c13a65",	"#67f7b6",	"#d3c504",	"#92e6f4",	"#e4f957",	"#d969e5",	"#a3fff2",	"#c6059d",	"#f6c9ff",	"cyan3",	"#175193",	"#568708",	"#3faaaa",	"#f4bd30",	"#db4a4c",	"#c8f77e",	"#db84ed",	"#57b2c6",	"#cd96f2",	"#f966f9",	"#a2f9dc",	"darkslategray1",	"#ead57e",	"darkblue",	"#f17cf9",	"darkcyan",	"#d77df2",	"cyan1",	"#efe57a",	"#f90429",	"#bc2044",	"lightskyblue1",	"#5968cc",	"#49026d",	"#ff7579",	"royalblue1",	"#ce14a9",	"#91ffc2",	"#52a30b",	"#ce4618",	"#2997a5",	"#fcb5cc",	"#3d0d9b",	"#8ba5f4",	"#c0ed65",	"#c91c6a",	"blue1",	"grey85")
actino_palette <- c("#fcc6a6",	"#daa5ef",	"darkgreen",	"#08e804",	"limegreen",	"#cfadff",	"chartreuse1",	"#90bbed",	"#ef94c9",	"#84e8cb",	"#f0b5fc",	"#f4abbb",	"#ddfcab",	"darkseagreen1",	"#20ea38",	"#277e96",	"#15bbc6",	"#986ae2",	"#82c440",	"#ffa67a",	"darkseagreen4",	"greenyellow",	"#109609",	"#a0c9f7",	"#88f7f3",	"#d66a17",	"#9e55cc",	"#f4ec7a",	"#383da0",	"#fcd01e",	"#39db2e",	"#e7fc94",	"#bc9bdd",	"#e59cd6",	"#f9b8c3",	"olivedrab4",	"#bffc5d",	"#b6bf07",	"#5e32b5",	"#ed55ea",	"#5eff5e",	"#62bcc4",	"#e27175",	"#e2bf0d",	"#0fafe0",	"#38c4b6",	"#52b717",	"#ebef0e",	"#d1fc99",	"#63e500",	"#beacf9",	"#aa52d3",	"#35eac9",	"#a093ff",	"#b82edb",	"mediumseagreen",	"#06447a",	"#c8ff91",	"#dcc6ff",	"#82ea6e",	"#fcc7f0",	"#8df449",	"#f4abd7",	"#d6c50e",	"#2dacb7",	"#ce1069",	"#f2b3d4",	"#4b51aa",	"#f0c4ff",	"#d6972a",	"#4fef79",	"#81d358",	"#ed97d0")

no_la_palette <- c("red4", "royalblue1", "navy", "gold", "springgreen", "goldenrod4", 'goldenrod', "darkorange1", "orangered2", "royalblue3")
harvest_palette <- c("grey0", "brown")
paired_palette <- c("grey0", "brown","red4", "royalblue1", "navy", "gold", "springgreen", "goldenrod4", 'goldenrod', "darkviolet", "darkorange1", "orangered2", "royalblue3")

foregut_palete <- c("red4", "orangered2" )
SI_palette <- c("gold", "goldenrod4", "goldenrod")
LI_palette <- c("royalblue1", "navy", "royalblue3")
other_palette <- c("springgreen", "darkviolet","darkorange")

no_la_paired_pallette <- c("grey0", "brown","red4", "royalblue1", "navy", "gold", "springgreen", "goldenrod4", 'goldenrod', "darkorange1", "orangered2", "royalblue3")



# Import Data -------------------------------------------------------------


changeSILVAtaxa <- function(x) {
  # remove the D__ etc...
  tax.clean <- data.frame(row.names = row.names(x),
                          Kingdom = str_replace(x[,1], "k__",""),
                          Phylum = str_replace(x[,2], "p__",""),
                          Class = str_replace(x[,3], "c__",""),
                          Order = str_replace(x[,4], "o__",""),
                          Family = str_replace(x[,5], "f__",""),
                          Genus = str_replace(x[,6], "g__",""),
                          Species = str_replace(x[,7],"s_",""),
                          stringsAsFactors = FALSE)}

#read in microbiome data
fgidata <- import_biom('table-with-taxonomy.biom', 'tree.nwk', 'dna-sequences.fasta')

#Write Sample names to create metadata
#write.csv(sample_names(fgidata), "samplenames1.csv")

#read in metadata
metadata <- import_qiime_sample_data('samplenames.txt')

#merge metadata
data <- merge_phyloseq(fgidata, metadata)
data #193 samples 77896 taxa


# Data Exploration --------------------------------------------------------

#names pre processing
#rank then paste better names
rank_names(data)
colnames(tax_table(data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
rank_names(data)

#change the silva style names
tax.data <- data.frame(tax_table(data)) 
#changed this from changesilvanames, because it was reverting back to domain instead of kingdom
tax.data.names <- changeSILVAtaxa(tax.data)

#change NA's to better names
#convert to charcters
for(i in 1:7){ tax.data.names[,i] <- as.character(tax.data.names[,i])}
#replace with empty string
tax.data.names[is.na(tax.data.names)] <- ""

head(tax.data)

# now filling in the empty slots with the highest assigned taxonomy
for (i in 1:nrow(tax.data.names)){
  if (tax.data.names[i,2] == ""){
    domain <- paste("unclassified ", tax.data.names[i,1], sep = "")
    tax.data.names[i, 2:7] <- domain
  } else if (tax.data.names[i,3] == ""){
    kingdom <- paste("unclassified ", tax.data.names[i,2], sep = "")
    tax.data.names[i, 3:7] <- kingdom
  } else if (tax.data.names[i,4] == ""){
    phylum <- paste("unclassified ", tax.data.names[i,3], sep = "")
    tax.data.names[i, 4:7] <- phylum
  } else if (tax.data.names[i,5] == ""){
    class <- paste("unclassified ", tax.data.names[i,4], sep = "")
    tax.data.names[i, 5:7] <- class
  } else if (tax.data.names[i,6] == ""){
    order <- paste("unclassified ", tax.data.names[i,5], sep = "")
    tax.data.names[i, 6:7] <- order
  } else if (tax.data.names[i,7] == ""){
    tax.data.names$Genus[i] <- paste("unclassified ",tax.data.names$Family[i], sep = "")
  }
}

#check headers re-insert matrix and check tail
head(tax.data.names)
tax_table(data) <- as.matrix(tax.data.names)
tail(tax_table(data))

#prune for taka within Eukaryota
data <- subset_taxa(data, Domain!="Eukaryota")
data


#split samples and controls
samples <- subset_samples(data, Sample_Type == "Sample")
samples <- prune_taxa(taxa_sums(samples) > 0, samples)
samples #77777 taxa 176 samples

controls <- subset_samples(data, Sample_Type =="Control")
controls <- prune_taxa(taxa_sums(controls) > 0, controls)
controls #1175 taxa 17 sample
sort(sample_sums(controls)) #NTC3  NTC5 NTC13 NTC14  NTC4  NTC1 NTC15 NTC10  NTC2 NTC17 NTC12  NTC7 NTC16 NTC18  NTC8   EB1 NTC11 
                            #21    22    23    23    30    35    38    80    91    92    94   116   197   374  1272  2250 45221 
controls <- prune_samples(sample_sums(controls)> 100, controls)

#investigate the ntc with reads
controls_genus <-  tax_glom(controls, taxrank = "Genus", NArm = F) %>% 
  psmelt()
length(unique(controls_genus$Genus)) #220 genera
ggplot(controls_genus, aes(x=Sample_names, y=Abundance, fill= Genus))+
  geom_bar(stat = "identity")


#look at number of reads per sample and distributions
sample_sum_df <- data.frame(sum= sample_sums(samples)) 
ggplot(sample_sum_df, aes(x= sum))+
  geom_histogram(color= "black", fill = "indianred", binwidth = 5000)+
  ggtitle("Distribution of sample sequencing depth")+ 
  xlab("Read Counts")+
  theme(axis.title.y = element_blank())

mean(sample_sums(samples)) #1,176,341
median(sample_sums(samples)) #1,223,189
sort(sample_sums(samples)) #several 1's in the duodenum and lots below 10k

#cut off below 300,000 
samples <- prune_samples(sample_sums(samples) > 300000, samples)
samples #cuts us down to 163 samples
samples <- prune_taxa(taxa_sums(samples) > 0, samples)
samples #74565 taxa and 163 samples remain
min(sample_sums(samples)) # 312634
max(sample_sums(samples))#2,363,766
mean(sample_sums(samples))#1,267,216
sort(sample_sums(samples))
sum(sample_sums(samples))

# Comparing Between Samples -----------------------------------------------
sample_sum_df_new <- data.frame(ASV_count = sample_sums(samples))
metadata.df <- as(sample_data(samples), "data.frame") 
seqDepth_metadata <- cbind(metadata.df, sample_sum_df_new)

#comparing treatment groups for mootral
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$TRT) # NS

#comparing across liver abscess (not for severity)
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Liver_Abscess) #NS

#comparing for block
kruskal.test(seqDepth_metadata$ASV_count, seqDepth_metadata$Sample_Time)#NS


# Comparing sequencing depth ----------------------------------------------

#orals
oral_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Oral"),]
#by trt
kruskal.test(oral_seq_depth$ASV_count, oral_seq_depth$TRT) #NS
#by liver
kruskal.test(oral_seq_depth$ASV_count, oral_seq_depth$Liver_Abscess)#NS

#rumen
rumen_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Rumen"),]
#by trt
kruskal.test(rumen_seq_depth$ASV_count, rumen_seq_depth$TRT) #NS
#by liver
kruskal.test(rumen_seq_depth$ASV_count, rumen_seq_depth$Liver_Abscess)#NS

#abomasum
abo_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Abomasum"),]
#by trt
kruskal.test(abo_seq_depth$ASV_count, abo_seq_depth$Liver_Abscess)#NS
#by liver
kruskal.test(abo_seq_depth$ASV_count, abo_seq_depth$Liver_Abscess)#NS

#duo
duo_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Duodenum"),]
#by trt
kruskal.test(duo_seq_depth$ASV_count, duo_seq_depth$TRT) #NS
#by liver
kruskal.test(duo_seq_depth$ASV_count, duo_seq_depth$Liver_Abscess)#NS

#jej
jej_seq_dep <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Jejunum"),]
#by trt
kruskal.test(jej_seq_dep$ASV_count, jej_seq_dep$TRT)#NS
#by liver
kruskal.test(jej_seq_dep$ASV_count, jej_seq_dep$Liver_Abscess) #NS

#ile
ile_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Ileum"),]
#by trt
kruskal.test(ile_seq_depth$ASV_count, ile_seq_depth$TRT)#NS
#by liver
kruskal.test(ile_seq_depth$ASV_count, ile_seq_depth$Liver_Abscess)#NS

#cec
cec_seq_depth <-  seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Cecum"),]
#by trt
kruskal.test(cec_seq_depth$ASV_count, cec_seq_depth$TRT)#NS
#by liver
kruskal.test(cec_seq_depth$ASV_count, cec_seq_depth$Liver_Abscess)#NS

#sprC
sprc_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Spiral_Colon"),]
#by trt
kruskal.test(sprc_seq_depth$ASV_count, sprc_seq_depth$TRT)#NS
#by liver
kruskal.test(sprc_seq_depth$ASV_count, sprc_seq_depth$Liver_Abscess)#NS

#disC
disc_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Distal_Colon"),]
#by trt
kruskal.test(disc_seq_depth$ASV_count, disc_seq_depth$TRT)#NS
#by liver
kruskal.test(disc_seq_depth$ASV_count, disc_seq_depth$Liver_Abscess)#NS

#fecal
fecal_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Fecal"),]
#by trt
kruskal.test(fecal_seq_depth$ASV_count, fecal_seq_depth$TRT)# p= 0.059
#by liver
kruskal.test(fecal_seq_depth$ASV_count, fecal_seq_depth$Liver_Abscess)#NS

#livers
la_seq_depth <- seqDepth_metadata[which(seqDepth_metadata$Body_Site=="Liver_Abscess"),]
#by trt
kruskal.test(la_seq_depth$ASV_count, la_seq_depth$TRT)#NS

# Rarefaction Curve -------------------------------------------------------
#subset ASV table
otu_matrix <- otu_table(samples)
otu_matrix <- as.matrix(otu_matrix)


#Calculate rarefaction curve
set.seed(42)

calculate_rarefaction_curves <- function(otu_matrix, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(otu_matrix, measures, depth) {
    if(max(sample_sums(otu_matrix)) < depth) return()
    otu_matrix <- prune_samples(sample_sums(otu_matrix) >= depth, otu_matrix)
    
    rarified_otu_matrix <- rarefy_even_depth(otu_matrix, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_otu_matrix, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, otu_matrix = otu_matrix, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

#calculate observed and Shannon at different depths (caps at 1 mil)
rarefaction_curve_data <- calculate_rarefaction_curves(otu_matrix, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(samples)), by.x = 'Sample', by.y = 'row.names')


View(rarefaction_curve_data_summary_verbose)



# Create the paired plot
ggplot(
  data = rarefaction_curve_data_summary_verbose,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour = Body_Site,
    group = Sample
  )
) +
  geom_line() +
  geom_pointrange() +
  facet_wrap(
    facets = ~ Measure,
    scales = 'free_y'
  ) +
  scale_colour_manual(values = bodysite_palette) +  # Apply the custom color palette
  theme_bw() +
  labs(y = "", x = "", title = "") 
        theme(
    # Customize other plot elements as needed
    plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),
    panel.border = element_rect(colour = "black", size = 1.7),
    axis.ticks = element_line(size = 1, colour = "black"),
    plot.title = element_text(size = 30),
    axis.title.y = element_text(size = 30, vjust = 2.5),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

#subset observed data
sub_rare_data <- subset(rarefaction_curve_data_summary_verbose, Measure == 'Observed')

#only plot observed
rare_plot <- ggplot(data = sub_rare_data, 
                    mapping = aes(
                      x = Depth,
                      y = Alpha_diversity_mean,
                      ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                      ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                      colour = Body_Site, 
                      group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = bodysite_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)))+
  labs(y = "", x = "", title = "")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

rare_plot

# Rarefaction by region ---------------------------------------------------


#subset foregut
foregut_rare <- subset(sub_rare_data, Body_Site %in% c("Abomasum", "Rumen"))

#only plot foregut
foregut_rare_plot <- ggplot(data = foregut_rare, 
                    mapping = aes(
                      x = Depth,
                      y = Alpha_diversity_mean,
                      ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                      ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                      colour = Body_Site, 
                      group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = foregut_palete)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)))+
  ylim(0,6000)+
  labs(y = "", x = "", title = "")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

foregut_rare_plot

#subset SI
SI_rare <- subset(sub_rare_data, Body_Site %in% c("Jejunum", "Duodenum", "Ileum"))

#only plot foregut
SI_rare_plot <- ggplot(data = SI_rare, 
                            mapping = aes(
                              x = Depth,
                              y = Alpha_diversity_mean,
                              ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                              ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                              colour = Body_Site, 
                              group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = SI_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)))+
  ylim(0,6000)+
  labs(y = "", x = "", title = "")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

SI_rare_plot

#subset LI
LI_rare <- subset(sub_rare_data, Body_Site %in% c("Cecum", "Distal_Colon", "Spiral_Colon"))

#only plot foregut
LI_rare_plot <- ggplot(data = LI_rare, 
                            mapping = aes(
                              x = Depth,
                              y = Alpha_diversity_mean,
                              ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                              ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                              colour = Body_Site, 
                              group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = LI_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)))+
  ylim(0,6000)+
  labs(y = "", x = "", title = "")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

LI_rare_plot

#subset other
other_rare <- subset(sub_rare_data, Body_Site %in% c("Oral", "Fecal", "Liver_Abscess"))

#only plot other
other_rare_plot <- ggplot(data = other_rare, 
                       mapping = aes(
                         x = Depth,
                         y = Alpha_diversity_mean,
                         ymin = Alpha_diversity_mean - Alpha_diversity_sd,
                         ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                         colour = Body_Site, 
                         group = Sample)) +
  geom_line(size = 1) +
  scale_colour_manual(values = other_palette)+
  theme_bw() +
  scale_x_continuous(breaks = round(seq(0,max(sub_rare_data$Depth), by = 25000)))+
  ylim(0,6000)+
  labs(y = "", x = "", title = "")+
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75,0.75,0.75), "cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

other_rare_plot

# Checking classification efficiency --------------------------------------

# First we have to extract the tax_table and convert it to a "data.frame"
taxa.df <- as.data.frame(tax_table(samples))
taxa.df
#KINGDOM
# Now let's search for which taxa start with "unclassified"
unclassified_kingdom.df <- taxa.df %>% filter(grepl('Unassigned', Kingdom))
# Pull out the unique 
unclassified_kingdom <- row.names(unclassified_kingdom.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_kingdom.ps = pop_taxa(samples, unclassified_kingdom)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_kingdom.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_kingdom.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_kingdom.ps)) / sum(sample_sums(samples)) * 100

#Phylum
# Now let's search for which taxa start with "unclassified"
unclassified_phy.df <- taxa.df %>% filter(grepl('unclassified', Phylum))
# Pull out the unique 
unclassified_phy <- row.names(unclassified_phy.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_phy.ps = pop_taxa(samples, unclassified_phy)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_phy.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_phy.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_phy.ps)) / sum(sample_sums(samples)) * 100

#CLASS
# Now let's search for which taxa start with "unclassified"
unclassified_class.df <- taxa.df %>% filter(grepl('unclassified', Class))
# Pull out the unique 
unclassified_class <- row.names(unclassified_class.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_class.ps = pop_taxa(samples, unclassified_class)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_class.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_class.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_class.ps)) / sum(sample_sums(samples)) * 100

#ORDER
# Now let's search for which taxa start with "unclassified"
unclassified_order.df <- taxa.df %>% filter(grepl('unclassified', Order))
# Pull out the unique 
unclassified_order <- row.names(unclassified_order.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_order.ps = pop_taxa(samples, unclassified_order)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_order.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_order.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_order.ps)) / sum(sample_sums(samples)) * 100

#FAMILY
# Now let's search for which taxa start with "unclassified"
unclassified_family.df <- taxa.df %>% filter(grepl('unclassified', Family))
# Pull out the unique 
unclassified_family <- row.names(unclassified_family.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_family.ps = pop_taxa(samples, unclassified_family)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_family.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_family.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_family.ps)) / sum(sample_sums(samples)) * 100


#GENUS
# Now let's search for which taxa start with "unclassified"
unclassified_genus.df <- taxa.df %>% filter(grepl('unclassified', Genus))
# Pull out the unique 
unclassified_genus <- row.names(unclassified_genus.df)

# We can run the "pop_taxa" function here using the genus phyloseq object and the list of features we want to remove
trimmed_microbiome_genus.ps = pop_taxa(samples, unclassified_genus)

# Sum the counts for each sample
sample_sums(trimmed_microbiome_genus.ps)

# To sum across all samples, we can place that command inside the "sum()" function
sum(sample_sums(trimmed_microbiome_genus.ps))

# Now, we calculate the proportion of reads mapped to the species level, out of all microbiome mapped reads
sum(sample_sums(trimmed_microbiome_genus.ps)) / sum(sample_sums(samples)) * 100




# Alpha Diversity ---------------------------------------------------------

alphadiv1 <- estimate_richness(samples, measures = c("Observed", "Shannon", "Simpson", "InvSimpson"))
alphadiv2 <- estimate_pd(samples)

alpha_div <- cbind(alphadiv1, alphadiv2)
alpha_div

alpha_div <- alpha_div[,c(1:5)]
alpha_div
alpha_div_meta <- cbind(metadata.df,alpha_div)
alpha_div_meta

#pairwise for observed by body site
obs_anova <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$Body_Site, p.adjust.method = "BH")
#superscript p-values
obs_pvalues <- obs_anova$p.value %>% 
  fullPTable()
multcompLetters(obs_pvalues)

#pairwise for observed harvest group
pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$Sample_Time, p.adjust.method = "BH") #p < 0.001

##being influenced by LA look with no LA
#subset
alpha_no_la <- subset(alpha_div_meta, Body_Site != "Liver_Abscess")
alpha_no_la

#pairwise for observed
pairwise.wilcox.test(alpha_no_la$Observed, alpha_no_la$Sample_Time, p.adjust.method = "BH") #p = 0.001

#shannon
pairwise.wilcox.test(alpha_no_la$Shannon, alpha_no_la$Sample_Time, p.adjust.method = "BH")#p=0.009

#pairwise for shannon by body Site
shan_anova <- pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$Body_Site, p.adjust.method = "BH")
shan_anova

#superscript p values
shan_pvalues <- shan_anova$p.value %>% 
  fullPTable()
multcompLetters(shan_pvalues)

#pairwise for shannon by harvest group
pairwise.wilcox.test(alpha_div_meta$Shannon, alpha_div_meta$Sample_Time, p.adjust.method = "BH") #p =0.002

#plots
#body site
microbiome_richness_boxplot <- ggplot(alpha_div_meta, aes(x= Body_Site, y= Observed, fill = Body_Site)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = bodysite_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

microbiome_richness_boxplot

#set factor
alpha_div_meta$Sample_Time <- as.factor(alpha_div_meta$Sample_Time)
#harvest group
harvest_richness_boxplot <- ggplot(alpha_div_meta, aes(x= Sample_Time, y= Observed, fill = Sample_Time)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = harvest_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
harvest_richness_boxplot

#by bodysite
microbiome_shannon_boxplot <- ggplot(alpha_div_meta, aes(x= Body_Site, y= Shannon, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  scale_y_continuous(limits = c(2,6.5))+
  scale_fill_manual(values= bodysite_palette)+
  geom_boxplot() +
  geom_point()+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

microbiome_shannon_boxplot

#by harvest group
harvest_shannon_boxplot <- ggplot(alpha_div_meta, aes(x= Sample_Time, y= Shannon, fill = Sample_Time)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  #scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  #geom_text(label = alpha_div_meta$Sample_names)+
  scale_fill_manual(values = harvest_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
harvest_shannon_boxplot

plot_grid(microbiome_richness_boxplot,microbiome_shannon_boxplot, align = "h", axis = "rl", ncol = 2)

microbiome_PD_boxplot <- ggplot(alpha_div_meta, aes(x= Body_Site, y= PD, fill = Body_Site)) +
  theme_bw() + labs(y= "PD", x="", title = "Faiths Distance") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+ 
  scale_fill_manual(values = bodysite_palette)+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank())
microbiome_PD_boxplot


# Beta Diversity ----------------------------------------------------------

####CSS Transform###
any(taxa_sums(samples)==0)
samples.css <- phyloseq_transform_css(samples, log = F)
samples.css.df <- as(sample_data(samples.css), "data.frame")

#split by body part
#oral samples
oral.css <- subset_samples(samples.css, Body_Site == "Oral",)
oral.css <- prune_taxa(taxa_sums(oral.css)>0, oral.css)
oral.css#18 samples 18212 taxa
oral.css.df <- as(sample_data(oral.css), "data.frame")

#rumen samples
rumen.css <- subset_samples(samples.css, Body_Site =="Rumen")
rumen.css <- prune_taxa(taxa_sums(rumen.css)>0, rumen.css)
rumen.css #18 samples 14573 taxa
rumen.css.df <- as(sample_data(rumen.css), "data.frame")


#abo
abo.css <- subset_samples(samples.css, Body_Site=="Abomasum")
abo.css <- prune_taxa(taxa_sums(abo.css)>0, abo.css)
abo.css #18 samples 15120 taxa
abo.css.df <- as(sample_data(abo.css), "data.frame")

#duo
duo.css <- subset_samples(samples.css, Body_Site=="Duodenum")
duo.css <- prune_taxa(taxa_sums(duo.css)>0, duo.css)
duo.css #11 samples 15332 taxa
duo.css.df <- as(sample_data(duo.css), "data.frame")

#jej
jej.css <- subset_samples(samples.css, Body_Site=="Jejunum")
jej.css <- prune_taxa(taxa_sums(jej.css)>0, jej.css)
jej.css #16 samples 17815 taxa
jej.css.df <- as(sample_data(jej.css), "data.frame")

#ile
ile.css <- subset_samples(samples.css, Body_Site=="Ileum")
ile.css <- prune_taxa(taxa_sums(ile.css)>0, ile.css)
ile.css #13 samples 17491 taxa
ile.css.df <- as(sample_data(ile.css), "data.frame")

#cec
cec.css <- subset_samples(samples.css, Body_Site=="Cecum")
cec.css <- prune_taxa(taxa_sums(cec.css)>0, cec.css)
cec.css #17 samples and 17083 taxa
cec.css.df <- as(sample_data(cec.css), "data.frame")

#sprC
sprC.css <- subset_samples(samples.css, Body_Site=="Spiral_Colon")
sprC.css <-  prune_taxa(taxa_sums(sprC.css)>0, sprC.css)
sprC.css #15 samples 19510 taxa
sprC.css.df <- as(sample_data(sprC.css), "data.frame")

#disC 
disC.css <- subset_samples(samples.css, Body_Site=="Distal_Colon")
disC.css <- prune_taxa(taxa_sums(disC.css)>0, disC.css)
disC.css #18 samples 22810 taxa
disC.css.df <- as(sample_data(disC.css), "data.frame")


#fecal
fecal.css <- subset_samples(samples.css, Body_Site=="Fecal")
fecal.css <- prune_taxa(taxa_sums(fecal.css)>0, fecal.css)
fecal.css #15 samples 17860 taxa
fecal.css.df <- as(sample_data(fecal.css), "data.frame")

#LA
la.css <- subset_samples(samples.css, Body_Site=="Liver_Abscess")
la.css <- prune_taxa(taxa_sums(la.css)>0, la.css)
la.css #4 samples 992 taxa
la.css.df <- as(sample_data(la.css), "data.frame")


# Unifrac and Ordination by body site -------------------------------------

samples.dist <- gunifrac(samples.css)
samples.ord <- ordinate(samples.css, method = "NMDS", distance = samples.dist)



# Plots -------------------------------------------------------------------
#Liver Score
plot_ordination(samples.css, samples.ord, type = "samples", color = "Liver_Score")+
  theme_bw()+
  labs(title = "Liver Score", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  stat_ellipse(geom = "polygon", aes(fill= Liver_Score), alpha = 0.3, lty= 2, size= 1)+
  scale_colour_manual(values = distinctColorPalette(20))+
  scale_fill_manual(values = distinctColorPalette(20))+
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))

adonis2(samples.dist ~ Liver_Score, samples.css.df)# p = 0.001
samples.disper <-  betadisper(samples.dist, samples.css.df$Liver_Score)
plot(samples.disper)
samples.permdisp <- permutest(samples.disper, permutations = 9999)
samples.permdisp #p = 0.2231

#GIT location
plot_ordination(samples.css, samples.ord, type = "samples", color = "Body_Site")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= Body_Site), alpha = 0.3, lty= 2, size= 1)+
  scale_colour_manual(values = bodysite_palette)+
  scale_fill_manual(values = bodysite_palette)+
  #(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  theme(legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24))

adonis2(samples.dist ~ Body_Site, samples.css.df)#p =0.001
samples.disper <-  betadisper(samples.dist, samples.css.df$Body_Site)
plot(samples.disper)
samples.permdisp <- permutest(samples.disper, permutations = 9999)
samples.permdisp #p = 0.00001
#pairwise
ord_anova <- pairwise.adonis2(samples.dist ~ Body_Site, samples.css.df, nperm = 9999, p.adjust.m= "BH")


pairdisp <- betadisper(samples.dist, samples.css.df$Body_Site)
plot(pairdisp)
permutest(pairdisp, permutations = 9999, pairwise = T)


#by sample group
samples.css@sam_data$Sample_Time <- as.factor(samples.css@sam_data$Sample_Time)

plot_ordination(samples.css, samples.ord, type = "samples", color = "Sample_Time")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= Sample_Time), alpha = 0.3, lty= 2, size= 1)+
  scale_colour_manual(values = harvest_palette)+
  scale_fill_manual(values = harvest_palette)+
  #(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

adonis2(samples.dist ~ Sample_Time, samples.css.df)#p = 0.001
samples.disper <-  betadisper(samples.dist, samples.css.df$Sample_Time)
plot(samples.disper)
samples.permdisp <- permutest(samples.disper, permutations = 9999)
samples.permdisp #p = 0.11


#check for effect of animal in the model
adonis2(samples.dist ~ Body_Site+Sample_Time, by = "margin", samples.css.df)

adonis2(samples.dist ~ Body_Site+Animal, by = "margin", samples.css.df)

samples.disper <-  betadisper(samples.dist, samples.css.df$Animal)
plot(samples.disper)
samples.permdisp <- permutest(samples.disper, permutations = 9999)
samples.permdisp

plot_ordination(samples.css, samples.ord, type = "Animal", color = "Animal")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= Animal), alpha = 0.3, lty= 2, size= 1)+
  scale_colour_manual(values = randomColor(20))+
  scale_fill_manual(values = randomColor(20))+
  #(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

# Make plots with no LA ---------------------------------------------------

no_la <- subset_samples(samples.css, Body_Site %in% c("Oral", "Rumen", "Abomasum","Duodenum","Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon", "Fecal"))
no_la_samples.dist <- gunifrac(no_la)
no_la.ord <- ordinate(no_la, method = "NMDS", distance = no_la_samples.dist)
no_la.df <- as.data.frame(as(sample_data(no_la),"matrix"))

#check and see what effect of animal is 
adonis2(no_la_samples.dist~ Body_Site+Animal, by = "margin", no_la.df)

#ordination plot no LA
plot_ordination(no_la, no_la.ord, type = "samples", color = "Body_Site")+
  theme_bw()+
  labs(title = "", x= "NMDS1", y= "NMDS2")+
  geom_point(size= 5, shape= 18)+
  #geom_text(label = samples.css@sam_data$Sample_names)+
  stat_ellipse(geom = "polygon", aes(fill= Body_Site), alpha = 0.3, lty= 2, size= 1)+
  scale_colour_manual(values = no_la_palette)+
  scale_fill_manual(values = no_la_palette)+
  #(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24))

#Dendro no la
#cluster
no_la_samples.hclust <- hclust(no_la_samples.dist, method = "ward.D2")
no_la_samples.dendro.plot <- plot(no_la_samples.hclust)
no_la_samples.dendro <- as.dendrogram(no_la_samples.hclust)
no_la_samples.dendro.data <- dendro_data(no_la_samples.dendro, type = "rectangle")
no_la_metadata_for_dendro.samples <- as_tibble(no_la@sam_data)
no_la_samples.dendro.data$labels <- no_la_samples.dendro.data$labels %>% 
  left_join(no_la_metadata_for_dendro.samples, by = c("label" = "Sample_names"))
no_la_samples.dendro.data$labels
no_la_samples.dendro.order <- no_la_samples.dendro.data$labels$label
no_la_samples.dendro.order

#make sample time a factor
no_la_samples.dendro.data$labels$Sample_Time <- as.factor(no_la_samples.dendro.data$labels$Sample_Time)

no_la_dendro <- ggplot(no_la_samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = no_la_samples.dendro.data$labels, aes(x=x, y=y, colour = Body_Site), size = 9, shape = 15, position = position_nudge())+
  geom_point(data = no_la_samples.dendro.data$labels, aes(x=x,y=y, colour = Sample_Time), size = 9, shape = 15, position = position_nudge(y=-.2)) +
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_x_discrete(limits = no_la_samples.dendro.order, expand = c(0.01,0,0,0))+
  scale_colour_manual(values = no_la_paired_pallette)+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
no_la_dendro

#calculate
no_la_ra_samples <- prune_taxa(taxa_sums(no_la) > 0, no_la)
no_la_ra_samples#74193 taxa, 159 samples
no_la_samples_phy <- tax_glom(no_la_ra_samples, taxrank = "Phylum", NArm = F)

#look at phylum
no_la_samples_phy #49 taxa 159 samples
no_la_samples_phy_filt <- merge_low_abundance(no_la_samples_phy, threshold = 0.1)
no_la_samples_phy_filt #13 taxa
no_la_samples_phy_melt <- psmelt(no_la_samples_phy_filt)


#calcualte abundance
no_la_datphy <- no_la_samples_phy_melt %>% group_by(Phylum) %>% 
  summarize(avg=mean(Abundance, na.rm=TRUE),
            sd=sd(Abundance, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))


#RA plot
no_la_ra_phy <- ggplot(no_la_samples_phy_melt, aes(x= Sample, y= Abundance, fill = Phylum)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = no_la_samples.dendro.order, expand = c(0.01,0,0,0)) +
  scale_fill_manual(values = phylum_palette ) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
no_la_ra_phy

#merge plots
ra_plot_new <- plot_grid(no_la_dendro, no_la_ra_phy, align = "hv", axis = "btrl", ncol = 1)
ra_plot_new

family <- plot_grid(firm_plot, bac_plot, act_plot, align = "hv", axis = "btrl", ncol = 1)
family

co_fam <- plot_grid(firm_co_plot, bac_co_plot, act_co_plot,align = "hv", axis = "btrl", ncol = 1 )
co_fam


# Relative Abundance ------------------------------------------------------
#calculations
rel_abund <- transform_sample_counts(samples.css, function(x) {x/sum(x)} * 100)
ra_family <- tax_glom(rel_abund, taxrank = "Family", NArm = F)
ra_family_filt <- merge_low_abundance(ra_family, threshold = 0.001) # 551 taxa 166 samples
ra_family_filt #163 samples and 226 taxa

ra_family_palette <- distinctColorPalette(192)
write.csv(ra_family_palette, "famil_palette.csv")
write.csv(tax_table(ra_family_filt), "familytax.csv")
ra_genus <- tax_glom(rel_abund, taxrank = "Genus", NArm = F)

#test for a difference
microbiome_df = as.data.frame(as(sample_data(samples.css),"matrix"))
adonis_microbiome <- adonis2(samples.dist ~Body_Site, microbiome_df)
adonis_microbiome


#pairwise
adonis_microbiome_body_site <- pairwise.adonis2(samples.dist ~ Body_Site, microbiome_df)
# view the results
adonis_microbiome_body_site

#check dispersion
permutest(betadisper(samples.dist, microbiome_df$Body_Site), pairwise = TRUE)

##kill group
adonis2(samples.dist ~ Sample_Time, microbiome_df)#p = 0.001
permutest(betadisper(samples.dist, microbiome_df$Sample_Time), pairwise = T)#no diff p = 0.123


# Check individual sample RA for odd things by body site ------------------

#orals
ra_oral <- subset_samples(rel_abund, Body_Site=="Oral")
ra_oral <- prune_taxa(taxa_sums(ra_oral) > 0, ra_oral)
ra_oral #18212 taxa, 18 samples
oral_phy <- tax_glom(ra_oral, taxrank = "Phylum", NArm = F)
oral_class <- tax_glom(ra_oral, taxrank = "Class", NArm = F)
oral_family <- tax_glom(ra_oral, taxrank = "Family", NArm = F)
oral_genus <- tax_glom(ra_oral, taxrank = "Genus", NArm = F)

oral_family_filt <- merge_low_abundance(oral_family, threshold = 0.1)
oral_family_filt #68 taxa 
oral_family_melt <- psmelt(oral_family_filt)
oral_family_melt

oralPlot <- ggplot(oral_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = oral_family_palette ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))

oralPlot

#rumen
ra_rumen <- subset_samples(rel_abund, Body_Site=="Rumen")
ra_rumen <- prune_taxa(taxa_sums(ra_rumen) > 0, ra_rumen)
ra_rumen #14573 taxa, 18 samples
rumen_phy <- tax_glom(ra_rumen, taxrank = "Phylum", NArm = F)
rumen_class <- tax_glom(ra_rumen, taxrank = "Class", NArm = F)
rumen_family <- tax_glom(ra_rumen, taxrank = "Family", NArm = F)
rumen_genus <- tax_glom(ra_rumen, taxrank = "Genus", NArm = F)

rumen_family_filt <- merge_low_abundance(rumen_family, threshold = 0.1)
rumen_family_filt #51 taxa 
rumen_family_melt <- psmelt(rumen_family_filt)

rumenPlot <- ggplot(rumen_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = rumen_family_palette ) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))

rumenPlot

#abo
ra_abo <- subset_samples(rel_abund, Body_Site=="Abomasum")
ra_abo <- prune_taxa(taxa_sums(ra_abo) > 0, ra_abo)
ra_abo #15120 taxa, 18 samples
abo_phy <- tax_glom(ra_abo, taxrank = "Phylum", NArm = F)
abo_class <- tax_glom(ra_abo, taxrank = "Class", NArm = F)
abo_family <- tax_glom(ra_abo, taxrank = "Family", NArm = F)
abo_genus <- tax_glom(ra_abo, taxrank = "Genus", NArm = F)

abo_family_filt <- merge_low_abundance(abo_family, threshold = 0.1)
abo_family_filt #48 taxa 
abo_family_melt <- psmelt(abo_family_filt)

aboPlot <- ggplot(abo_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))

aboPlot

#duo
ra_duo <- subset_samples(rel_abund, Body_Site=="Duodenum")
ra_duo <- prune_taxa(taxa_sums(ra_duo) > 0, ra_duo)
ra_duo #15332 taxa, 11 samples
duo_phy <- tax_glom(ra_duo, taxrank = "Phylum", NArm = F)
duo_class <- tax_glom(ra_duo, taxrank = "Class", NArm = F)
duo_family <- tax_glom(ra_duo, taxrank = "Family", NArm = F)
duo_genus <- tax_glom(ra_duo, taxrank = "Genus", NArm = F)

duo_family_filt <- merge_low_abundance(duo_family, threshold = 0.1)
duo_family_filt #36 taxa 
duo_family_melt <- psmelt(duo_family_filt)

duoPlot <- ggplot(duo_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))

duoPlot

#jej
ra_jej <- subset_samples(rel_abund, Body_Site=="Jejunum")
ra_jej <- prune_taxa(taxa_sums(ra_jej) > 0, ra_jej)
ra_jej # 17185 taxa, 16 samples
jej_phy <- tax_glom(ra_jej, taxrank = "Phylum", NArm = F)
jej_class <- tax_glom(ra_jej, taxrank = "Class", NArm = F)
jej_family <- tax_glom(ra_jej, taxrank = "Family", NArm = F)
jej_genus <- tax_glom(ra_jej, taxrank = "Genus", NArm = F)

jej_family_filt <- merge_low_abundance(jej_family, threshold = 0.1)
jej_family_filt #42 taxa 
jej_family_melt <- psmelt(jej_family_filt)

jejPlot <- ggplot(jej_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))

jejPlot


#ile
ra_ile <- subset_samples(rel_abund, Body_Site=="Ileum")
ra_ile <- prune_taxa(taxa_sums(ra_ile) > 0, ra_ile)
ra_ile #17491 taxa, 13 samples
ile_phy <- tax_glom(ra_ile, taxrank = "Phylum", NArm = F)
ile_class <- tax_glom(ra_ile, taxrank = "Class", NArm = F)
ile_family <- tax_glom(ra_ile, taxrank = "Family", NArm = F)
ile_genus <- tax_glom(ra_ile, taxrank = "Genus", NArm = F)

ile_family_filt <- merge_low_abundance(ile_family, threshold = 0.1)
ile_family_filt #52 taxa 
ile_family_melt <- psmelt(ile_family_filt)

ilePlot <- ggplot(ile_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))

ilePlot

#cec
ra_cec <- subset_samples(rel_abund, Body_Site=="Cecum")
ra_cec <- prune_taxa(taxa_sums(ra_cec) > 0, ra_cec)
ra_cec #17083 taxa, 17 samples
cec_phy <- tax_glom(ra_cec, taxrank = "Phylum", NArm = F)
cec_class <- tax_glom(ra_cec, taxrank = "Class", NArm = F)
cec_family <- tax_glom(ra_cec, taxrank = "Family", NArm = F)
cec_genus <- tax_glom(ra_cec, taxrank = "Genus", NArm = F)

cec_family_filt <- merge_low_abundance(cec_family, threshold = 0.1)
cec_family_filt #51 taxa 
cec_family_melt <- psmelt(cec_family_filt)

cecPlot <- ggplot(cec_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))

cecPlot

#sprC
ra_sprC <- subset_samples(rel_abund, Body_Site=="Spiral_Colon")
ra_sprC <- prune_taxa(taxa_sums(ra_sprC) > 0, ra_sprC)
ra_sprC #19510 taxa, 15 samples
sprC_phy <- tax_glom(ra_sprC, taxrank = "Phylum", NArm = F)
sprC_class <- tax_glom(ra_sprC, taxrank = "Class", NArm = F)
sprC_family <- tax_glom(ra_sprC, taxrank = "Family", NArm = F)
sprcfam_melt <- psmelt(sprC_family)
sprC_genus <- tax_glom(ra_sprC, taxrank = "Genus", NArm = F)

sprC_family_filt <- merge_low_abundance(sprC_family, threshold = 0.1)
sprC_family_filt #52 taxa 
sprC_family_melt <- psmelt(sprC_family_filt)
sprC_family_melt


sprCPlot <- ggplot(sprC_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(318) ) +
  theme(legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
sprCPlot

#disC
ra_disC <- subset_samples(rel_abund, Body_Site=="Distal_Colon")
ra_disC <- prune_taxa(taxa_sums(ra_disC) > 0, ra_disC)
ra_disC #22810 taxa, 18 samples
disC_phy <- tax_glom(ra_disC, taxrank = "Phylum", NArm = F)
disC_class <- tax_glom(ra_disC, taxrank = "Class", NArm = F)
disC_family <- tax_glom(ra_disC, taxrank = "Family", NArm = F)
disC_genus <- tax_glom(ra_disC, taxrank = "Genus", NArm = F)

disC_family_filt <- merge_low_abundance(disC_family, threshold = 0.1)
disC_family_filt #55 taxa 
disC_family_melt <- psmelt(disC_family_filt)

disCPlot <- ggplot(disC_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))

disCPlot


#fecal
ra_fecal <- subset_samples(rel_abund, Body_Site=="Fecal")
ra_fecal <- prune_taxa(taxa_sums(ra_fecal) > 0, ra_fecal)
ra_fecal #17860 taxa, 15 samples
fecal_phy <- tax_glom(ra_fecal, taxrank = "Phylum", NArm = F)
fecal_class <- tax_glom(ra_fecal, taxrank = "Class", NArm = F)
fecal_family <- tax_glom(ra_fecal, taxrank = "Family", NArm = F)
fecal_genus <- tax_glom(ra_fecal, taxrank = "Genus", NArm = F)

fecal_family_filt <- merge_low_abundance(fecal_family, threshold = 0.1)
fecal_family_filt #49 taxa 
fecal_family_melt <- psmelt(fecal_family_filt)

fecalPlot <- ggplot(fecal_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))

fecalPlot

#LA
ra_la <- subset_samples(rel_abund, Body_Site=="Liver_Abscess")
ra_la <- prune_taxa(taxa_sums(ra_la) > 0, ra_la)
ra_la #992 taxa, 4 samples
la_phy <- tax_glom(ra_la, taxrank = "Phylum", NArm = F)
la_class <- tax_glom(ra_la, taxrank = "Class", NArm = F)
la_family <- tax_glom(ra_la, taxrank = "Family", NArm = F)
la_genus <- tax_glom(ra_la, taxrank = "Genus", NArm = F)

la_genus_melt <- psmelt(la_genus)

la_family_filt <- merge_low_abundance(la_family, threshold = 0.1)
la_family_filt #5 taxa 
la_family_melt <- psmelt(la_family_filt)

laPlot <- ggplot(la_genus_melt, aes(x= Sample, y= Abundance, fill = Genus)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = randomColor(113) ) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))

laPlot


# Relative Abundance along the GIT ----------------------------------------
ra_samples <- prune_taxa(taxa_sums(samples.css) > 0, samples.css)
ra_samples#74565 taxa, 163 samples
samples_phy <- tax_glom(ra_samples, taxrank = "Phylum", NArm = F)
samples_class <- tax_glom(ra_samples, taxrank = "Class", NArm = F)
samples_genus <- tax_glom(ra_samples, taxrank = "Genus", NArm = F)


#cluster
samples.hclust <- hclust(samples.dist, method = "ward.D2")
samples.dendro.plot <- plot(samples.hclust)
samples.dendro <- as.dendrogram(samples.hclust)
samples.dendro.data <- dendro_data(samples.dendro, type = "rectangle")
metadata_for_dendro.samples <- as_tibble(samples.css@sam_data)
samples.dendro.data$labels <- samples.dendro.data$labels %>% 
  left_join(metadata_for_dendro.samples, by = c("label" = "Sample_names"))
samples.dendro.data$labels
samples.dendro.order <- samples.dendro.data$labels$label
samples.dendro.order

#look at phylum
samples_phy #49 taxa 163 samples
samples_phy_filt <- merge_low_abundance(samples_phy, threshold = 0.1)
samples_phy_filt #14 taxa
samples_phy_melt <- psmelt(samples_phy_filt)
samples_phy_melt
View(samples_phy_melt)

#calcualte abundance
datphy <- samples_phy_melt %>% group_by(Phylum) %>% 
  summarize(avg=mean(Abundance, na.rm=TRUE),
            sd=sd(Abundance, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))

write.csv(datphy, "phylum abundance.csv")

#build custom pallete
phy_palette <- distinctColorPalette(14)
write.csv(phy_palette, "Phy_palette.csv")
write.csv(tax_table(samples_phy_filt), "RA_Phy_filt.csv")



#look at family
ra_family#651 taxa 163 samples
samples_family_filt <- merge_low_abundance(ra_family, threshold = 0.1)
samples_family_filt #120 taxa 
#write.csv(tax_table(samples_family_filt), "allsamplesfilt.csv")
samples_family_melt <- psmelt(samples_family_filt)
samples_family_melt

#calculate abundance at family
datfamily <- samples_family_melt %>% group_by(Body_Site,Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))

samples.dendro.data$labels$Sample_Time <- as.factor(samples.dendro.data$labels$Sample_Time)

#Plots
microbiome_dendro_plot <- ggplot(samples.dendro.data$segments) +
  theme_minimal() +
  labs(y= "") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data = samples.dendro.data$labels, aes(x=x, y=y, colour = Body_Site), size = 9, shape = 15, position = position_nudge())+
  geom_point(data = samples.dendro.data$labels, aes(x=x,y=y, colour = Sample_Time), size = 9, shape = 15, position = position_nudge(y=-.2)) +
  #geom_text(aes(x=x, y=y, label=label), data = samples.dendro.data$labels)+
  scale_colour_manual(values = paired_palette)+
  theme(plot.title = element_text(size = 20),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(size = 0.75, colour = "black"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.7))
microbiome_dendro_plot

##joined plot
body_site_col <- data.frame(Body_Site = c(""))
time_col <- data.frame(Sample_Time = c(""))
samples.dendro.data$labels <- cbind(samples.dendro.data$labels, body_site_col, time_col)
samples.dendro.data

joined_dendro_fig<-ggplot(samples.dendro.data$segments) + theme_bw() +
  labs(y= "Ward's Distance") +
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend), linewidth =.5, lineend = "round", linejoin = "round") +
  geom_point(data = samples.dendro.data$labels, aes(x=x,y=y, color = Body_Site), 
             size =8, position = position_nudge(x=-0.08), shape = 15) +
  geom_text(data = samples.dendro.data$labels,
            aes(x, y, label = Body_Site),  # replace label_variable with your actual label variable
            size = 4, position = position_nudge(y = -0.05)) +
  scale_colour_manual(values = bodysite_palette) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size = 0.75),
        axis.title.y = element_text(size = 36),
        axis.text.y = element_text(size = 20, colour = "black"))
joined_dendro_fig

ra_phy <- ggplot(samples_phy_melt, aes(x= Sample, y= Abundance, fill = Phylum)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0)) +
  scale_fill_manual(values = phylum_palette ) +
  theme(legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
ra_phy


allsamples_ra <- ggplot(samples_family_melt, aes(x= Sample, y= Abundance, fill = Family)) +
  theme_minimal() + #facet_wrap(~Body_Site, scales = "free_x")+
  labs(y= "Relative Abundance (%)") +
  geom_bar(stat = "identity", colour = "black") +
  scale_x_discrete(limits = samples.dendro.order, expand = c(0.035,0,0,0)) +
  scale_fill_manual(values = randomColor(66) ) +
  theme(legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(size = 0.75, colour = "black"),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", size = 0.7))
allsamples_ra


plot_grid(microbiome_dendro_plot, allsamples_ra, align = "v", axis = "btrl", ncol = 1)


# Compare Firmicutes vs Bacteroidetes -------------------------------------
#Firmicutes
firmicutes_sub <- subset_taxa(rel_abund, Phylum=="Firmicutes")

firmicutes_phyla <- tax_glom(firmicutes_sub, taxrank = "Phylum", NArm = F) %>%
  psmelt()
firmicutes_family <- tax_glom(firmicutes_sub, taxrank = "Family", NArm = F) %>%
  psmelt() #118

#make a custom palette
firmicutes_family1 <- tax_glom(firmicutes_sub, taxrank = "Family", NArm = F)
firm_pallete <- randomColor(118)
write.csv(firm_pallete, "RA_firm_palette.csv")
write.csv(tax_table(firmicutes_family1), "RA_firm_filt.csv")

firm_plot <- ggplot(firmicutes_family, aes(x= Body_Site, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y= "", x = "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(firmicutes_phyla, mapping= aes(x= Body_Site, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = firm_palette) +
  scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"),
                   labels = c("","","","","","","","","","","")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 28, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))
firm_plot

kruskal.test(firmicutes_phyla$Abundance, firmicutes_phyla$Body_Site, block = firmicutes_phyla$Animal)
firm_anova <- pairwise.wilcox.test(firmicutes_phyla$Abundance, firmicutes_phyla$Body_Site, block = firmicutes_phyla$Animal, p.adjust.method = "BH")
firm_anova

#By harvest group
firmicutes_phyla$Sample_Time <- as.factor(firmicutes_phyla$Sample_Time)
pairwise.wilcox.test(firmicutes_phyla$Abundance, firmicutes_phyla$Sample_Time, p.adjust.method = "BH")

firmicutes_family$Sample_Time <- factor(firmicutes_family$Sample_Time)

firm_co_plot <- ggplot(firmicutes_family, aes(x= Sample_Time, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y= "", x = "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(firmicutes_phyla, mapping= aes(x= Sample_Time, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = firm_palette) +
  #scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"),
  #                 labels = c("","","","","","","","","","","")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size = 0.9, colour = "black"))
firm_co_plot

#superscript
firm_pvalues <- firm_anova$p.value %>% 
  fullPTable()
multcompLetters(firm_pvalues)

datfirm <- firmicutes_family %>% group_by(Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datfirm,"RA Firm.csv")
#Bacteroidetes
bacteroidetes_sub <- subset_taxa(rel_abund, Phylum=="Bacteroidota")

bacteroidetes_phyla <- tax_glom(bacteroidetes_sub, taxrank = "Phylum", NArm = F) %>%
  psmelt()
bacteroidetes_family <- tax_glom(bacteroidetes_sub, taxrank = "Family", NArm = F) %>%
  psmelt()

#build palette
bacteroidetes_family1 <- tax_glom(bacteroidetes_sub, taxrank = "Family", NArm = F)
bac_pallete <- randomColor(57)
write.csv(bac_pallete, "RA_bac_palette.csv")
write.csv(tax_table(bacteroidetes_family1), "RA_bac_filt.csv")

bac_plot <- ggplot(bacteroidetes_family, aes(x= Body_Site, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y="", x= "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(bacteroidetes_phyla, mapping= aes(x= Body_Site, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = bac_palette) +
  scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"),
                   labels = c("","","","","","","","","","","")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black"),
        axis.ticks = element_line(size = 0.9, colour = "black"))
bac_plot
bac_anova <- pairwise.wilcox.test(bacteroidetes_phyla$Abundance, bacteroidetes_phyla$Body_Site, p.adjust.method = "BH")
#superscript
bac_pvalues <- bac_anova$p.value %>% 
  fullPTable()
multcompLetters(bac_pvalues)

datbacteroidetes <- bacteroidetes_family %>% group_by(Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datbacteroidetes, "RA Bacter.csv")

#by harvest group
bacteroidetes_phyla$Sample_Time <- as.factor(bacteroidetes_phyla$Sample_Time)
pairwise.wilcox.test(bacteroidetes_phyla$Abundance, bacteroidetes_phyla$Sample_Time, p.adjust.method = "BH")

bacteroidetes_family$Sample_Time <- factor(bacteroidetes_family$Sample_Time)

bac_co_plot <- ggplot(bacteroidetes_family, aes(x= Sample_Time, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y= "", x = "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(bacteroidetes_family, mapping= aes(x= Sample_Time, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = bac_palette) +
  #scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"),
  #                 labels = c("","","","","","","","","","","")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size = 0.9, colour = "black"))
bac_co_plot

#Actinobacteriota
actino_sub <- subset_taxa(rel_abund, Phylum=="Actinobacteriota")

actino_phy <- tax_glom(actino_sub, taxrank = "Phylum", NArm = F) %>% 
  psmelt()
actino_phy
actino_family <- tax_glom(actino_sub, taxrank = "Family", NArm = F) %>% 
  psmelt()

#make palette
actino_pallete <- randomColor(73)
actino_family1 <- tax_glom(actino_sub, taxrank = "Family", NArm = F)
write.csv(actino_pallete, "RA_actino_palette.csv")
write.csv(tax_table(actino_family1), "RA_actino_filt.csv")


act_plot <- ggplot(actino_family, aes(x= Body_Site, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y="", x= "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(actino_phy, mapping= aes(x= Body_Site, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = actino_palette) +
  scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"),
                   labels = c("", "", "", "", "", "", "", "", "","", "")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 28, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, colour = "black"),
        axis.ticks = element_line(size = 0.9, colour = "black"))
act_plot
act_anova <- pairwise.wilcox.test(actino_phy$Abundance, actino_phy$Body_Site, p.adjust.method = "BH")
act_anova
#superscript
act_pvalues <- act_anova$p.value %>% 
  fullPTable()
multcompLetters(act_pvalues)

datactino <- actino_family %>% group_by(Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datactino, "RA Actino.csv")

#by harvest group
actino_phy$Sample_Time <- as.factor(actino_phy$Sample_Time)
pairwise.wilcox.test(actino_phy$Abundance, actino_phy$Sample_Time, p.adjust.method = "BH")

actino_family$Sample_Time <- as.factor(actino_family$Sample_Time)

act_co_plot <- ggplot(actino_family, aes(x= Sample_Time, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y="", x= "") +
  geom_bar(aes(fill= Family), stat = "summary", colour = "black") +
  geom_errorbar(actino_phy, mapping= aes(x= Sample_Time, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = actino_palette) +
  #scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal", "Liver_Abscess"),
                   #labels = c("", "", "", "", "", "", "", "", "","", "")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  ylim(0,100)+
  theme(legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 28, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_line(size = 0.9, colour = "black"))
act_co_plot
datactino1 <- actino_family %>% group_by(Sample_Time,Family) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
# Search For Important Genera ----------------------------------------------
#Fuso
fuso <- subset_taxa(ra_genus, Genus=="Fusobacterium")
fuso
fuso.melt <- psmelt(fuso)


fusoboxplot <- ggplot(fuso.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = bodysite_palette )+
  coord_cartesian(ylim = c(0,5))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
fusoboxplot

kruskal.test(fuso.melt$Abundance,fuso.melt$Body_Site)#p = 2.0733-14
fuso_anova <- pairwise.wilcox.test(fuso.melt$Abundance,fuso.melt$Body_Site, p.adjust.m = "BH")
#superscript
fuso_pvalues <- fuso_anova$p.value %>% 
  fullPTable()
multcompLetters(fuso_pvalues)
#write abundance values
datfuso <- fuso.melt %>% group_by(Body_Site) %>% 
  summarize(avg=mean(Abundance, na.rm=TRUE),
            sd=sd(Abundance, na.rm = T),
            lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
            upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n()))

write.csv(datfuso, "RA Fuso.csv")

#bact family 
bact <- subset_taxa(ra_family, Family=="Bacterroidaceae")
bact <- psmelt(bact)
bact
ggplot(bact, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "Abundance", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())

#Bacteroides
bacteroides <- subset_taxa(ra_genus, Genus=="Bacteroides")
bacteroides
bac.melt <- psmelt(bacteroides)
bac_boxplot <- ggplot(bac.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  scale_fill_manual(values = bodysite_palette )+
  coord_cartesian(ylim = c(0,20))+
  geom_point()+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
bac_boxplot
kruskal.test(bac.melt$Abundance,bac.melt$Body_Site) #p= 2.2e-16
bac_g_anova <- pairwise.wilcox.test(bac.melt$Abundance, bac.melt$Body_Site, p.adjust.method = "BH")

#superscript
bac_g_pvalues <- bac_g_anova$p.value %>% 
  fullPTable()
multcompLetters(bac_g_pvalues)

#write abundacne
datbac <- bac.melt %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(datbac, "RA Bact.csv")

#Bifido
bifido <- subset_taxa(ra_genus, Genus=="Bifidobacterium")
bifido
bifido.melt <- psmelt(bifido)
bifido_plot <- ggplot(bifido.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  scale_y_continuous(breaks = c(0, 5, 10, 15),
                     labels = c("0", "5", "10", "15"))+
  scale_fill_manual(values = bodysite_palette)+
  geom_point()+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
bifido_plot
kruskal.test(bifido.melt$Abundance,bifido.melt$Body_Site)
bif_anova <- pairwise.wilcox.test(bifido.melt$Abundance,bifido.melt$Body_Site, p.adjust.method = "BH")
#superscript
bif_pvalues <- bif_anova$p.value %>% 
  fullPTable()
multcompLetters(bif_pvalues)

datbifido <- bifido.melt %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=T),
            sd(Abundance, na.rm = T))

write.csv(datbifido, "Bifido RA.csv")

#Porphormonas
porph <- subset_taxa(ra_genus, Genus=="Porphyromonas")
porph
porph.melt <- psmelt(porph)
porph_boxplot <- ggplot(porph.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values =  bodysite_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
porph_boxplot

kruskal.test(porph.melt$Abundance,porph.melt$Body_Site)
pairwise.wilcox.test(porph.melt$Abundance,porph.melt$Body_Site, p.adjust.method = "BH")

#superscript
porph_pvalues <- porph_anova$p.value %>% 
  fullPTable()
multcompLetters(porph_pvalues)

datporph <- porph.melt %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=T),
            sd(Abundance,na.rm = T))
write.csv(datporph, "RA Porph.csv")

#truperella
trup <- subset_taxa(ra_genus, Genus=="Trueperella")
trup
trup.melt <- psmelt(trup)
trup_boxplot <- ggplot(trup.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "", x= "", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values =  bodysite_palette)+
  coord_cartesian(ylim = c(0,5))+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
trup_boxplot

kruskal.test(trup.melt$Abundance,trup.melt$Body_Site)
trup_anova <- pairwise.wilcox.test(trup.melt$Abundance,trup.melt$Body_Site, p.adjust.method = "BH")
#superscript
trup_pvalues <- trup_anova$p.value %>% 
  fullPTable()
multcompLetters(trup_pvalues)

dattrup <- trup.melt %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))
write.csv(dattrup, "RA Trup.csv")



#entercoccus
enter <- subset_taxa(ra_genus, Genus=="Enterococcus")
enter
enter.melt <- psmelt(enter)
enterplot <- ggplot(enter.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "Relative Abundance (%)", x="", title = "Enterococcus Abundance") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal", "Liver_Abscess"))+
  ylim(0,1)+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = bodysite_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
enterplot
#test for difference
kruskal.test(enter.melt$Abundance,enter.melt$Body_Site)
enter_anova <- pairwise.wilcox.test(enter.melt$Abundance,enter.melt$Body_Site,p.adjust.method = "BH")

#superscript
enter_pvalues <- enter_anova$p.value %>% 
  fullPTable()
multcompLetters(enter_pvalues)

datenter <- enter.melt %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))


#Lactobacillus
lact <- subset_taxa(ra_genus, Genus=="Lactobacillus")
lact
lact.melt <- psmelt(lact)
lactplot <- ggplot(lact.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_y_continuous(breaks = c(0,0.05,0.10,.15),
                     labels = c("0","0.05","0.10","0.15"))+
  scale_fill_manual(values = bodysite_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
lactplot
#test for difference
kruskal.test(lact.melt$Abundance,lact.melt$Body_Site)#p=0.09 trend

datlact <- lact.melt %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))


#salmonella
sal <- subset_taxa(ra_genus, Genus=="Salmonella")


#enterobacteriaceae
entero_sub <- subset_taxa(rel_abund, Family=="Enterobacteriaceae")

enter_genus <- tax_glom(entero_sub, taxrank = "Genus", NArm = F) %>% 
  psmelt()

entero_family <- tax_glom(entero_sub, taxrank = "Family", NArm = F) %>% 
  psmelt()

#entero <- subset_taxa(ra_family, Family=="Enterobacteriaceae")
#entero
#entero.melt <- psmelt(entero)
ggplot(enter_genus, aes(x= Body_Site, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y="", x= "") +
  geom_bar(aes(fill= Genus), stat = "summary", colour = "black") +
  geom_errorbar(entero_family, mapping= aes(x= Body_Site, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = randomColor(31)) +
  scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal", "Liver_Abscess"),
                   labels = c("", "", "", "", "", "", "", "", "","", "")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
    theme(#legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 28, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))


#test for difference
kruskal.test(entero.melt$Abundance,entero.melt$Body_Site)
entero_anova <- pairwise.wilcox.test(entero.melt$Abundance, entero.melt$Body_Site, p.adjust.method = "BH")

#superscript
entero_pvalues <- entero_anova$p.value %>% 
  fullPTable()
multcompLetters(entero_pvalues)


datentero <- enter_genus%>% group_by(Genus) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))


#Escherica
esch <- subset_taxa(enter_genus, Genus=="Escherichia-Shigella")
esch
esch.melt <- psmelt(esch)
eschplot <- ggplot(esch.melt, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "Relative Abundance (%)", x="", title = "Escherichia-Shigella Abundance") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal", "Liver_Abscess"))+
  geom_boxplot() +
  geom_point()+
  scale_fill_manual(values = bodysite_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
eschplot

esch_anova <- pairwise.wilcox.test(esch.melt$Abundance,esch.melt$Body_Site,p.adjust.method = "BH")

#superscript
esch_pvalues <- esch_anova$p.value %>% 
  fullPTable()
multcompLetters(esch_pvalues)

datesch <- esch.melt %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=TRUE),
            sd(Abundance, na.rm = T))

#unclassified enterobacteriacea
un <- subset_taxa(enter_genus, Genus=="Unclasssified Enterobacteriaceae")


#bacillus
bacill <- subset_taxa(ra_genus, Genus=="Bacillus") %>% 
  psmelt()


bacill_boxplot <- ggplot(bacill, aes(x= Body_Site, y= Abundance, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1),
                     labels = c("0","0.25","0.50","0.75","1"))+
  scale_fill_manual(values = bodysite_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
bacill_boxplot

kruskal.test(bacill$Abundance, bacill$Body_Site)
bacil_anova <- pairwise.wilcox.test(bacill$Abundance,bacill$Body_Site, p.adjust.method = "BH")
#superscript
bacil_pvalues <- bacil_anova$p.value %>% 
  fullPTable()
multcompLetters(bacil_pvalues)

datbacil <- bacill %>% group_by(Body_Site) %>% 
  summarize(mean(Abundance, na.rm=T),
            sd(Abundance, na.rm = T))

write.csv(datbacil, "Bacillus RA.csv")

#gama proteobacteria

gama <- subset_taxa(rel_abund, Class =="Gammaproteobacteria") 
  psmelt()

gama_class <- tax_glom(gama, taxrank = "Class", NArm = F) %>% 
  psmelt()
gama_genus <- tax_glom(gama, taxrank = "Genus", NArm = F) %>% 
  psmelt()
gama_family <- tax_glom(gama, taxrank = "Family", NArm = F) %>% 
  psmelt()

gamma_palette <- randomColor(211)

ggplot(gama_genus, aes(x= Body_Site, y= Abundance)) +
  theme_bw() + #coord_flip() +
  labs(y="", x= "") +
  geom_bar(aes(fill= Genus), stat = "summary", colour = "black") +
  geom_errorbar(gama_class, mapping= aes(x= Body_Site, y= Abundance), stat = "summary", width = 0.55, size = 0.6) +
  scale_fill_manual(values = gamma_palette) +
  scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal", "Liver_Abscess"),
                   labels = c("", "", "", "", "", "", "", "", "","", "")) +
  scale_y_continuous(expand = c(0.002,0,0.1,0)) +
  theme(legend.position = "none",
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 28, vjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 28, colour = "black"),
    axis.ticks = element_line(size = 0.9, colour = "black"))
datgamma <- gama_genus%>% group_by(Genus) %>% 
  summarize(mean(Abundance, na.rm=T),
            sd(Abundance, na.rm = T))

# Combining plots ---------------------------------------------------------

plot_grid(fusoboxplot,trup_boxplot, bac_boxplot, porph_boxplot, align = "v", axis = "btrl", ncol = 2)

plot_grid(firmtobac_boxplot, bacill_boxplot, bifido_plot,lactplot, align = "v", axis = "btrl", ncol = 2)

plot_grid(enterplot, enteroplot, align = "h", axis = "btrl", ncol = 2)


# Firmicutes : bacteroides ratio ------------------------------------------
#pop out the excel
write.csv(bacteroidetes_phyla, "Bac phyla.csv")
write.csv(firmicutes_phyla, "Firm phyla.csv")

#read in new excel
ratio <- read_excel("Ratio_doc.xlsx")

#histogram of ratio
histogram(ratio$firm_ratio)

firmtobac_boxplot <- ggplot(ratio, aes(x= Body_Site, y= firm_ratio, fill = Body_Site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_y_log10()+
  scale_fill_manual(values = bodysite_palette)+
  theme(legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
firmtobac_boxplot
#test for differences
kruskal.test(ratio$firm_ratio,ratio$Body_Site)
pairtest <- pairwise.wilcox.test(ratio$firm_ratio, ratio$Body_Site, p.adjust.method = "BH")

ratio_pvalues <- pairtest$p.value %>% 
  fullPTable()
multcompLetters(ratio_pvalues)

#write abundance values
datratio <- ratio %>% group_by(Body_Site) %>% 
  summarize(mean(firm_ratio, na.rm=TRUE),
            sd(firm_ratio, na.rm = T))
write.csv(datratio, "ratio.csv")


# Core Taxa ---------------------------------------------------------------
#transform data to check for core

core_ra <- transform(samples_genus,"compositional")
core_ra_genus <- aggregate_taxa(core_ra,"Genus")
core_ra_genus_filt <- core(core_ra_genus, detection= 0.001, prevalence =.9)
core_ra_genus_filt@tax_table

##bubble plot
#subset core taxa
core_taxa <- core_ra_genus_filt %>% psmelt() %>% 
  group_by(Body_Site, OTU) %>% 
  summarize(mean(Abundance, na.rm = T))
#paste a better name
names(core_taxa)[3] <- paste("avg")
#transform long ways
core_taxa <- core_taxa %>% spread(OTU, avg)
#transform for plot
core_taxa_long <- melt(core_taxa, id =c("Body_Site"))

test <- randomColor(6)
#plot, if transformed right variable will include all OTUs on the y axis
#adjust limits, for min and max, range for size of the dots, and breaks for the legned 
plot1 <- ggplot(core_taxa_long, aes(x = Body_Site, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, 1), range = c(1,17), breaks = c(0.001, 0.01, 0.1, 1)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")+ 
theme(legend.key=element_blank(), 
      axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
      axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
      legend.text = element_text(size = 10, face ="bold", colour ="black"), 
      legend.title = element_text(size = 12, face = "bold"), 
      panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
      legend.position = "right") +  
  scale_fill_manual(values = test, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(core_taxa_long$variable))) 


plot1

#plotb heatmap
# core with compositionals
prevalences <- seq(0.1,1,0.1)
detections <- 10^seq(log10(1e-3), log10(0.1), length = 6)

detections
detections <- trunc(detections*10^3)/10^3

p1 <- plot_core(core_ra_genus,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p1 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))


#core no LA
#subset and calculate core taxa
core_no <- subset_samples(core_ra_genus, Body_Site %in% c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon", "Fecal"))
core_no_filt <- core(core_no, detection = 0.001, prevalence = 0.9)
core_no_filt@tax_table

#subsset core taxa and calulate average at each body site
core_taxa_no <- core_no_filt %>% psmelt() %>% 
  group_by(Sample_names, OTU, Body_Site)  
  summarize(mean(Abundance, na.rm = T))
names(core_taxa_no)[4] <- paste("avg")


#transform long ways
core_taxa_no <- core_taxa_no %>% spread(Body_Site, avg)
#transform for plot
core_taxa_no_long <- melt(core_taxa_no, id =c("Sample_names","Body_Site"))
#build a palette
test2 <- randomColor(10)

#plot, if transformed right variable will include all OTUs on the y axis
#adjust limits, for min and max, range for size of the dots, and breaks for the legend
plot2 <- ggplot(core_taxa_no_long, aes(x = variable, y = value)) + 
  geom_point(aes(size = value, fill = Body_Site), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, .2), range = c(5,50), breaks = c(0.0001,0.001, 0.01, 0.1)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")+ 
  #scale_x_discrete(limits = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"),
  #                 labels = c("Oral", "Rumen", "Abomasum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal")) +
  theme_bw()+
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +  
  scale_fill_manual(values = test2, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(core_taxa_no_long$variable))) 


plot2


plot3 <- ggplot(core_taxa_no, aes(x= OTU, y = avg, fill = Body_Site))+
  geom_point(aes(),size = 5, shape = 21, alpha = .7, position = position_jitterdodge())+
  scale_fill_manual(values = no_la_palette )+
  labs(x = "", y = "Relative Abundnance (%)" )+
  scale_y_continuous(breaks = seq(0, max(core_taxa_no$avg), by = 0.05))+
  theme_bw()+
  theme(#legend.position = "none",
       panel.border = element_rect(colour = "black", size = 1),
       title = element_text(size = 28),
       axis.ticks = element_line(colour = "black", size = 0.75),
       axis.text = element_text(colour = "black", size = 12),
       axis.title = element_text(size = 24),
       axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))

plot3
View(core_taxa_no)
plot6 <- ggplot(core_taxa_no, aes(x= OTU, y = avg, fill = Body_Site))+
  geom_point(aes(),stat = "identity", position = "jitter")+
  scale_fill_manual(values = no_la_palette )+
  labs(x = "", y = "Relative Abundnance (%)" )+
  #scale_y_continuous(breaks = seq(0, max(core_taxa_no$avg), by = 0.05))+
  theme_bw()+
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))

plot6

#new bubble for avg
core_taxa_no1 <- core_no_filt %>% psmelt() %>% 
  group_by(OTU, Body_Site) %>%   
  summarize(mean(Abundance, na.rm = T),
            sd(Abundance, na.rm =T))
names(core_taxa_no1)[3] <- paste("avg")
names(core_taxa_no1)[4] <- paste("sd")
View(core_taxa_no1)



plot4 <- ggplot(core_taxa_no1, aes(x= OTU, y = avg, fill = Body_Site))+
  geom_point(aes(),position = position_jitterdodge(), size = 5, shape = 21, alpha = .7)+
  scale_fill_manual(values = no_la_palette )+
  #geom_linerange(mapping =aes(ymin=avg-sd,ymax=avg+sd),position = "identity", linetype = 1)+
  stat_summary(fun.ymax = core_taxa_no1$avg+sd, fun.ymin= core_taxa_no1$avg-sd,geom = "line",position = "identity")+
  labs(x = "", y = "Relative Abundnance (%)" )+
  scale_y_continuous(breaks = seq(0, max(core_taxa_no1$avg), by = 0.05))+
  theme_bw()+
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))
plot4


#heat map no LA
p7 <- plot_core(core_no,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p7 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))

# Starting over -----------------------------------------------------------
#new bubble for avg
core_taxa_no1 <- core_no_filt %>% psmelt() %>% 
  group_by(OTU, Body_Site) %>% 
  summarize(avg =mean(Abundance, na.rm = T),
            sd = sd(Abundance, na.rm =T))

# Calculate the confidence interval for mean abundance
core_taxa_no1 <- core_no_filt %>%
  psmelt() %>%
  group_by(OTU, Body_Site) %>%
  summarize(
    avg = mean(Abundance, na.rm = TRUE),
    sd = sd(Abundance, na.rm = TRUE),
    lower_ci = avg - qt(0.975, df = n() - 1) * sd / sqrt(n()),
    upper_ci = avg + qt(0.975, df = n() - 1) * sd / sqrt(n())
  )
View(core_taxa_no3)





#plot averages dodge and jitter
set.seed(10)
plot4 <- ggplot(core_taxa_no1, aes(x= OTU, y = avg, fill = Body_Site))+
  geom_point(aes(),position = position_jitterdodge(jitter.width = .15), size = 5, shape = 21, alpha = .7)+
  scale_fill_manual(values = no_la_palette )+
  geom_linerange(mapping =aes(ymin=avg-sd,ymax=avg+sd),position = position_jitterdodge(jitter.width = .0001), linetype = 1)+
  labs(x = "", y = "Relative Abundnance (%)" )+
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05))+
  theme_bw()+
  theme(legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))
plot4


#just dodge
plot14 <- ggplot(core_taxa_no1, aes(x= OTU, y = avg, fill = Body_Site))+
  geom_point(aes(),position = position_dodge(width = .5), size = 5, shape = 21, alpha = .7)+
  scale_fill_manual(values = no_la_palette )+
  geom_linerange(mapping =aes(ymin=lower_ci,ymax=upper_ci),position = position_dodge(width = .5), linetype = 1)+
  labs(x = "", y = "" )+
  scale_y_continuous(breaks = seq(0, 0.25, by = 0.05))+
  theme_bw()+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", size = 1),
        title = element_text(size = 28),
        axis.ticks = element_line(colour = "black", size = 0.75),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(size = 24),
        axis.text.x = element_blank())
plot14

#create individual
core_taxa_no2 <- core_no_filt %>% psmelt() %>% 
  group_by(OTU, Body_Site)   
  summarize(mean(Abundance, na.rm = T),
            sd(Abundance, na.rm =T))
names(core_taxa_no1)[3] <- paste("avg")
names(core_taxa_no1)[4] <- paste("sd")
View(core_taxa_no2)

plot7 <- ggplot(core_taxa_no2, aes(x= OTU, y = Abundance, fill = Body_Site))+
  geom_jitter(aes(),position = position_jitterdodge(jitter.width = 1, dodge.width = 1), size = 5, shape = 21, alpha = .7)+
  scale_fill_manual(values = no_la_palette )+
  #stat_summary(geom = "line",position = "identity")+
  labs(x = "", y = "Relative Abundnance (%)" )+
  scale_y_continuous(breaks = seq(0, max(core_taxa_no2$Abundance), by = 0.01))+
  scale_x_discrete("Cut")+
  theme_bw()+
  theme(legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    title = element_text(size = 28),
    axis.ticks = element_line(colour = "black", size = 0.75),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 24),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1))
plot7





# Core taxa by region -----------------------------------------------------

oral_core <- subset_samples(core_ra_genus, Body_Site=="Oral")
core_ra_oral <- transform(oral_core, "compositional")
core_oral_filt <- core(core_ra_oral, detection= 0.001, prevalence =.9)
core_oral_filt

p2 <- plot_core(core_ra_oral,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p2 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        legend.position = "none",
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text("none"),
        axis.title.x = element_text(size = 40, vjust = -0.75))

#foregut
foregut_core <- subset_samples(core_ra_genus, Body_Site %in% c("Oral", "Rumen", "Abomasum"))
foregut_ra <- transform(foregut_core, "compositional")
core_foregut_filt <- core(foregut_ra, detection= 0.001, prevalence =.9)
core_foregut_filt

p3 <- plot_core(foregut_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p3 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text(colour = "black", size = 24),
        axis.title.x = element_text(size = 40, vjust = -0.75))
#Rumen and abo
fgut_core <- subset_samples(core_ra_genus, Body_Site %in% c("Rumen", "Abomasum"))
fgut_ra <- transform(fgut_core, "compositional")
core_fgut <- core(fgut_ra, detection = 0.001, prevalence = 0.90)
core_fgut

p6 <- plot_core(fgut_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p6 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        legend.position = "none",
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text("none"),
        axis.title.x = element_text(size = 40, vjust = -0.75))


#SI
si_core <- subset_samples(core_ra_genus, Body_Site %in% c("Duodenum", "Jejunum", "Ileum"))
si_ra <- transform(si_core, "compositional")
core_si <- core(si_ra, detection = 0.001, prevalence = 0.90)
core_si

p4 <- plot_core(si_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p4 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        legend.position = "none",
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text("none"),
        axis.title.x = element_text(size = 40, vjust = -0.75))

#hindgut
hindgut_core <- subset_samples(core_ra_genus, Body_Site %in% c("Cecum", "Distal_Colon", "Spiral_Colon", "Fecal"))
hindgut_ra <- transform(hindgut_core, "compositional")
core_hindgut <- core(hindgut_ra, detection = 0.001, prevalence = 0.90)
core_hindgut

p5 <- plot_core(hindgut_ra,
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences,
                detections = detections, min.prevalence = .90)
p5 + theme_bw() +
  labs(y= "", x= "Detection Threshold \n (Relative Abundance (%))", title = "") +
  theme(legend.key.size = unit(4, "lines"),
        legend.position = "none",
        plot.title = element_text(size = 24),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.border = element_rect(colour = "black", size = 1),
        axis.ticks = element_line(size = 0.75, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 26),
        axis.text.x = element_text("none"),
        axis.title.x = element_text(size = 40, vjust = -0.75))


hindgut_bubble <- core_hindgut %>% psmelt() %>% 
  group_by(Body_Site, OTU) %>% 
  summarize(mean(Abundance, na.rm = T)) 
  
names(hindgut_bubble)[3] <- paste("avg")
hindgut_bubble  <- hindgut_bubble %>% spread(OTU, avg)
hindgut_bubble <- melt(hindgut_bubble, id = c("Body_Site"))

test3 <- randomColor(39)

plot3 <- ggplot(hindgut_bubble, aes(x = Body_Site, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.0001, .2), range = c(5,50), breaks = c(0.0001,0.001, 0.01, 0.1)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")+ 
  scale_x_discrete(limits = c("Cecum", "Spiral_Colon", "Distal_Colon","Fecal"),
                   labels = c("Cecum", "Spiral_Colon", "Distal_Colon","Fecal")) +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +  
  scale_fill_manual(values = test3, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(hindgut_bubble$variable))) 
plot3
