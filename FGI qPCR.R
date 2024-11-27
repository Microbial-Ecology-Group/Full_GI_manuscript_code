# Library -----------------------------------------------------------------
library(btools)
library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(metagenomeSeq)
library(metagMisc)
library(pairwiseAdonis)
library(phyloseq)
library(picante)
library(randomcoloR)
library(scales)
library(stringr)
library(vegan)
library(devtools)
library(ggdendro)
library(multcompView)
library(rcompanion)


#set wd
setwd("C:/Users/danie/OneDrive - West Texas A and M University/Research/Full GI/qPCR/")

# Read in data ------------------------------------------------------------
dat <- read.csv("Raw data.csv")

summary(dat$mean_copy)

# Analysis ----------------------------------------------------------------

#shapirowilk test - Make sure your excel file has your means listed as numbers
shapiro.test(dat$mean_copy) #p= 4.017e-13 not normal

#log transform
dat$mean_log <- log(dat$mean_copy)

#sharpirowilk
shapiro.test(dat$mean_log)

#kruskal test for location
kruskal.test(dat$mean_copy,dat$body_site) #p< 2.2e-16 highly significant

#pairwise test
qpcr_anova <- pairwise.wilcox.test(dat$mean_copy, dat$body_site, p.adjust.method = "BH")
#assign superscripts
qpcr_pvalues <- qpcr_anova$p.value %>% 
  fullPTable()
multcompLetters(qpcr_pvalues)

qpcr_anova_log <- pairwise.wilcox.test(dat$mean_log, dat$body_site, p.adjust.method = "BH")
#assign superscripts
qpcr_log_pvalues <- qpcr_anova_log$p.value %>% 
  fullPTable()
multcompLetters(qpcr_log_pvalues)

#checking plate and LAYoN, and harvest group
kruskal.test(dat$mean_copy,dat$LAYoN)#ns
kruskal.test(dat$mean_copy,dat$plate)#ns
kruskal.test(dat$mean_copy,dat$harvest_group)#ns


# Plots -------------------------------------------------------------------
#palate
bodysite_palette <- c("red4", "royalblue1", "navy", "gold", "springgreen", "goldenrod4", 'goldenrod', "darkviolet", "darkorange1", "orangered2", "royalblue3")


ggplot(dat, aes(x= body_site, y= mean_copy, fill = body_site)) +
  theme_bw() + labs(y= "", x="", title = "") +
  scale_x_discrete(limits=c("Oral", "Rumen", "Abomassum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
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

#calculate Averages
dat2 <- dat%>% group_by(body_site) %>% 
  summarize(mean(mean_copy, na.rm=F),
            sd(mean_copy, na.rm = F))

#plot harvest group
dat$harvest_group <- factor(dat$harvest_group)
ggplot(dat, aes(x= harvest_group, y= mean_copy, fill = harvest_group)) +
  theme_bw() + labs(y= "", x="", title = "") +
  #scale_x_discrete(limits=c("Oral", "Rumen", "Abomassum", "Duodenum", "Jejunum", "Ileum", "Cecum", "Spiral_Colon", "Distal_Colon","Fecal"))+
  geom_boxplot() +
  geom_point()+
  scale_y_log10()+
  scale_fill_manual(values = c("grey0", "brown"))+
  theme(#legend.position = "none",
        plot.margin = unit(c(0.75,0.75,0.75,0.75),"cm"),
        panel.border = element_rect(colour = "black", size = 1.7),
        axis.ticks = element_line(size = 1, colour = "black"),
        plot.title = element_text(size = 30),
        axis.title.y = element_text(size = 30, vjust = 2.5),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_text(size = 24, colour = "black"),
        panel.grid.major.x = element_blank())

