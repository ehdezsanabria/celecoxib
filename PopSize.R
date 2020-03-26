#### load libraries ####

library(readxl)
library(tidyverse)
library(phyloseq)
library(RColorBrewer)
library(vegan)
library(grid)
library(gridExtra)
library(egg)
library(Phenoflow)
source("functions.R")
source("normalization.R")
source("highqualgraphr.R")
#### load data ####

dataset <- readxl::read_excel("Cell_count_Diversity.xlsx",sheet = "ForR",col_types = c(rep("numeric",4),"text"))
dilcords <- read.csv2("DilutionCorrected.csv")
levels(dilcords$Evenness.tag2)[levels(dilcords$Evenness.tag2)=="H"] <- "High"
levels(dilcords$Evenness.tag2)[levels(dilcords$Evenness.tag2)=="M"] <- "Medium"
levels(dilcords$Evenness.tag2)[levels(dilcords$Evenness.tag2)=="L"] <- "Low"
dilcords$Evenness.tag2 <- factor(dilcords$Evenness.tag2 ,levels(dilcords$Evenness.tag2 )[c(1,3,2)])

dilcords$timepoint <- as.factor(dilcords$timepoint)
levels(dilcords$timepoint)[levels(dilcords$timepoint)=="1"] <- "Transfer 1"
levels(dilcords$timepoint)[levels(dilcords$timepoint)=="2"] <- "Transfer 2"
levels(dilcords$timepoint)[levels(dilcords$timepoint)=="3"] <- "Transfer 3"
levels(dilcords$timepoint)[levels(dilcords$timepoint)=="4"] <- "Transfer 4"
levels(dilcords$timepoint)[levels(dilcords$timepoint)=="5"] <- "Transfer 5"

#### plotting ####

p1 <- dilcords %>%  ggplot(aes(x=as.factor(timepoint),y=CellsCount,fill=Evenness.tag2,colour=Evenness.tag2))
p1 + geom_boxplot() + facet_wrap(~Evenness.tag2) + ylab("Total cell count (events/mL)") + xlab("timepoint")


p2 <- dilcords %>%  ggplot(aes(x=Inverse,y=CellsCount,
                              fill=Evenness.tag2,pch=Evenness.tag2))
# p2 + geom_point() + facet_wrap(~timepoint,scales = "free_y") + 
#      ylab("Total cell count (events/mL)") + xlab("D2") + 
#      scale_color_discrete(name="Initial Evenness", breaks=c("H","M","L")) +
#      scale_shape_discrete(name="Initial Evenness", breaks=c("H","M","L"))
# p2 + geom_point() + facet_wrap(~timepoint) + 
#      ylab("Total cell count (events/mL)") + xlab("D2") + 
#      scale_color_discrete(name="Initial Evenness", breaks=c("H","M","L")) +
#      scale_shape_discrete(name="Initial Evenness", breaks=c("H","M","L"))
# p2 + geom_point() + facet_grid(~timepoint) + 
#      ylab("Total cell count (events/mL)") + xlab("D2") + 
#      scale_color_discrete(name="Initial Evenness", breaks=c("H","M","L")) +
#      scale_shape_discrete(name="Initial Evenness", breaks=c("H","M","L")) + theme_bw()

p2.save <- p2 + geom_point(shape = 21, size =2.75, alpha = 0.6) + facet_grid(~timepoint) + 
  ylab("Total cell density (Cells/mL)") + xlab("Phenotypic diversity (a.u.)") + theme_bw() +
  labs(y = bquote("Total cell density (cells mL"^{-1}*")"), fill = "Initial evenness")+
  scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62"))+
  labs(y = bquote("Total cell density (cells mL"^{-1}*")"))+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        axis.text = element_text(size=14),title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2))
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )

pdf(file = "Fig2_v2.pdf", width = 12, height = 6)
p2.save
dev.off()

png(file = "Fig2_v2.png", width = 12, height = 6, res = 500, units = "in")
p2.save
dev.off()

# highqualgraphR(p2.save,filename = "TotalvsD2grid",extension = c("pdf","png"))


#### Comparison with Locey et al 2016 (PNAS) ####

# p2 + geom_point() + facet_wrap(~timepoint,scales = "free") + 
#   ylab("log10 Total cell count (events/mL)") + xlab(expression(log[10]~D[2])) + 
#   scale_x_log10() + scale_y_log10() + annotation_logticks() +
#   scale_color_discrete(name="Initial Evenness", breaks=c("H","M","L")) +
#   scale_shape_discrete(name="Initial Evenness", breaks=c("H","M","L")) + theme_bw() 

#### Generalized regression for Figure 1####
library(mgcv)
library(nlme)
library(lmtest)
library(car)
library(multcomp)
library(ade4)
ctrl <- lmeControl(opt='optim')

# Import raw analysed data including biological replicates label
data.lme <- read.csv2("results_raw.csv")
# Make sure to have more than 1,000 cells per sample
data.lme <- data.lme[data.lme$Total>1000, ]
data.lme <- data.lme[data.lme$Evenness.tag2 != "B",]
# Format this data a bit
data.lme <- data.lme[, c(1,3,12,13,14,15)]
data.lme <- data.lme[data.lme$timepoint<6, ]; data.lme <- droplevels(data.lme)
data.lme$Sample <- gsub(data.lme$Sample,pattern=" .fcs",replacement=".fcs")
data.lme$rep <- as.factor(data.lme$rep); data.lme$Evenness.tag <- as.factor(data.lme$Evenness.tag)

levels(data.lme$Evenness.tag2)[levels(data.lme$Evenness.tag2)=="H"] <- "High"
levels(data.lme$Evenness.tag2)[levels(data.lme$Evenness.tag2)=="M"] <- "Medium"
levels(data.lme$Evenness.tag2)[levels(data.lme$Evenness.tag2)=="L"] <- "Low"
data.lme$Evenness.tag2 <- factor(data.lme$Evenness.tag2 ,levels(data.lme$Evenness.tag2 )[c(1,3,2)])
data.lme$timepoint <- as.factor(data.lme$timepoint)
levels(data.lme$timepoint)[levels(data.lme$timepoint)=="1"] <- "Transfer 1"
levels(data.lme$timepoint)[levels(data.lme$timepoint)=="2"] <- "Transfer 2"
levels(data.lme$timepoint)[levels(data.lme$timepoint)=="3"] <- "Transfer 3"
levels(data.lme$timepoint)[levels(data.lme$timepoint)=="4"] <- "Transfer 4"
levels(data.lme$timepoint)[levels(data.lme$timepoint)=="5"] <- "Transfer 5"

# Run model
data.lme$TimeEv <- interaction(data.lme$timepoint, data.lme$Evenness.tag2)
lme.el <- lme(Inverse ~ TimeEv,random=~1|rep,
              correlation = corAR1(form=~1|rep/Evenness.tag), 
              data=data.lme,
              method="REML", control = ctrl)


lm.resid <- lm(residuals(lme.el,type="normalized")~predict(lme.el))

# Test for initial interaction effect
anova(lme.el)

# Run posthoc analysis to test for difference between timepoints/transfers
posthoc_test <- summary(glht(lme.el, linfct = mcp(TimeEv = "Tukey")))

# plot model diagnostics
png("FigSI_1.elham.png",width=15,height=5,units="in",pointsize=12,res=500)
par(mfrow=c(1,3),mar=c(7,7,7,7))
qqPlot(residuals(lme.el,type="normalized"),col="blue",pch=21,ylab="",cex=1.5,xlab="",las=1)
mtext(side=2,"Studentized residuals",cex=2,line=3.75)
mtext(side=1,"Normal distribution quantiles",cex=2,line=3.25)
plot(residuals(lme.el,type="normalized")~predict(lme.el),ylab="",las=1,cex.axis=1.5,xlab="")
mtext(side=2,"Studentized residuals",cex=2,line=3.5)
mtext(side=1,"Fitted values",cex=2,line=3.75)
lines(x=predict(lme.el),y=predict(lm.resid),col="red")
plot(predict(lme.el),data.lme$Inverse,las=1,ylab="",cex.axis=1.5,xlab="")
mtext(side=2,"Diversity",cex=2,line=3.75)
mtext(side=1,"Fitted values",cex=2,line=3.25)
lines(x=data.lme$Inverse,y=data.lme$Inverse,col="red")
lines(col="red",x=predict(lme.el),y=predict(lm.resid))
dev.off()

# Update dataframe
data.lme <- data.frame(data.lme, model_fit = predict(lme.el))

# Plot data + regression
p3 <- data.lme %>%  ggplot(aes(y = Inverse, x = timepoint))+
  geom_boxplot(aes(y = Inverse, x = timepoint,
                   alpha = TimeEv), notch = FALSE, outlier.shape = NA)+
  geom_point(aes(fill = Evenness.tag2),shape = 21, size = 3, alpha = 0.4, 
             position = position_jitterdodge(jitter.width = 0.2)) +
  # geom_point(aes(x = timepoint, y = model_fit, color = Evenness.tag2), shape = "-", size = 20)+
  ylab("Phenotypic diversity (a.u.)") + xlab("") + theme_bw() +
  labs(fill = "Initial evenness") +
  scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62"))+
  scale_color_manual(values=c("#88419d","#a6cee3","#fc8d62"))+
  theme(axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_text(size=15), legend.text=element_text(size=14),
        axis.text = element_text(size=14), title=element_text(size=20),
        strip.background=element_rect(fill=adjustcolor("lightgray",0.2)),
        #,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_blank()
  )+
  geom_crossbar(aes(x = timepoint, y = model_fit, fill = Evenness.tag2), ymin=0, ymax=0, 
                position = position_jitterdodge(),width = 0.5, size = 1.05)+ 
  geom_crossbar(aes(x = timepoint, y = model_fit, color = Evenness.tag2), ymin=0, ymax=0, 
                position = position_jitterdodge(), width = 0.4)+ 
  guides(color = FALSE, alpha = FALSE)+
  facet_grid(~timepoint, scales = "free_x")
  

  
png("Fig1_B.png",width=12,height = 6,units="in",pointsize=10,res=500)
p3
dev.off()

pdf(file = "Fig1_B.pdf", width = 12, height = 6)
p3
dev.off()

#### Create 16s --- Figure 3 ####
# Import rarefied data
data_16S <- read.csv("Indices_ALL_evenness_paper.csv")

# Format for Phyloseq
meta_16S <- data_16S[,c(1:3)]; rownames(meta_16S) <- data_16S$Sample
levels(meta_16S$tpt)[levels(meta_16S$tpt)=="t0"] <- "Transfer 0"
levels(meta_16S$tpt)[levels(meta_16S$tpt)=="t5"] <- "Transfer 5"
meta_16S$Evenness.tag <- as.numeric(gsub(meta_16S$Sample, pattern = ".*n", 
                                         replacement = ""))

otu_16S <- data_16S[,c(4:ncol(data_16S))]; rownames(otu_16S) <- data_16S$Sample
TAX_16S <- matrix(colnames(otu_16S)); rownames(TAX_16S) <- colnames(otu_16S)

# Make phyloseq object
physeq <- phyloseq(otu_table(otu_16S, taxa_are_rows = FALSE), sample_data(meta_16S), tax_table(as.matrix(TAX_16S)))

# Make PCoA
pcoa <- ordinate(
  physeq = physeq, 
  method = "PCoA", 
  distance = "bray",
  correction = "lingoes",
  k=2
)
pcoa.df <- data.frame(pcoa$vectors[,1:2], sample_data(physeq))

# Make beta-diversity plot
myColours2 <- brewer.pal(n=12,"Paired"); myColours2 <- myColours2[c(2,6,7,9)]

# Create distance matrix for PERMANOVA
dist.seq <- vegan::vegdist(otu_table(physeq))

# Check homogeneity of variances
betadisper(dist.seq, group=pcoa.df$tpt)

# Run PERMANOVA
permanova <- adonis(dist.seq ~ tpt * evenness, data = pcoa.df)


my_grob = grobTree(textGrob(bquote(paste(r[Transfer]^2 == 
                                           .(round(100 * permanova$aov.tab[1, 5], 1)), 
                                         "%")), x = 0.72, y = 0.95, hjust = 0, gp = gpar(col = "black", 
                                                                                        fontsize = 12, fontface = "italic")))
my_grob2 = grobTree(textGrob(bquote(paste(r[Evenness]^2 == 
                                            .(format(round(100 * permanova$aov.tab[2, 5], 
                                                           1), nsmall = 1)), "%")), x = 0.72, y = 0.87, 
                             hjust = 0, gp = gpar(col = "black", fontsize = 12, 
                                                  fontface = "italic")))
my_grob3 = grobTree(textGrob(bquote(paste(r[Transfer:Evenness]^2 == 
                                            .(round(100 * permanova$aov.tab[3, 5], 1)), 
                                          "%")), x = 0.72, y = 0.79, hjust = 0, gp = gpar(col = "black", 
                                                                                         fontsize = 12, fontface = "italic")))

var <- round(pcoa$values$Eigenvalues/sum(pcoa$values$Eigenvalues)*100,1)

beta.pcoa <- ggplot(data=pcoa.df, aes(x=Axis.1, y=Axis.2, fill=evenness))+
  geom_point(alpha=0.7, size=7, aes(shape=tpt))+
  scale_shape_manual(values=c(21,24))+
  theme_bw()+
  scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62"))+
  scale_colour_manual(values=myColours2)+
  labs(x = paste0("PCoA axis 1 (",var[1], "%)"), y = paste0("PCoA axis 2 (",var[2], "%)"), fill="Initial Evenness", shape="Transfer", colour="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  guides(fill = guide_legend(override.aes = list(shape = 22)))+
  annotation_custom(my_grob)+
  annotation_custom(my_grob2)+
  annotation_custom(my_grob3)+
  labs(title="(A)")


png("Fig3.png",width=9, height = 5,units="in", pointsize=10,res=500)
beta.pcoa
dev.off()

pdf(file = "Fig3.pdf", width = 9, height = 5)
beta.pcoa
dev.off()

# Make stacked barplots for t0 and t5
# Convert abundances to percentages
physeq_rel <- transform_sample_counts(physeq, function(x) 100*x/sum(x))

# Convert phyloseq to dataframe
df_physeq_rel <- psmelt(physeq_rel)

# Calculate mean and sd relative abundance for each timepoint and initial evenness
# class
df_physeq_rel_subset <- aggregate(Abundance~tpt/evenness/OTU, df_physeq_rel, mean)
df_physeq_rel_subset_sd <- aggregate(Abundance~tpt/evenness/OTU, df_physeq_rel, sd)

df_physeq_rel_subset <- data.frame(df_physeq_rel_subset, sd_Abundance = df_physeq_rel_subset_sd[, 4])
df_physeq_rel_subset$evenness <-  factor(df_physeq_rel_subset$evenness, levels = c("Low", "Medium", "High"))

p5 <- ggplot(aes(x = tpt, y = Abundance, fill = OTU), data = df_physeq_rel_subset)+
  geom_bar(stat='identity')+
  facet_grid(~evenness)+
  ylab("Relative Abundance (%)")+ xlab("Transfer")+
  scale_fill_manual(values=brewer.pal(n=12,"Paired"))+
  theme_bw()+  
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
                     title=element_text(size=20), legend.text=element_text(size=16))

# Only compare t5
df_physeq_rel_subset_t5 <- df_physeq_rel_subset[df_physeq_rel_subset$tpt=="T5",]

p6 <- ggplot(aes(x = evenness, y = Abundance, fill = OTU), data = df_physeq_rel_subset_t5)+
  geom_bar(stat='identity')+
  ylab("Relative Abundance (%)")+ xlab("Initial Evenness")+
  scale_fill_manual(values=brewer.pal(n=12,"Paired"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16))+
  labs(title="(B)")

png("Fig4.png", width=9, height = 10,units="in", pointsize=10, res=500)
grid.arrange(beta.pcoa, p6, nrow=2)
dev.off()

pdf(file = "Fig4.pdf", width = 9, height = 10)
grid.arrange(beta.pcoa,p6, nrow=2)
dev.off()

#### Compare FCM/actual diversity data with observed 16S  ####
div_16S <- Diversity_16S(physeq, R = 999, brea = FALSE,
                         parallel = TRUE, ncores = 2)
div_16S <- data.frame(sample = rownames(div_16S), div_16S)
div_16S <- left_join(div_16S, meta_16S, by = c("sample" = "Sample"))
div_16S$TimeEv <- interaction(div_16S$tpt, div_16S$Evenness.tag)

### Compare 16S diversity with phenotypic diversity at t5 ###

# for the fcm data we first have to average over measured replicates
dilcords$TimeEv <- interaction(dilcords$timepoint,dilcords$Evenness.tag)
dilcordsDav <- data.frame(D2.fcm = aggregate(Inverse~TimeEv, dilcords, mean),
                          D1.fcm = aggregate(exp(Shannon)~TimeEv, dilcords, mean)[,2],
                          D0.fcm = aggregate(Observed_R1~TimeEv, dilcords, mean)[,2]
)
dilcordsDsd <- data.frame(D2.fcm.sd = aggregate(Inverse~TimeEv, dilcords, sd)[,2],
                          D1.fcm.sd = aggregate(exp(Shannon)~TimeEv, dilcords, sd)[,2],
                          D0.fcm.sd = aggregate(Observed_R1~TimeEv, dilcords, sd)[,2]
)

div_FCM <- cbind(dilcordsDav, dilcordsDsd)
colnames(div_FCM)[1:2] <- c("TimeEv","D2.fcm")

# Add evenness label (Evenness.tag2)
div_FCM <- inner_join(div_FCM, dilcords, by = "TimeEv")
# div_FCM <- div_FCM[, c(1:3, 17)]

# extract joint samples at t5
df_16S_fcm <- inner_join(div_16S, div_FCM, by = "TimeEv")

# Only retain data of interest
df_16S_fcm <- df_16S_fcm[, c(1:3,8:21)]

# Adjust level order
df_16S_fcm$evenness <- factor(df_16S_fcm$evenness, levels = c("High","Medium","Low"))


# Plot the data
p_ref_fcm_D2 <- ggplot(data=df_16S_fcm,aes(y=D2.fcm,x=D2, fill=evenness))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=evenness))+
  theme_bw()+labs(x=expression('Taxonomic diversity - D'[2]),y=expression('Phenotypic diversity - D'[2]), fill="Initial Evenness")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_x_continuous(trans='log2', breaks = seq(1,10, 1),minor_breaks =NULL) +
  scale_y_continuous(trans='log2', breaks = seq(1000,4250,250),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",formula=y~x)+
  geom_errorbar(aes(ymin = D2.fcm - D2.fcm.sd, ymax = D2.fcm + D2.fcm.sd))

p_ref_fcm_D1 <- ggplot(data=df_16S_fcm,aes(y=D1.fcm,x=D1, fill=evenness))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=evenness))+
  theme_bw()+labs(x=expression('Taxonomic diversity - D'[1]),y=expression('Phenotypic diversity - D'[1]), fill="Initial Evenness")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_x_continuous(trans='log2', breaks = seq(1,10, 1),minor_breaks =NULL) +
  scale_y_continuous(trans='log2', breaks = seq(1000,5000,500),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",formula=y~x)+
  geom_errorbar(aes(ymin = D1.fcm - D1.fcm.sd, ymax = D1.fcm + D1.fcm.sd))

p_ref_fcm_D0 <- ggplot(data=df_16S_fcm, aes(y=D0.fcm,x=D0, fill=evenness))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=evenness))+
  theme_bw()+labs(x=expression('Taxonomic diversity - D'[0]),y=expression('Phenotypic diversity - D'[0]), fill="Initial Evenness")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_x_continuous(trans='log2', breaks = seq(1,10, 1),minor_breaks =NULL) +
  scale_y_continuous(trans='log2', breaks = seq(1000,17000,1000),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",formula=y~x)+
  geom_errorbar(aes(ymin = D0.fcm - D0.fcm.sd, ymax = D0.fcm + D0.fcm.sd))

png("FigSI_X.png", width=9, height = 5,units="in", pointsize=10, res=500)
grid_arrange_shared_legend(p_ref_fcm_D0, p_ref_fcm_D2, ncol = 2)
dev.off()
png("FigSI_X.pdf", width = 9, height = 5)
grid_arrange_shared_legend(p_ref_fcm_D0, p_ref_fcm_D2, ncol = 2)
dev.off()

### Compare 16S diversity with actual diversity at t0 ###
div_ref <- read.csv2("Initial_evenness_t0.csv")
div_ref_16S <- inner_join(div_16S[div_16S$tpt == "Transfer 0",], div_ref, by = c("Evenness.tag" = "Sample"))
div_ref_16S <- div_ref_16S[div_ref_16S$tpt == "Transfer 0",]
# div_ref_16S <- div_ref_16S[, c(1,10,11,12,13,14,19)]
div_ref_16S$Evenness.tag <- as.numeric(div_ref_16S$Evenness.tag)
div_ref_16S$Evenness.tag2[div_ref_16S$Evenness.tag<41] <- "Low"
div_ref_16S$Evenness.tag2[div_ref_16S$Evenness.tag>40 & div_ref_16S$Evenness.tag<61] <- "Medium"
div_ref_16S$Evenness.tag2[div_ref_16S$Evenness.tag>60] <- "High"
div_ref_16S$Evenness.tag2 <- factor(div_ref_16S$Evenness.tag2, levels = c("High", "Medium", "Low"))

# Plot the data
p_ref_initE_D2 <- ggplot(data=div_ref_16S,aes(x=D2,y=InverseSimpson, fill=Evenness.tag2))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=Evenness.tag2))+
  theme_bw()+labs(y=expression('Initial diversity - D'[2]),x=expression('16S diversity - D'[2]), fill="Initial Evenness"
                  , title = "at Transfer 0")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(trans='log2', breaks = seq(1,10, 1),minor_breaks =NULL) +
  scale_x_continuous(trans='log2', breaks = seq(1,10, 1),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",formula=y~x)

p_ref_initE_D0 <- ggplot(data=div_ref_16S,aes(x=Evenness.tag2,y=D0, fill=Evenness.tag2))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_boxplot(alpha=0.6,aes(fill=Evenness.tag2), outlier.shape = NA)+
  geom_jitter(shape=21,size=6,alpha=0.6,aes(fill=Evenness.tag2), width = 0.3)+
  # geom_point(shape=21,size=6,alpha=0.6,aes(fill=Evenness.tag2))+
  theme_bw()+labs(y=expression('Taxonomic diversity - D'[0]), fill="Initial Evenness", x = "", title = "at Transfer 0")+
  theme(axis.text=element_text(size=15),axis.title.y =element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_y_continuous(breaks = seq(1,10, 1),minor_breaks =NULL) +
  geom_hline(yintercept = 10, linetype = 2, size = 2)

p_ref_fcm_D0_2 <- ggplot(data=df_16S_fcm, aes(y=D0.fcm,x=D0, fill=evenness))+ scale_fill_manual(values=c("#88419d","#a6cee3","#fc8d62")) +
  geom_point(shape=21,size=6,alpha=0.6,aes(fill=evenness))+
  theme_bw()+labs(x=expression('Taxonomic diversity - D'[0]),y=expression('Phenotypic diversity - D'[0]), fill="Initial Evenness",
                  title = "at Transfer 5")+
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=15),
        legend.title=element_text(size=16),strip.text.x = element_text(size = 22))+
  scale_x_continuous(trans='log2', breaks = seq(1,10, 1),minor_breaks =NULL) +
  scale_y_continuous(trans='log2', breaks = seq(1000,17000,1000),minor_breaks =NULL) +
  geom_smooth(method="lm",color="black",formula=y~x)+
  geom_errorbar(aes(ymin = D0.fcm - D0.fcm.sd, ymax = D0.fcm + D0.fcm.sd))


png("FigSI_2.png", width=9, height = 5,units="in", pointsize=10, res=500)
grid_arrange_shared_legend(p_ref_initE_D2, p_ref_fcm_D2, ncol = 2)
dev.off()

pdf(file = "FigSI_2.pdf", width = 9, height = 5)
grid_arrange_shared_legend(p_ref_initE_D2, p_ref_fcm_D2, ncol = 2)
dev.off()

png("FigSI_2_D0.png", width=9, height = 5,units="in", pointsize=10, res=500)
grid_arrange_shared_legend(p_ref_initE_D0, p_ref_fcm_D0_2, ncol = 2)
dev.off()

pdf(file = "FigSI_2_D0.pdf", width = 9, height = 5)
grid_arrange_shared_legend(p_ref_initE_D0, p_ref_fcm_D0_2, ncol = 2)
dev.off()
