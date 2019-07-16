library(ggplot2)
#library(tidyr)
library(reshape2)
library(dplyr)
library(vegan)
library(gridExtra)


##Input data and manipulate format
#Read in data, change if data stored somewhere else
phylatest <- read.table("/home/jenny/Desktop/wls_resistance/bothmetagenome/relabundforR_p.tsv",header=TRUE)
speciestable <-read.table("/home/jenny/Desktop/wls_resistance/metagenome2019/humann/Rinput/fullmetaphlanRELAB_s.tsv",header=TRUE)

#need to add something to mask R1/2 duplication, Humans vs mouse
#can also add filters later where appropriate
phylatestR1 <- filter(phylatest, R1==TRUE)

#make a matrix version from long data
wide_data <- dcast(phylatest,R1+Host+Donor+SampleID~Taxa,value.var="Percent")
head(wide_data)
wide_dataR1 <- dcast(phylatestR1,daysto10+R1+Host+Donor+SampleID~Taxa,value.var="Percent")

#remove non-numeric names or metadata
onlyabundanceALL <- wide_data[5:ncol(wide_data)]
onlyabundance <- wide_dataR1[6:ncol(wide_data)]
wide_dataR1_10 <- filter(wide_dataR1,daysto10!='NA')
onlyabundance10 <- wide_dataR1_10[6:ncol(wide_data)]


##Calculate Diversity
shannon<-diversity(onlyabundance, index = "shannon", MARGIN = 1, base = exp(1))
#diversity(df2, index = "simpson", MARGIN = 1, base = exp(1))
invsimpson<-diversity(onlyabundance, index = "invsimpson", MARGIN = 1, base = exp(1))

#link to matrix with all metadata, 
wide_data_with_diversity <- cbind(shannon,wide_data)
wide_data_with_diversity <-cbind(invsimpson,wide_data_with_diversity)
#write.csv(shannon,"~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/shannon.csv")
#write.csv(invsimpson,"~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/invsimpson.csv")

#calculate with error bars for each group, statitically test differences
wide_data_with_diversityR1 <- filter(wide_data_with_diversity,R1==TRUE)
wide_data_with_diversityR1 <- filter(wide_data_with_diversityR1,Host=='Mouse')
hist(wide_data_with_diversityR1$shannon)
#TO DO: get shorter names or angle names, scatter points instead of box, print HSD above, n=? above
#TO DO: use sp or genus level here?
#Plot in one figure
p1<-ggplot(wide_data_with_diversityR1, aes(x=Donor,y=shannon)) + geom_boxplot()
p2<-ggplot(wide_data_with_diversityR1, aes(x=Donor,y=invsimpson)) + geom_boxplot()
grid.arrange(p1,p2,nrow=2)

#Test to see if diversity is statiscally different, will only say if there is a difference
#https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
anova_diversity <- aov(shannon~Donor,wide_data_with_diversityR1)
summary(anova_diversity)
#Tukey's Honest Signifcant Difference (HSD) for pairwise comparsions of means
library(agricolae)
tukey_diversity<-HSD.test(anova_diversity,"Donor",group=TRUE)
print(tukey_diversity)

##Richness plot + calculate if sampling 
#can calculate in mothur/QIIME instead, might be more informative doing it then that at the end
#plot(specaccum(onlyabundance,groups=wide_data$Donor)) #doesn't work
#TO DO: only accepts integers (count data)
rarecurve(onlyabundance)



##Stacked Bar Charts
#plot for each sample
ggplot(data=phylatest,aes(x=SampleID,y=Percent,fill=Taxa)) +geom_bar(stat="identity")+
  scale_fill_manual(values=paltest2) + labs(x="sample",y="Relative Abundance")

#Split-apply-combine based on group of interest
testgrouping <- group_by(phylatest, Donor,Taxa) %>% summarize(mean(Percent))
testgrouping <- testgrouping %>% rename('Percent'='mean(Percent)')

#plot based on split groups
ggplot(data=testgrouping,aes(x=Donor,y=Percent,fill=Taxa)) +geom_bar(stat="identity")+
  scale_fill_manual(values=paltest2) + labs(x="sample",y="Relative Abundance")


##making PCA plot
#get coordinates
myPCA<- prcomp(onlyabundance, center=TRUE, scale.=TRUE)
df_out <- as.data.frame(myPCA$x)
ggplot(df_out,aes(x=PC1,y=PC2,color=wide_data$Donor)) + scale_shape_manual(values=c(15, 16, 17, 18)) +  geom_point(aes(shape=wide_data$Host),size=3) +
    labs(x="PC1 (39.3%)", y="PC2 (23.5%)", color="Donor",shape="Host") + coord_fixed(1)
summary(myPCA)
#where broken stick is above variance, should use that componenet, where below not so important
#Broken stick has been recommended as a stopping rule in principal component analysis (Jackson 1993): principal components should be retained as long as observed eigenvalues are higher than corresponding random broken stick components.
screeplot(myPCA,main="PCA Variances",bstick=TRUE)
library(factoextra)
fviz_screeplot(myPCA)
#barplot of ro/column contributions, dashed line corresponds to expected value if contribution were uniform(anything above considered"important")
fviz_contrib(myPCA, choice="var")
#get proportion of variance with summary(myPCA2) to change x labels
#recommend to adjust aspect ratio based on Nguyen et al 2019 10 tips DR
#add scree plot to be more concious about number of variables chosen

##making nmds plot
myNMDS <- metaMDS(onlyabundance,distance='bray',k=2,trymax=1000)
#Look for stress below 0.2 according to Dill-McFarland
stressplot(myNMDS)
#large scatter around the line suggests original dissimilarities are not well preserved

ordiplot(myNMDS,type="n")
ordihull(myNMDS,groups=wide_dataR1$Host,draw="polygon",col=c("blue","green"),label=T)
orditorp(myNMDS,display="species",col="red",air=0.01)
orditorp(myNMDS,display="sites",air=0.01)

#plot with lines showing continuous data
#need to remove NAs
myNMDS10 <- metaMDS(onlyabundance10,distance='bray',k=2,trymax=1000)
plot(myNMDS10,type="n",main="Bray-Curtis")
points(myNMDS10,pch=20,display="sites",col=c("blue","red")[wide_dataR1$Host])
fitTime10 <- envfit(myNMDS10,wide_dataR1_10$daysto10)
plot(fitTime10)
#note that can plot only significant arrows if wanted

#Madison statistically test beta diversity
#PERMANOVA --does not assume normality but does assume equal beta dispersion
#if beta dispersion signifcantly different, use ANOSIM

testdist <- vegdist(onlyabundance,distance="bray")
dispersionHost <- betadisper(testdist,wide_dataR1$Host)
permutest(dispersionHost,permutations=1000)

adonis(testdist ~ wide_data$R1,data=wide_data,permutations=1000)
adonis(testdist ~ wide_data$Donor,data=wide_data,permutations=1000)
adonis(testdist ~ wide_data$Host,data=wide_data,permutations=1000)

#ANOSIM only allows 1 simple variable model, no interactions, and only categorical

##Taxa over-underabundance
#SIMPER identifies OTUs that most contribute to beta-diversity (most abundandant and/or most variable)
#Only for categorical variables

simper(onlyabundance,wide_dataR1$Donor,permutations=100)
kruskal.test(onlyabundance$Other ~ meta$Donor)
#Andrew Steinberger SIMPER script to perform SIMPER
#Then takes output of above to run kruskal-wallis on each OTU
#can this be replaced with apply?
source("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/simper_pretty.R")
source("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/R_krusk.R")
meta = wide_dataR1[0:5]
meta=wide_dataR1
simper.pretty(onlyabundance,meta,c("Donor"),perc_cutoff=1,low_cutoff='y',low_val=0.01,'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/Donor')
simper.results = data.frame(read.csv("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/Donor_clean_simper.csv"))
kruskal.pretty(onlyabundance,meta,simper.results,c("Donor"),'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/Donork')
KW.results<-data.frame(read.csv('/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/Donork_krusk_simper.csv'))
KW.results.signif = KW.results[KW.results$fdr_krusk_p.val < 0.05,]
KW.results.signif = KW.results.signif[with(KW.results.signif,order(OTU)),]
head(KW.results.signif)
#make boxplot for signicant results
#make phyloseq object for other comparisons
