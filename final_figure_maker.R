
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')


library("phyloseq")
library(ggplot2)
#library(tidyr)
library(reshape2)
library(dplyr)
library(vegan)
library(gridExtra)
library(factoextra)
library(ggsci)
library(car)
library(plyr)

paltest2 <- c("#fabebe","#d2f53c",
              "#f032e6","#46f0f0","#911eb4","#ffe119",
              "#3cb44b","#f58231","#0082c8","#e6194b",
              "#808080","#fc2c00","#ffd8b1","#808000",
              "#aaffc3","#800000","#000080","#aa6e28",
              "#e6beff","#008080")


speciestable <-read.table("/home/jenny/Desktop/wls_resistance/metagenome2019/humann/Rinput/fullmetaphlanRELABnocutoffMINcutoff_s.tsv",header=TRUE)
#need to add something to mask R1/2 duplication
speciestableR1 <- filter(speciestable, R1==TRUE)
speciestableRESEQ <-filter(speciestableR1,multipleseq=='YES')
speciestableOnlyMiceR1 <-filter(speciestableR1,Host=='Mouse')
speciestableOnlyMiceR1time <-filter(speciestableOnlyMiceR1,daysto10!='NA')
speciestableOnlyMiceR1_39 <-filter(speciestableOnlyMiceR1,Donor=='WLS39')
#make a matrix version from long data
wide_data <- dcast(speciestable,cage+multipleseq+set+cagediscordant+daysto10+R1+Host+Donor+SampleID~Taxa,value.var="Percent")
wide_dataR1 <- dcast(speciestableR1,cage+multipleseq+set+cagediscordant+daysto10+R1+Host+Donor+SampleID~Taxa,value.var="Percent")
wide_dataRESEQ <-dcast(speciestableRESEQ,cage+multipleseq+set+cagediscordant+daysto10+R1+Host+Donor+SampleID~Taxa,value.var="Percent")
wide_dataOnlyMiceR1 <-dcast(speciestableOnlyMiceR1,cage+multipleseq+set+cagediscordant+daysto10+R1+Host+Donor+SampleID~Taxa,value.var="Percent")
wide_dataOnlyMiceR1time <-dcast(speciestableOnlyMiceR1time,cage+multipleseq+set+cagediscordant+daysto10+R1+Host+Donor+SampleID~Taxa,value.var="Percent")
wide_dataOnlyMiceR1_39 <-dcast(speciestableOnlyMiceR1_39,cage+multipleseq+set+cagediscordant+daysto10+R1+Host+Donor+SampleID~Taxa,value.var="Percent")



onlyabundanceALL <- wide_dataR1[10:ncol(wide_data)]
onlyabundanceR1 <- wide_dataR1[10:ncol(wide_dataR1)]
onlyabundanceRESEQ <-wide_dataRESEQ[10:ncol(wide_dataRESEQ)]
onlyabundanceOnlyMiceR1 <-wide_dataOnlyMiceR1[10:ncol(wide_dataRESEQ)]
onlyabundanceOnlyMiceR1time <-wide_dataOnlyMiceR1time[10:ncol(wide_dataRESEQ)]
onlyabundanceOnlyMiceR1_39 <-wide_dataOnlyMiceR1_39[10:ncol(wide_dataRESEQ)]


#phyloseq object?
rownames(onlyabundanceRESEQ) <- wide_dataRESEQ$SampleID
otumat<-t(onlyabundanceRESEQ)
OTU <- otu_table(otumat,taxa_are_rows=TRUE)
TAXtable <-read.table("/home/jenny/Desktop/wls_resistance/metagenome2019/humann/Rinput/taxa_table.txt",header=TRUE,row.names=1)
TAX <-tax_table(as.matrix(TAXtable))
meta_info <- wide_dataRESEQ[1:9]
rownames(meta_info) <-wide_dataRESEQ$SampleID
sampledat <- sample_data(meta_info)

physeq <- phyloseq(OTU,TAX,sampledat)
plot_bar(physeq,fill="Family")
plot_heatmap(physeq)
testord<-ordinate(physeq,"PCoA","bray")
p1 <- plot_ordination(physeq,testord,type="sample",color="Donor")
p1

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="SampleType")
}, physeq, dist)

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
plot(x=testord$vectors['Axis.1'],y=testord$vectors['Axis.2'])


##Calculate Diversity
shannon<-diversity(onlyabundanceR1, index = "shannon", MARGIN = 1, base = exp(1))
#diversity(df2, index = "simpson", MARGIN = 1, base = exp(1))
invsimpson<-diversity(onlyabundanceR1, index = "invsimpson", MARGIN = 1, base = exp(1))
#link to matrix with all metadata, 
wide_data_with_diversity <- cbind(shannon,wide_dataR1)
wide_data_with_diversity <-cbind(invsimpson,wide_data_with_diversity)


###Looking for differences based on set (batch effects)
#wide_data_with_diversityR1 <- filter(wide_data_with_diversityR1,Host=='Mouse')
hist(wide_data_with_diversity$shannon)
p1<-ggplot(wide_data_with_diversity, aes(x=set,y=shannon)) + geom_boxplot()
p2<-ggplot(wide_data_with_diversity, aes(x=set,y=invsimpson)) + geom_boxplot()
grid.arrange(p1,p2,nrow=2)

#this might be wrong, should be bray-curtis?
myPCA<- prcomp(onlyabundanceR1, center=TRUE, scale.=TRUE)
df_out <- as.data.frame(myPCA$x)
ggplot(df_out,aes(x=PC1,y=PC2,color=wide_dataR1$set)) + scale_shape_manual(values=c(15, 16, 17, 18)) +  geom_point(aes(shape=wide_dataR1$Host),size=3) +
  labs(x="PC1 (14.97%)", y="PC2 (10.57%)", color="Batch",shape="Host") + coord_fixed(1)
summary(myPCA)
screeplot(myPCA,main="PCA Variances",bstick=TRUE)
fviz_screeplot(myPCA)
#barplot of ro/column contributions, dashed line corresponds to expected value if contribution were uniform(anything above considered"important")
fviz_contrib(myPCA, choice="var")

##making nmds plot
myNMDS <- metaMDS(onlyabundanceR1,distance='bray',k=2,trymax=1000)
#Look for stress below 0.2 according to Dill-McFarland -- here too much stress
#use k=3 b/c stress >0.2, but haven't got a 3d display to work yet. Maybe doesn't matter?
myNMDS <- metaMDS(onlyabundanceR1,distance='bray',k=3,trymax=1000)
stressplot(myNMDS)

#Batch plot, R1 only
plot(myNMDS, type="n")
points(myNMDS, display="sites", pch=20, col=pal_startrek()(2)[wide_dataR1$set])
legend(0.9,-0.8, title="Batch",legend=c("first","second"), col=pal_startrek()(2), pch=20)
ordiellipse(myNMDS,groups=wide_dataR1$set,col=pal_startrek()(2),conf=0.95)

#RESEQ NMDS
myNMDS <- metaMDS(onlyabundanceRESEQ,distance='bray',k=2,trymax=1000)
plot(myNMDS, type="n")
points(myNMDS, display="sites", pch=20, col=pal_startrek()(2)[wide_dataRESEQ$set])
text(myNMDS,labels=wide_dataRESEQ$Donor,cex=0.5)
ordiellipse(myNMDS,groups=wide_dataRESEQ$set,col=pal_startrek()(2),conf=0.95)
ordiellipse(myNMDS,groups=wide_dataRESEQ$Donor,col=pal_tron()(5),conf=0.95)
legend(0.5,-0.8, title="Batch",legend=c("first","second"), col=pal_startrek()(2), pch=20)

BCxyz <-scores(myNMDS,display="sites")
#plot(BCxyz[,1],BCxyz[,2], main="Bray-Curtis",pch=20, col=pal_startrek()(2)[wide_dataR1$set],xlab="NMDS1",ylab="NMDS2")
#ordiellipse(myNMDS,groups=wide_dataR1$set,col=pal_startrek()(2),conf=0.95)
#plot(BCxyz[,1],BCxyz[,3], main="Bray-Curtis",pch=20, col=pal_startrek()(2)[wide_dataR1$set],xlab="NMDS1",ylab="NMDS3")
#legend(-0.5,-0.5, title="Batch",legend=c("first","second"), col=pal_startrek()(2), pch=20,cex=.5)

#display in ggplot #HOST
data.scores <-as.data.frame(scores(myNMDS))
data.scores$site <- rownames(data.scores)
ggplot() +geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=wide_dataR1$Host)) + scale_shape_manual(values=c(15, 21))+
   scale_color_manual(values=pal_startrek()(2)) + labs(shape="Host") +coord_fixed(1) +guides(fill=FALSE)
#+confidenceEllipse(data=data.scores,aes(x=NMDS1,y=NMDS2,color=wide_dataR1$Host))
#+stat_ellipse(data=data.scores,aes(x=NMDS1,y=NMDS2,color=wide_dataR1$Host)) #+theme_classic() 

#DONOR X HOST from R1
ggplot() +geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=wide_dataR1$Donor,shape=wide_dataR1$Host,size=3)) + scale_shape_manual(values=c(15, 16, 17, 18))+
  scale_color_manual(values=paltest2) + labs(col="Donor",shape="Host",size=NA) +coord_fixed(1) 

ggplot(data=data.scores) +geom_point(aes(x=NMDS1,y=NMDS2,color=wide_dataR1$Donor,bg=wide_dataR1$Donor,shape=wide_dataR1$Host)) +
  scale_shape_manual(values=c(15, 21))+ scale_fill_manual(values=paltest2) + 
  scale_color_manual(values=paltest2)+ labs(col="Donor",shape="Host",size=NA) +coord_fixed(1) +
  guides(fill=FALSE) #+stat_ellipse(aes(x=NMDS1,y=NMDS2,color=wide_dataR1$Donor))

##only mouse R1
myNMDS <- metaMDS(onlyabundanceOnlyMiceR1,distance='bray',k=2,trymax=1000)
data.scores <-as.data.frame(scores(myNMDS))
data.scores$site <- rownames(data.scores)
ggplot() +geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=wide_dataOnlyMiceR1$Donor)) + 
  scale_color_manual(values=paltest2) + labs(col="Donor") +coord_fixed(1)


##only mouse R1, color by timeto10
myNMDS <- metaMDS(onlyabundanceOnlyMiceR1time,distance='bray',k=2,trymax=1000)
data.scores <-as.data.frame(scores(myNMDS))
data.scores$site <- rownames(data.scores)
ggplot() +geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=wide_dataOnlyMiceR1time$Donor)) + 
  scale_color_manual(values=paltest2) + labs(col="Donor") +coord_fixed(1)

fit.onlymicer1time <- envfit(myNMDS,wide_dataOnlyMiceR1time) #doesn't do what I want yet
#fit.onlymicer1time <- envfit(myNMDS,wide_dataOnlyMiceR1time[,c("Host")])
ggplot() +geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=wide_dataOnlyMiceR1time$daysto10)) + 
    scale_color_gradient(low="red",high="yellow") + coord_fixed(1) + labs(col="Days to Death") #scale_color_manual(values=paltest2)


##only mice, R1, and 39
myNMDS <- metaMDS(onlyabundanceOnlyMiceR1_39,distance='bray',k=2,trymax=1000)
data.scores <-as.data.frame(scores(myNMDS))
data.scores$site <- rownames(data.scores)
ggplot() +geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,color=wide_dataOnlyMiceR1_39$SampleID)) + 
  coord_fixed(1) + labs(col="Days to Death") +geom_jitter() #scale_color_manual(values=paltest2)


ordiplot(myNMDS,type="n")
ordihull(myNMDS,groups=wide_dataR1$set,draw="polygon",col=c("blue","green"),label=T)
orditorp(myNMDS,display="species",col="red",air=0.0)
orditorp(myNMDS,display="sites",air=0.0001,cex=1.45)

testdist <- vegdist(onlyabundanceRESEQ,distance="bray")
dispersionHost <- betadisper(testdist,wide_dataRESEQ$set)
permutest(dispersionHost,permutations=1000)
adonis(testdist ~ wide_dataRESEQ$set,data=wide_dataRESEQ,permutations=1000)
adonis(testdist ~ wide_dataRESEQ$Donor,data=wide_dataRESEQ,permutations=1000)

testdist <- vegdist(onlyabundanceR1,distance="bray")
dispersionHost <- betadisper(testdist,wide_dataR1$set)
permutest(dispersionHost,permutations=1000)
adonis(testdist ~ wide_dataR1$set,data=wide_dataR1,permutations=1000)
adonis(testdist ~ wide_dataR1$Donor,data=wide_dataR1,permutations=1000)

testdist <- vegdist(onlyabundanceOnlyMiceR1,distance="bray")
dispersionHost <- betadisper(testdist,wide_dataOnlyMiceR1$set)
permutest(dispersionHost,permutations=1000)
adonis(testdist ~ wide_dataOnlyMiceR1$set,data=wide_dataOnlyMiceR1,permutations=1000)
adonis(testdist ~ wide_dataOnlyMiceR1$Donor,data=wide_dataOnlyMiceR1,permutations=1000)

testdist <- vegdist(onlyabundanceOnlyMiceR1_39,distance="bray")
dispersionHost <- betadisper(testdist,wide_dataR1$set)
permutest(dispersionHost,permutations=1000)
adonis(testdist ~ wide_dataOnlyMiceR1_39$set,data=wide_dataOnlyMiceR1_39,permutations=1000)
adonis(testdist ~ wide_dataOnlyMiceR1_39$dayto10,data=wide_dataOnlyMiceR1_39,permutations=1000)
adonis(testdist ~ wide_dataR1$Donor,data=wide_dataR1,permutations=1000)


#test which taxa are associated with seqeuncing batch? RESEQ
source("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/simper_pretty.R")
source("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/R_krusk.R")
meta=wide_dataRESEQ
simper.pretty(onlyabundanceRESEQ,meta,c("set"),perc_cutoff=1,low_cutoff='y',low_val=0.01,'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/RESEQset')
simper.results = data.frame(read.csv("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/RESEQset_clean_simper.csv"))
kruskal.pretty(onlyabundanceRESEQ,meta,simper.results,c("set"),'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/reseqset')
KW.results<-data.frame(read.csv('/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/reseqset_krusk_simper.csv'))
KW.results.signif = KW.results[KW.results$fdr_krusk_p.val < 0.05,]
KW.results.signif = KW.results.signif[with(KW.results.signif,order(OTU)),]
head(KW.results.signif)

##which are associated just with batch??
meta=wide_dataR1
simper.pretty(onlyabundanceR1,meta,c("set"),perc_cutoff=1,low_cutoff='y',low_val=0.01,'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/R1set')
simper.results = data.frame(read.csv("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/R1set_clean_simper.csv"))
kruskal.pretty(onlyabundanceR1,meta,simper.results,c("set"),'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/r1set')
KW.results<-data.frame(read.csv('/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/r1set_krusk_simper.csv'))
KW.results.signif = KW.results[KW.results$fdr_krusk_p.val < 0.05,]
KW.results.signif = KW.results.signif[with(KW.results.signif,order(OTU)),]
head(KW.results.signif)
write.table(KW.results.signif,file="/home/jenny/Desktop/wls_resistance/figures/batchsig.tsv",sep='\t')

#just39
meta=wide_dataR1
simper.pretty(onlyabundanceR1,meta,c("set"),perc_cutoff=1,low_cutoff='y',low_val=0.01,'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/R1set')
simper.results = data.frame(read.csv("/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/R1set_clean_simper.csv"))
kruskal.pretty(onlyabundanceR1,meta,simper.results,c("set"),'/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/r1set')
KW.results<-data.frame(read.csv('/home/jenny/Desktop/wls_resistance/figures/seq_scripts-masterSteinberger/r1set_krusk_simper.csv'))
KW.results.signif = KW.results[KW.results$fdr_krusk_p.val < 0.05,]
KW.results.signif = KW.results.signif[with(KW.results.signif,order(OTU)),]
head(KW.results.signif)
write.table(KW.results.signif,file="/home/jenny/Desktop/wls_resistance/figures/batchsig.tsv",sep='\t')

