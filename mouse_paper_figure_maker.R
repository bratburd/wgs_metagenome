#mouse paper figure maker

library(ggplot2)
library(gplots)
library(scales)
library(vegan)
library(gridExtra)


####
#antismash humsal day 0 vs day 3
#arranged in 3 columns: day, cluster type, percentage
#bgccol <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/test-mg-bgc-box-mapsepsalcolumn-smallvalues-larger.csv")
bgccol <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/anti-11-10-17copro_mouse_r_rawdivtot-onlyday0day3.tbl")
#copro_mouse_r_percentage_out2 <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/copro_mouse_r_percentage_out2.tbl")

#bbduk m3
bgccol <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/bbduk/anti-m3-4-5-18copro_mouse_r_rawdivtotnorm.tbl")


#bbduk m1
bgccol <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/bbduk/anti-m1-4-5-18copro_mouse_r_rawdivtotnorm-noday1-no0.tbl")
bgccol$time <- as.factor(bgccol$time)
#levels(birds$effect) <- gsub(" ", "\n", levels(birds$effect))
levels(bgccol$genome) <-gsub("-", "\n", levels(bgccol$genome))

ptest <- ggplot(bgccol, aes(x=0, y = raw, fill=time)) +
  geom_boxplot() +
  scale_y_sqrt(breaks=c(0,1,2,3))+
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) + 
  facet_grid(. ~ reorder(genome,-raw,mean))  +  
  #theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.text.x=element_text(size=8))+
  labs(y="Percentage of total reads", fill="Day") + scale_x_continuous(breaks = NULL) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(position=position_jitterdodge(dodge.width=1),aes(fill=time),color="black",pch=21)
ptest
#statisical tests for each group
lantipeptide0 <- bgccol$genome == 'lantipeptide' & bgccol$time == '0'
lantipeptide3 <- bgccol$genome == 'lantipeptide' & bgccol$time == '3'
wtest <- wilcox.test(bgccol[lantipeptide0,]$raw,bgccol[lantipeptide3,]$raw)
pvalue1 <- wtest$p.value

cluster0 <- bgccol$genome == 'putative' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'putative' & bgccol$time == '3'
wtest <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)
pvalue2 <- wtest$p.value

cluster0 <- bgccol$genome == 'arylpolyene' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'arylpolyene' & bgccol$time == '3'
wtest <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)
pvalue3 <- wtest$p.value

cluster0 <- bgccol$genome == 'fatty_acid' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'fatty_acid' & bgccol$time == '3'
pvalue4 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value

cluster0 <- bgccol$genome == 'NRPS' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'NRPS' & bgccol$time == '3'
pvalue5 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value

cluster0 <- bgccol$genome == 'terpene' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'terpene' & bgccol$time == '3'
pvalue6 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value #ties

cluster0 <- bgccol$genome == 'sactipeptide' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'sactipeptide' & bgccol$time == '3'
pvalue7 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value

#cluster0 <- bgccol$genome == 'resorcinol' & bgccol$time == 'day 0'
#cluster3 <- bgccol$genome == 'resorcinol' & bgccol$time == 'day 3'
#pvalue8 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value

cluster0 <- bgccol$genome == 'saccharide' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'saccharide' & bgccol$time == '3'
pvalue9 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value

cluster0 <- bgccol$genome == 'fatty_acid saccharide' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'fatty_acid saccharide' & bgccol$time == '3'
pvalue10 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value

cluster0 <- bgccol$genome == 'thiopeptide' & bgccol$time == '0'
cluster3 <- bgccol$genome == 'thiopeptide' & bgccol$time == '3'
pvalue11 <- wilcox.test(bgccol[cluster0,]$raw,bgccol[cluster3,]$raw)$p.value

all_p <- c(pvalue1, pvalue2, pvalue3,pvalue4,pvalue5,pvalue6,pvalue7,pvalue9,pvalue10,pvalue11)
pvaladj <- p.adjust(all_p, method="fdr")

#####
#relative abundance
copro_mouse_r_percentage_out2 <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/copro_mouse_r_percentage_out2.tbl")
copro_mouse_r_percentage_out2 <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/mchbplasmid_11-15-17nomouse_forR.tbl")
copro_mouse_r_percentage_out2 <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/mchb_plasmid_11-15-17min4copro_mouse_r_percent-nomousenorm.tbl")
copro_mouse_r_percentage_outbb <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/bbduk-mg-4-2-17copro_mouse_r_percent-nomousenorm.tbl")



#m1
copro_mouse_r_percentage_outbb <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/bbduk/mg-bbduk-m1-4-6-18copro_mouse_r_percent-nomousenorm.tbl")


counts_gather <- copro_mouse_r_percentage_outbb
counts_gather$genome <- factor(counts_gather$genome, levels=counts_gather$genome)
counts_gather$genome <- factor(counts_gather$genome, levels=rev(counts_gather$genome))
paltest2 <- c("#808080","#fc2c00","#ffd8b1","#808000",
              "#aaffc3","#800000","#000080","#aa6e28",
              "#e6beff","#008080","#fabebe","#d2f53c",
              "#f032e6","#46f0f0","#911eb4","#ffe119",
              "#3cb44b","#f58231","#0082c8","#e6194b")
              #"#d01c8b","#fb8072","#e41a1c","#377eb8") 
ggplot (data=counts_gather, aes(x=number,y=raw ,fill=genome)) + geom_bar(stat="identity") +
  scale_fill_manual(values=paltest2) + facet_grid(pathogen~time,labeller=label_both) +
  labs(x="Mouse", y="Relative Abundance",fill="Taxa") + scale_x_continuous(breaks = unique(counts_gather$number))

#display palette on pie chart
pie(rep(1, length(paltest2)), labels = sprintf("%d (%s)", seq_along(paltest2), 
                                          paltest2), col = paltest2)

altpal <- c("#d01c8b","#fb8072","#e41a1c","#377eb8", "#d01c8b","#fb8072","#e41a1c","#377eb8")
pie(rep(1, length(altpal)), labels = sprintf("%d (%s)", seq_along(altpal), 
                                               altpal), col = altpal)

#####
#pca plot
#need plot with genomes across top, samples down 1 column, extra data in first three columns
m <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/normraw_for_pca_withOUTpathogen2.csv",row.names=1)
m <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/bbduk/normraw-nomouse-nosal-nopli.csv",row.names=1)

#m <- rawreads_mouse_microbiome[2:45]
percentabund <- t(m)/colSums(m)

write.csv(percentabund,"~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/percentabund_forpcaWITHPATH.csv")

write.csv(percentabund,"~/GRAD SCHOOL/currie/mouse-paper-fig/bbduk/percentabund_forpcaWITHOUTPATH.csv")
#with pathogens
reads_mouse_microbiome <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/percentabund_forpcaWITHPATH.csv")

#without pathogens
#reads_mouse_microbiome <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/pca_nopli_nomousecandsal.csv")
reads_mouse_microbiome <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/percentabund_forpcaWITHOUTPATH.csv")
reads_mouse_microbiome <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/bbduk/percentabund_forpcaWITHOUTPATH.csv")


treatment <- reads_mouse_microbiome$treatment
d <-reads_mouse_microbiome$day_of_infection
day_of_infection <- as.factor(d)
abs_community <- reads_mouse_microbiome[4:93] #or 95 with pathogens, 93 without
myPCA3<- prcomp(abs_community, center=TRUE, scale.=TRUE)
df_out <- as.data.frame(myPCA3$x)
ggplot(df_out,aes(x=PC1,y=PC2,color=treatment)) + scale_shape_manual(values=c(15, 18, 17, 16)) +
  geom_point(aes(shape=day_of_infection),size=3) + labs(x="PC1 (31.5%)", y="PC2 (18.1%)", color="Treatment",shape="Day of Infection")

summary(myPCA3)
#get proportion of variance with summary(myPCA2) to change x labels


reads_mouse_microbiome <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/percentabund_forpcaWITHPATH.csv")
community <- reads_mouse_microbiome[4:93]
df <- subset(reads_mouse_microbiome, select = -c(treatment, day_of_infection))

df2 <- df[,-1]

rownames(df2) <- df[,1]

shannon<-diversity(df2, index = "shannon", MARGIN = 1, base = exp(1))
#diversity(df2, index = "simpson", MARGIN = 1, base = exp(1))
invsimpson<-diversity(df2, index = "invsimpson", MARGIN = 1, base = exp(1))

write.csv(shannon,"~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/shannon.csv")
write.csv(invsimpson,"~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/invsimpson.csv")

alldiversity <- read.csv("~/GRAD SCHOOL/currie/mouse-paper-fig/assorted_data/diversity-withPATH.csv")

newdata <- mydata[ which(gender=='F' & age > 65),]
day0and3 <- alldiversity[which(alldiversity$day_of_infection == 0 | alldiversity$day_of_infection == 3),]
day0and3 <- day0and3[which(day0and3$treatment == 'Salmonella' | day0and3$treatment=='Candida'),]
                               alldiversity$treatment == 'Salmonella' | alldiversity$treatment=='Candida'),]
                               #&alldiversity$treatment=='Salmonella'),]
#day0sal <- day0and3[which(day0and3$treatment == 'Salmonella' & day0and3$day_of_infection == '0'),]
day0and3cand <- alldiversity[which(alldiversity$day_of_infection == 0 | alldiversity$day_of_infection == 3 &alldiversity$treatment=='Candida'),]
day0and3sal <- alldiversity[which(alldiversity$day_of_infection == 0 | alldiversity$day_of_infection == 3 &alldiversity$treatment=='Salmonella'),]


day0and3$day_of_infection <- as.factor(day0and3$day_of_infection)
bp1 <- ggplot(day0and3, aes(x=day_of_infection, y=shannon, fill=day_of_infection)) + geom_boxplot() + facet_grid(. ~ treatment)+
  labs(x="Days Post Infection",y="Shannon Diversity")+ theme(legend.position="none")
#bp2 <- ggplot(day0and3cand, aes(x=day_of_infection, y=shannon)) + geom_boxplot()
bp3 <- ggplot(day0and3, aes(x=day_of_infection, y=invsimpson, fill=day_of_infection)) + geom_boxplot() + facet_grid(. ~ treatment)+
  labs(x="Days Post Infection",y="Inverse Simpson Index") + theme(legend.position="none")
#bp4 <- ggplot(day0and3cand, aes(x=day_of_infection, y=invsimpson)) + geom_boxplot()
grid.arrange(bp1, bp3, nrow = 2)
kruskal.test(shannon ~ day_of_infection, data=day0and3cand)
kruskal.test(shannon ~ day_of_infection, data=day0and3sal)

kruskal.test(invsimpson ~ day_of_infection, data=day0and3cand)
kruskal.test(invsimpson ~ day_of_infection, data=day0and3sal)

View(reads_mouse_microbiome)
newdata <- reads_mouse_microbiome[c(0:34,37:39,41,43:44),]
treat2 <- newdata$treatment
d2 <-newdata$day_of_infection
day_of_infection2 <- as.factor(d2)
abscomm<- newdata[4:93]
myPCAnoout<-prcomp(abscomm,center=TRUE)
df_out2 <- as.data.frame(myPCAnoout$x)
ggplot(df_out2)
ggplot(df_out2,aes(x=PC1,y=PC2,color=treat2)) + scale_shape_manual(values=c(15, 18, 17, 16)) +
  geom_point(aes(shape=day_of_infection2),size=3) + labs(x="PC1 (50.7%)", y="PC2 (16.2%)", color="Treatment",shape="Day of Infection")
summary(myPCAnoout)
