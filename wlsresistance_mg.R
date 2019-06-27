#plots from metagenomic data
library(ggplot2)
library(gplots)
library(scales)
library(vegan)
library(gridExtra)


copro_mouse_r_percentage_outbb <- read.delim("~/GRAD SCHOOL/currie/mouse-paper-fig/bbduk/mg-bbduk-m1-4-6-18copro_mouse_r_percent-nomousenorm.tbl")
wlsall <-read.delim("file:///C:/Users/Jbabe/Documents/metagenomes/combine201819relab.tsv")
taxa <- wlsall$X.SampleID
library(tidyr)
library(stringr)
library(broman)
library(dplyr)
test %>% separate(taxa, c("domain", "phylum", "class", "order", "family", "genus"), ",[a-z]|")
wlsall %>%
  stringr::str_extract("[k,p,c,o,f,g,s,t](?=__[\\w]+$)") %>%
  switchv(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", 
          f = "Family", g = "Genus", s = "Species", t = "Strain")

mutate(wlsall)

library(reshape2)
only28 <-read.csv("file:///C:/Users/Jbabe/Documents/metagenomes/28familyfilter.csv")
long <- melt(only28, id.vars = c("Taxa"),value.name="raw")
names(long)[names(long) == 'variable'] <- 'number'
#long <- melt(only28, varnames = c("Taxa","number"),value.name="raw")
View(long)

only39 <-read.csv("file:///C:/Users/Jbabe/Documents/metagenomes/39familyfilter.csv")
long <- melt(only39, id.vars = c("Taxa"),value.name="raw")
names(long)[names(long) == 'variable'] <- 'number'
#long <- melt(only28, varnames = c("Taxa","number"),value.name="raw")
View(long)
levels(long$Taxa)

donorphyla<-read.csv("file:///C:/Users/Jbabe/Documents/metagenomes/donorsphyla.csv")
long <- melt(donorphyla, id.vars = c("Taxa"),value.name="raw")
long$Taxa <- factor(long$Taxa, levels = rev(bruteorder))
names(long)[names(long) == 'variable'] <- 'number'
View(long)
long$Taxa<-reorder(long$Taxa, rowSums(long$raw))
bruteorder <- c("Firmicutes","Bacteroidetes","Verrucomicrobia","Actinobacteria","Proteobacteria","Other")
long$Taxa<-reorder(long$Taxa, bruteorder)
## To reorder the levels:
## note, if x is not a factor use levels(factor(x))
x = factor(x,levels(x)[c(4,5,1:3)])

long <- factor(long,levels(long$Taxa)[bruteorder])

phylatest <- read.table("file:///C:/Users/Jbabe/Documents/metagenomes/relabundforR_p.tsv", header=TRUE)
ggplot(data=phylatest,aes(x=SampleID,y=Percent,fill=Taxa)) +geom_bar(stat="identity")+
  scale_fill_manual(values=paltest2) + labs(x="sample",y="Relative Abundance")

mean(phylatest$Donor)

counts_gather <- long
counts_gather$Taxa <- factor(counts_gather$Taxa, levels=counts_gather$Taxa)
counts_gather$Taxa <- factor(counts_gather$Taxa, levels=rev(counts_gather$Taxa))
paltest2 <- c("#808080","#fc2c00","#ffd8b1","#808000",
              "#aaffc3","#800000","#000080","#aa6e28",
              "#e6beff","#008080","#fabebe","#d2f53c",
              "#f032e6","#46f0f0","#911eb4","#ffe119",
              "#3cb44b","#f58231","#0082c8","#e6194b")
#"#d01c8b","#fb8072","#e41a1c","#377eb8") 
ggplot (data=counts_gather, aes(x=number,y=raw ,fill=Taxa)) + geom_bar(stat="identity") +
  scale_fill_manual(values=paltest2) + #facet_grid(pathogen~time,labeller=label_both) +
  labs(x="Mouse", y="Relative Abundance",fill="Taxa") + scale_x_continuous(breaks = unique(counts_gather$number))

counts_gather %>% mutate(Taxa = factor(Taxa, levels = Taxa))
testorder<-factor(long$Taxa)

ggplot(data=counts_gather,aes(x=number,y=raw,fill=Taxa)) +geom_bar(stat="identity")+
  scale_fill_manual(values=paltest2) + labs(x="sample",y="Relative Abundance")


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data=counts_gather,aes(x=number,y=raw,fill=Taxa)) +geom_bar(stat="identity")+
   scale_fill_manual(values=cbPalette) +labs(x="Donor",y="Relative Abundance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

counts_gather %>%
  arrange(desc(Taxa)) %>%
  ggplot(aes(x=number, y=raw)) +
  geom_area(aes(fill= Taxa), position = 'stack') +
  #scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_fill_brewer(palette = "YlOrBr")



#display palette on pie chart
pie(rep(1, length(paltest2)), labels = sprintf("%d (%s)", seq_along(paltest2), 
                                               paltest2), col = paltest2)

altpal <- c("#d01c8b","#fb8072","#e41a1c","#377eb8", "#d01c8b","#fb8072","#e41a1c","#377eb8")
pie(rep(1, length(altpal)), labels = sprintf("%d (%s)", seq_along(altpal), 
                                             altpal), col = altpal)
