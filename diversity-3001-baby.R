rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename<- args[1]
filename="GSS3001.metaphlan.estimated_count.filtered.7"



#-- could not do beta analysis as some empty can I add 1 to all 

test <- read.table(filename, header = TRUE, row.names = 1, sep="\t")
test[is.na(test)]<-0
d=dim(test)
d
X=apply(test, 1, sum)
test_filter2=test[(X/d[2] >=10),]

test_d=data.frame(t(test_filter2))
dim(test_d)
min(apply(test_d, 1, sum))


library(vegan)
# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
test_a <- decostand(test_d, method = "total")
# check total abundance in each sample
#apply(test_a, 1, sum) ### this is just for confirmation
#------------------------------------------
write.table(t(test_a), file =paste0(filename, "recal.Percent.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

Cname=rownames(test_d)

Clen=length(Cname) ##there are 3 annotation columns

splitname<-strsplit(Cname, "_")
Collect=rep("NA", Clen)
Sid=rep("NA", Clen)
BodySite=rep("NA", Clen)
id2=rep("NA", Clen)
Visit=rep("NA", Clen)
for(mm in  1:Clen ){
  Collect[mm]=splitname[[mm]][5]
  Sid[mm]=splitname[[mm]][2]
  BodySite[mm]=splitname[[mm]][4]
  Visit[mm]=splitname[[mm]][3]
  id2[mm]=splitname[[mm]][6]
}

NewID=paste0(Collect, ".", BodySite, ".", Sid)

SiteSID=paste0(BodySite, ".", Sid)

SiteCollect=paste0(BodySite, Collect)
SiteCollectSid=paste0(BodySite, Collect, Sid)
SidCollect=paste0(Sid, Collect)
library(stringr)

library(ggplot2)

#https://stackoverflow.com/questions/30057765/histogram-ggplot-show-count-label-for-each-bin-for-each-category

shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(test_d, 2) - 1
## Fisher alpha
#fisher <- fisher.alpha(test_d)
####observed species
observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)
#unbias.simp
mydata1=data.frame(shannon,simp,invsimp,  observed, Menhinick_index, Margalef_index, BodySite, SiteSID, Collect, id2, Sid, SiteCollect, SiteCollectSid, SidCollect)
mydata=mydata1[Collect!="Control" & Sid != "SkinMC" & BodySite != "Blank",]                   

write.table(t(mydata1), file =paste0(filename, "alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

saveRDS(mydata1,paste0(filename, ".diversitydata"))
saveRDS(test_d,paste0(filename,".test_d"))



library(ggplot2)

#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
library(tidyverse)
library(ggpubr)
library(rstatix)

png(filename=paste0(filename,"shannon.Collect.ttest.png"))
stat.test <- mydata1 %>%
  #group_by(Country) %>%
  t_test(shannon ~ Collect) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 3.2)
bxp <- ggboxplot(mydata1, x = "Collect", y = "shannon", color = "Collect")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

png(filename=paste0(filename,"shannon.Collect.wilcox.png"),width=280, height = 480)
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ Collect) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 3.2)
bxp <- ggboxplot(mydata, x = "Collect", y = "shannon", color = "Collect")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

png(filename=paste0(filename,"shannon.BodySite.wilcox.png"), width=480, height = 480)
stat.test <- mydata1 %>%
  wilcox_test(shannon ~ BodySite) %>%
  mutate(y.position = 3.2)
bxp <- ggboxplot(mydata1, x = "BodySite", y = "shannon", color = "BodySite")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()

png(filename=paste0(filename,"shannon.Sid.wilcox.png"), width=480, height = 480)
stat.test <- mydata %>%
  wilcox_test(shannon ~ Sid) %>%
  mutate(y.position = 3.2)
bxp <- ggboxplot(mydata, x = "Sid", y = "shannon", color = "Sid")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()


png(filename=paste0(filename,"shannon.Collect.BodySite.ttest.png"), height =480,width=480)
bxp <- ggboxplot(mydata1, x = "Collect", y = "shannon", color = "Collect", facet.by = "BodySite")
bxp + stat_compare_means(comparisons = list(c("Swab1", "Tape1")), method = "t.test")
dev.off()

png(filename=paste0(filename,"shannon.Collect.BodySite.wilcox.png"), height = 600, width=600)
bxp <- ggboxplot(mydata1, x = "Collect", y = "shannon", color = "Collect", facet.by = "BodySite")
bxp + stat_compare_means(comparisons = list(c("Swab1", "Tape1")), method = "wilcox")

dev.off()

wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                        mydata$Collect, 
                                        p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- Collect wilcoxon adjustp")
tab.shannon


wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$BodySite, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- BodySite wilcoxon adjustp")
tab.shannon


wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$SiteCollect, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- BodySite & Collect wilcoxon adjustp")
tab.shannon

wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$Sid, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- PanelistID wilcoxon adjustp")
tab.shannon


wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$SidCollect, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- PanelistID & Collect wilcoxon adjustp")
tab.shannon



png(filename=paste0(filename,"observed.bodysite.ttest.png"), width=480, height=280)
stat.test <- mydata1 %>%
  t_test(observed ~ BodySite) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 75)
bxp <- ggboxplot(mydata1, x = "BodySite", y = "observed", color = "BodySite")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()


png(filename=paste0(filename,"observed.BodySite.wilcox.png"), width=480, height=280)
stat.test <- mydata1 %>%
  wilcox_test(observed ~ BodySite) %>%
  mutate(y.position = 75)
bxp <- ggboxplot(mydata1, x = "BodySite", y = "observed", color = "BodySite")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()

png(filename=paste0(filename,"observed.Collect.wilcox.png"), width=280, height=280)
stat.test <- mydata %>%
  wilcox_test(observed ~ Collect) %>%
  mutate(y.position = 150)
bxp <- ggboxplot(mydata, x = "Collect", y = "observed", color = "Collect")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()

png(filename=paste0(filename,"observed.Collect.Site.wilcox.png"), height = 600, width=600)
bxp <- ggboxplot(mydata1, x = "Collect", y = "observed", color = "Collect", facet.by = "BodySite")
bxp + stat_compare_means(comparisons = list(c("Swab1", "Tape1")), method = "wilcox")
dev.off()


########add in to show statistical table directly
wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
                                        mydata$Collect, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("observed -- Collect wilcoxon adjustp")
tab.observed

wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
                                        mydata$BodySite, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("observed -- Site wilcoxon adjustp")
tab.observed

wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
                                        mydata$SiteCollect, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("observed -- Site Swabber wilcoxon adjustp")
tab.observed

KP_Collect=kruskal.test(shannon ~ Collect, data = mydata)$p.value
KP_BodySite=kruskal.test(shannon ~BodySite, data = mydata)$p.value
KP_SiteCollect=kruskal.test(shannon ~ SiteCollect, data = mydata)$p.value
KP_Sid=kruskal.test(shannon ~Sid, data = mydata)$p.value
KP_SidCollect=kruskal.test(shannon ~SidCollect, data = mydata)$p.value
print(paste0("KP_Collect=", KP_Collect, "; KP_BodySite=", KP_BodySite,"; KP_SiteCollect=", KP_SiteCollect,"; KP_SidCollect=", KP_SidCollect, "; KP_Sid=", KP_Sid
     ))
#[1] "KP_Swabber=0.470663419146733; KP_BodySite=1.87746075779544e-13; KP_SiteSwab=7.7498086336419e-10; KP_Visit=0.625837099162698; KP_Sid=0.0284246903301365"

x<-test_d+1 ###to avoid 0 
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- cbind(as.data.frame(mds$points), mydata1)
mds_data$SampleID <- rownames(mds_data)

library(stringr)

library(ggplot2)
 
 png(filename=paste0(filename,".beta.Collect.png"),  width=380, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Collect)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Collect, fill=Collect), alpha=0.1, geom="polygon")
 print(p1)
 dev.off()
 
 png(filename=paste0(filename,".beta.Site.png"),  width=380, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = BodySite)) +geom_point(size=6)+theme_bw()+
   stat_chull(aes(color=BodySite, fill=BodySite), alpha=0.1, geom="polygon")
 print(p1)
 dev.off()
 
 png(filename=paste0(filename,".beta.Site.Collect.png"),  width=480, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = BodySite)) +geom_point(size=3)+theme_bw()+
   stat_chull(aes(color=BodySite, fill=BodySite), alpha=0.05, geom="polygon")+ facet_grid(. ~ Collect)
 print(p1)
 dev.off()
 
 png(filename=paste0(filename,".beta.Collect.Site.png"),  width=800, height=280)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Collect)) +geom_point(size=3)+theme_bw()+stat_chull(aes(color=Collect, fill=Collect), alpha=0.1, geom="polygon")+ facet_grid(. ~ BodySite)
 print(p1)
 dev.off()
######################################
bc.dist <- vegdist(x, method = "bray")
jaccard.dist <- vegdist(x, method = "jaccard")
mds2 <- metaMDS(jaccard.dist)
mds_data2 <- as.data.frame(cbind(mds2$points, mydata1))
mds_data2$SampleID <- rownames(mds_data2)


png(filename=paste0(filename,".beta.Collect.BodySite.jaccard.png"),  width=800, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = Collect)) +geom_point(size=6)+theme_bw()+
  stat_chull(aes(color=Collect, fill=Collect), alpha=0.1, geom="polygon")+ facet_grid(. ~ BodySite)
print(p1)
dev.off()

# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram
png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=800, height=660)
plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
dev.off()

print("adonis result bc.dist ~ Collect")
adonis2(bc.dist ~ Collect)
#p=0.015 *

test.adonis <- adonis(bc.dist ~ BodySite)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis
#p=0.001
###################################################################
#PAIRWISE PERMANOVA
cbn <- combn(x=unique(SiteCollect), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(SiteCollect %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(SiteCollect %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$SiteCollect)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
# Output the results to a file
write.table(p.table, "beta_site_swab.txt", sep = "\t", quote = FALSE)


cbn <- combn(x=unique(BodySite), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(BodySite %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(BodySite %in% cbn[,i]), ]
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$BodySite)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
write.table(p.table, "beta_bodysite.txt", sep = "\t", quote = FALSE)


cbn <- combn(x=unique(Collect), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(Collect %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(Collect %in% cbn[,i]), ]
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$Collect)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
write.table(p.table, ".beta_Collect.txt", sep = "\t", quote = FALSE)
###########################Jaccard test
print("adonis result jaccard.dist ~ BodySite")
adonis2(jaccard.dist ~ BodySite)
test.adonis <- adonis(jaccard.dist ~ BodySite)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

test.adonis <- adonis(jaccard.dist ~ Collect)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

#