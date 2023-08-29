rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename<- args[1]
filename="metaphlan.estimated_count.filtered.7"



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

Panelist=rep("NA", Clen)
Visit=rep("NA", Clen)
Swabber=rep("NA", Clen)
studyID=rep("NA", Clen)
BodySite=rep("NA", Clen)
flowcell=rep("NA", Clen)
scanOrder=rep("NA", Clen)

for(mm in  1:Clen ){
  Panelist[mm]=splitname[[mm]][3]
  Visit[mm]=splitname[[mm]][4]
  Swabber[mm]=splitname[[mm]][5]
  BodySite[mm]=splitname[[mm]][6]
  scanOrder[mm]=splitname[[mm]][7]
  flowcell[mm]=splitname[[mm]][2]
  studyID[mm]=splitname[[mm]][1]
}
SiteName=BodySite
SiteName[BodySite==1]="InnerThigh"
SiteName[BodySite==2]="MonsPubis"
SiteName[BodySite==3]="BikiniLegCrease"
SiteName[BodySite==4]="OutsideLabiaMajora"
SiteName[BodySite==5]="OutsideLabiaMinora"
SiteName[BodySite==6]="Perineum"
SiteName[BodySite==7]="Perianal"
SiteName[BodySite==8]="Introitus"
SiteName[BodySite==9]="Vaginal"
SiteName[Panelist=="Zymo"]="Zymo"
SiteName[is.na(Swabber)]="Other"
Sid=paste0(Panelist, "_", Visit)
SidSwab=as.factor(paste0(Sid, "_", Swabber))
SidSwabSite=as.factor(paste0(Sid, "_", Swabber, "_", SiteName))
SiteSwab=as.factor(paste0(SiteName, "_", Swabber))
VisitSwab=as.factor(paste0(Visit, "_", Swabber))
Swabber[is.na(Swabber)]="Other"


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
mydata1=data.frame(shannon,simp,invsimp,  observed, Menhinick_index, Margalef_index, SiteName,Swabber, BodySite,Panelist, Visit, scanOrder, flowcell, studyID, Sid, SidSwab, SidSwabSite, SiteSwab, VisitSwab)
mydata=mydata1[Swabber!="S73.metaphlan4",]                   

write.table(t(mydata), file =paste0(filename, "alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

saveRDS(mydata,paste0(filename, ".diversitydata"))
saveRDS(test_d,paste0(filename,".test_d"))



library(ggplot2)

#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
library(tidyverse)
library(ggpubr)
library(rstatix)

png(filename=paste0(filename,"shannon.swabber.ttest.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  t_test(shannon ~ Swabber) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 4.2)
bxp <- ggboxplot(mydata, x = "Swabber", y = "shannon", color = "Swabber")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

png(filename=paste0(filename,"shannon.Swabber.wilcox.png"))
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ Swabber) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 4.2)
bxp <- ggboxplot(mydata, x = "Swabber", y = "shannon", color = "Swabber")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
dev.off()

png(filename=paste0(filename,"shannon.BodySite.wilcox.png"), width=1200, height = 1200)
stat.test <- mydata %>%
  #group_by(Country) %>%
  wilcox_test(shannon ~ SiteName) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 4.2)
bxp <- ggboxplot(mydata, x = "SiteName", y = "shannon", color = "SiteName")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()

png(filename=paste0(filename,"shannon.Swabber.Site.ttest.png"), height = 1200, width=1200)
bxp <- ggboxplot(mydata, x = "Swabber", y = "shannon", color = "Swabber", facet.by = "SiteName")
bxp + stat_compare_means(comparisons = list(c("Nurse", "Self")), method = "t.test")
dev.off()

png(filename=paste0(filename,"shannon.Swabber.Site.wilcox.png"), height = 1200, width=1200)
bxp <- ggboxplot(mydata, x = "Swabber", y = "shannon", color = "Swabber", facet.by = "SiteName")
bxp + stat_compare_means(comparisons = list(c("Nurse", "Self")), method = "wilcox")
dev.off()

wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                        mydata$Swabber, 
                                        p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- Swabber wilcoxon adjustp")
tab.shannon


wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$SiteName, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- BodySite wilcoxon adjustp")
tab.shannon


wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$SiteSwab, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- BodySite & Swabber wilcoxon adjustp")
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
                                       mydata$SidSwab, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- PanelistID & Swabberwilcoxon adjustp")
tab.shannon

wilcox.shannon <- pairwise.wilcox.test(mydata$shannon, 
                                       mydata$Visit, 
                                       p.adjust.method = "BH")
tab.shannon <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("shannon -- Visit adjustp")
tab.shannon


png(filename=paste0(filename,"observed.SiteName.ttest.png"), width=1200, height=1200)
stat.test <- mydata %>%
  t_test(observed ~ SiteName) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 250)
bxp <- ggboxplot(mydata, x = "SiteName", y = "observed", color = "SiteName")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()


png(filename=paste0(filename,"observed.SiteName.wilcox.png"), width=1200, height=1200)
stat.test <- mydata %>%
  wilcox_test(observed ~ SiteName) %>%
  #adjust_pvalue() %>%
  mutate(y.position = 250)
bxp <- ggboxplot(mydata, x = "SiteName", y = "observed", color = "SiteName")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.08)
dev.off()


png(filename=paste0(filename,"observed.Swabber.Site.wilcox.png"), height = 1200, width=1200)
bxp <- ggboxplot(mydata, x = "Swabber", y = "observed", color = "Swabber", facet.by = "SiteName")
bxp + stat_compare_means(comparisons = list(c("Nurse", "Self")), method = "wilcox")
dev.off()


########add in to show statistical table directly
wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
                                        mydata$Swabber, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("observed -- Swabber wilcoxon adjustp")
tab.observed

wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
                                        mydata$SiteName, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("observed -- Site wilcoxon adjustp")
tab.observed

wilcox.observed <- pairwise.wilcox.test(mydata$observed, 
                                        mydata$SiteSwab, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.observed$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
print("observed -- Site Swabber wilcoxon adjustp")
tab.observed

KP_Swabber=kruskal.test(shannon ~ Swabber, data = mydata)$p.value
KP_SiteName=kruskal.test(shannon ~SiteName, data = mydata)$p.value
KP_SiteSwab=kruskal.test(shannon ~ SiteSwab, data = mydata)$p.value
KP_Sid=kruskal.test(shannon ~Sid, data = mydata)$p.value
KP_Visit=kruskal.test(shannon ~Visit, data = mydata)$p.value
print(paste0("KP_Swabber=", KP_Swabber, "; KP_BodySite=", KP_SiteName,"; KP_SiteSwab=", KP_SiteSwab,"; KP_Visit=", KP_Visit, "; KP_Sid=", KP_Sid
     ))
#[1] "KP_Swabber=0.470663419146733; KP_BodySite=1.87746075779544e-13; KP_SiteSwab=7.7498086336419e-10; KP_Visit=0.625837099162698; KP_Sid=0.0284246903301365"

x<-test_d+1 ###to avoid 0 
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)
mds_data$SampleID <- rownames(mds_data)
splitname<-strsplit(mds_data$SampleID, "[._]")

mds_data$Panelist=rep("NA", Clen)
mds_data$Visit=rep("NA", Clen)
mds_data$Swabber=rep("NA", Clen)
mds_data$studyID=rep("NA", Clen)
mds_data$BodySite=rep("NA", Clen)
mds_data$flowcell=rep("NA", Clen)
mds_data$scanOrder=rep("NA", Clen)

for(mm in  1:Clen ){
  mds_data$Panelist[mm]=splitname[[mm]][3]
  mds_data$Visit[mm]=splitname[[mm]][4]
  mds_data$Swabber[mm]=splitname[[mm]][5]
  mds_data$BodySite[mm]=splitname[[mm]][6]
  mds_data$scanOrder[mm]=splitname[[mm]][7]
  mds_data$flowcell[mm]=splitname[[mm]][2]
  mds_data$studyID[mm]=splitname[[mm]][1]
}
mds_data$SiteName=mds_data$BodySite
mds_data$SiteName[BodySite==1]="InnerThigh"
mds_data$SiteName[BodySite==2]="MonsPubis"
mds_data$SiteName[BodySite==3]="BikiniLegCrease"
mds_data$SiteName[BodySite==4]="OutsideLabiaMajora"
mds_data$SiteName[BodySite==5]="OutsideLabiaMinora"
mds_data$SiteName[BodySite==6]="Perineum"
mds_data$SiteName[BodySite==7]="Perianal"
mds_data$SiteName[BodySite==8]="Introitus"
mds_data$SiteName[BodySite==9]="Vaginal"
mds_data$SiteName[Panelist=="Zymo"]="Zymo"
mds_data$SiteName[is.na(Swabber)]="Other"
mds_data$Sid=paste0(mds_data$Panelist, "_", mds_data$Visit)
mds_data$SidSwab=as.factor(paste0(mds_data$Sid, "_", mds_data$Swabber))
mds_data$SidSwabSite=as.factor(paste0(mds_data$Sid, "_", mds_data$Swabber, "_", mds_data$SiteName))
mds_data$SiteSwab=as.factor(paste0(mds_data$SiteName, "_", mds_data$Swabber))
mds_data$VisitSwab=as.factor(paste0(mds_data$Visit, "_", mds_data$Swabber))
mds_data$Swabber[is.na(mds_data$Swabber)]="Other"


library(stringr)

library(ggplot2)
 
 png(filename=paste0(filename,".beta.Swabber.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Swabber)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Swabber, fill=Swabber), alpha=0.1, geom="polygon")
 print(p1)
 dev.off()
 
 png(filename=paste0(filename,".beta.Site.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = SiteName)) +geom_point(size=6)+theme_bw()+
   stat_chull(aes(color=SiteName, fill=SiteName), alpha=0.1, geom="polygon")
 print(p1)
 dev.off()
 
 png(filename=paste0(filename,".beta.Site.Swab.png"),  width=1200, height=280)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = SiteName)) +geom_point(size=3)+theme_bw()+
   stat_chull(aes(color=SiteName, fill=SiteName), alpha=0.05, geom="polygon")+ facet_grid(. ~ Swabber)
 print(p1)
 dev.off()
 
 png(filename=paste0(filename,".beta.Swabber.Site.png"),  width=800, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Swabber)) +geom_point(size=3)+theme_bw()+stat_chull(aes(color=Swabber, fill=Swabber), alpha=0.1, geom="polygon")+ facet_grid(. ~ SiteName)
 print(p1)
 dev.off()
######################################
bc.dist <- vegdist(x, method = "bray")
jaccard.dist <- vegdist(x, method = "jaccard")
mds2 <- metaMDS(jaccard.dist)
mds_data2 <- as.data.frame(mds2$points)
mds_data2$SampleID <- rownames(mds_data2)

mds_data2$Panelist=mds_data$Panelist
mds_data2$Visit=mds_data$Visit
mds_data2$Swabber=mds_data$Swabber
mds_data2$SiteName=mds_data$SiteName
mds_data2$SiteSwab=mds_data$SiteSwab
mds_data2$Sid=mds_data$Sid
mds_data2$SidSwab=mds_data$SidSwab
mds_data2$SidSwabSite=mds_data$SidSwabSite

png(filename=paste0(filename,".beta.Swabber.SiteName.jaccard.png"),  width=800, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = Swabber)) +geom_point(size=6)+theme_bw()+
  stat_chull(aes(color=Swabber, fill=Swabber), alpha=0.1, geom="polygon")+ facet_grid(. ~ SiteName)
print(p1)
dev.off()

# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram
png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=8000, height=2000)
plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
dev.off()

print("adonis result bc.dist ~ Swabber")
adonis2(bc.dist ~ Swabber)
#p=0.015 *

test.adonis <- adonis(bc.dist ~ SiteName)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis
#p=0.001
###################################################################
#PAIRWISE PERMANOVA
cbn <- combn(x=unique(SiteSwab), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(SiteSwab %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(SiteSwab %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$SiteSwab)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
# Output the results to a file
write.table(p.table, "beta_site_swab.txt", sep = "\t", quote = FALSE)




cbn <- combn(x=unique(SiteName), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(SiteName %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(SiteName %in% cbn[,i]), ]
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$SiteName)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
write.table(p.table, "beta_sitename.txt", sep = "\t", quote = FALSE)


cbn <- combn(x=unique(Swabber), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(Swabber %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(Swabber %in% cbn[,i]), ]
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$Swabber)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
write.table(p.table, "beta_swabber.txt", sep = "\t", quote = FALSE)
###########################Jaccard test
print("adonis result jaccard.dist ~ TreatVisit")
adonis2(jaccard.dist ~ SiteName)
test.adonis <- adonis(jaccard.dist ~ SiteName)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

test.adonis <- adonis(jaccard.dist ~ Swabber)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis
adonis2(jaccard.dist ~ SiteSwab)
#