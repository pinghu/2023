rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
filename<- args[1]
filename="metaphlan.estimated_count.filtered.7"

test <- read.table(filename, header = TRUE, row.names = 1, sep="\t")
test[is.na(test)]<-0
d=dim(test)
d
##848,161
X=apply(test, 1, sum)
test_filter2=test[(X/d[2] >=10),]

test_d=data.frame(t(test))

min(apply(test_d, 1, sum)) ###162338


library(vegan)

test_a <- decostand(test_d, method = "total")

Cname=rownames(test_d)
Clen=length(Cname) ##there are 3 annotation columns
Site=rep("NA", Clen)
Country=rep("NA", Clen)
Sex=rep("NA", Clen)
SID=rep("NA", Clen)
Age=rep("NA",Clen)
RashScore=rep("NA", Clen)
SkinPH=rep("NA", Clen)
TEWL=rep("NA", Clen)
ID=rep("NA", Clen)
OID=rep("NA", Clen)
ScoreCat=rep("NA", Clen)
SiteCountry=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")
PHCat=rep("NA", Clen)

for(mm in  1:Clen ){
       Site[mm]=splitname[[mm]][1]
       Age[mm]=splitname[[mm]][5]
       Country[mm]=splitname[[mm]][2]
       Sex[mm]=splitname[[mm]][4]
       SID[mm]=splitname[[mm]][3]
       RashScore[mm]=splitname[[mm]][6]
       TEWL[mm]=splitname[[mm]][8]
       SkinPH[mm]=splitname[[mm]][7]
       OID[mm]=splitname[[mm]][9]
        
       if(RashScore[mm]=="NA"){ScoreCat[mm]="NoRash"}
       else if(RashScore[mm] ==0){ScoreCat[mm]="NoRash"}
       else if(RashScore[mm] >0){ScoreCat[mm]="Rash"}
       if(SkinPH[mm]=="NA"){PHCat[mm]="NoRash"}
       else if(SkinPH[mm] <=550){PHCat[mm]="LowPH"}
       else if(SkinPH[mm] >550){PHCat[mm]="HighPH"}
       
}
SkinPH=as.numeric(SkinPH)/100
TEWL=as.numeric(TEWL)/10
RashScoreN=as.numeric(RashScore)/10
RashScore <- factor(RashScore,levels = c("0", "5", "10", "15", "20"))
RashCat=ScoreCat
RashCat[RashScoreN>1]="HighRash(>=1.5)"
RashCat[RashCat=="Rash"]="MildRash(0.5,1)"

RashPH=paste0(RashCat, ".", PHCat)
ScorePH=paste0(ScoreCat, ".", PHCat)
library(ggplot2)

#https://stackoverflow.com/questions/30057765/histogram-ggplot-show-count-label-for-each-bin-for-each-category

shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")

observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)

mydata=data.frame(shannon,simp,invsimp, observed, Menhinick_index, Margalef_index, 
                   Site,as.numeric(Age),RashScore, RashScoreN,ScoreCat,
                  Sex, Country, 
                 TEWL,SkinPH,PHCat, RashCat,RashPH, ScorePH,
                  SID, OID)

write.table(t(mydata), file =paste0(filename, "alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

saveRDS(mydata,paste0(filename, ".diversitydata"))
saveRDS(test_d,paste0(filename,".test_d"))

mydata<-mydata[mydata$Site!="Site",]
mydata2<-mydata[mydata$RashScore !="NA",]
mydata2$RashScore <- factor(mydata2$RashScore,levels = c("0", "5", "10", "15", "20"))

library(ggplot2)

#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
library(tidyverse)
library(ggpubr)
library(rstatix)



png(filename=paste0(filename,"RashScoreN.PHCat.png"), width=380, height=600)
stat.test <- mydata %>%
  t_test(RashScoreN ~ PHCat) %>%
  mutate(y.position = 2.5)
bxp <- ggboxplot(mydata, x = "PHCat", y = "RashScoreN", color = "PHCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"RashScoreN.RashPH.png"), width=800, height=600)
stat.test <- mydata %>%
  t_test(RashScoreN ~ RashPH) %>%
  mutate(y.position = 2.5)
bxp <- ggboxplot(mydata, x = "RashPH", y = "RashScoreN", color = "RashPH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"RashScoreN.ScorePH.png"), width=600, height=600)
stat.test <- mydata %>%
  t_test(RashScoreN ~ ScorePH) %>%
  mutate(y.position = 2.5)
bxp <- ggboxplot(mydata, x = "ScorePH", y = "RashScoreN", color = "ScorePH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()


png(filename=paste0(filename,"SkinPH.RashPH.png"), width=800, height=600)
stat.test <- mydata %>%
  t_test(SkinPH ~ RashPH) %>%
  mutate(y.position = 7.5)
bxp <- ggboxplot(mydata, x = "RashPH", y = "SkinPH", color = "RashPH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"SkinPH.ScorePH.png"), width=600, height=600)
stat.test <- mydata %>%
  t_test(SkinPH ~ ScorePH) %>%
  mutate(y.position = 7.5)
bxp <- ggboxplot(mydata, x = "ScorePH", y = "SkinPH", color = "ScorePH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()


png(filename=paste0(filename,"SkinPH.RashScoreCat.png"), width=380, height=600)
stat.test <- mydata %>%
  t_test(SkinPH ~ ScoreCat) %>%
  mutate(y.position = 7.5)
bxp <- ggboxplot(mydata, y = "SkinPH", x = "ScoreCat", color = "ScoreCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"SkinPH.RashCat.png"), width=380, height=600)
stat.test <- mydata %>%
  t_test(SkinPH ~ RashCat) %>%
  mutate(y.position = 7.5)
bxp <- ggboxplot(mydata, y = "SkinPH", x = "RashCat", color = "RashCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()


##########################################################
png(filename=paste0(filename,"shannon.PHCat.png"),  width=380, height=600)
stat.test <- mydata %>%
  t_test(shannon ~ PHCat) %>%
  mutate(y.position = 4.6)
bxp <- ggboxplot(mydata, x = "PHCat", y = "shannon", color = "PHCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()


png(filename=paste0(filename,"shannon.ScoreCat.png"),  width=380, height=600)
stat.test <- mydata %>%
  t_test(shannon ~ ScoreCat) %>%
  mutate(y.position = 4.6)
bxp <- ggboxplot(mydata, x = "ScoreCat", y = "shannon", color = "ScoreCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"shannon.RashCat.png"),  width=380, height=600)
stat.test <- mydata %>%
  t_test(shannon ~ RashCat) %>%
  mutate(y.position = 4.6)
bxp <- ggboxplot(mydata, x = "RashCat", y = "shannon", color = "RashCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"shannon.RashPH.png"),  width=800, height=600)
stat.test <- mydata %>%
  t_test(shannon ~ RashPH) %>%
  mutate(y.position = 4.6)
bxp <- ggboxplot(mydata, x = "RashPH", y = "shannon", color = "RashPH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"shannon.ScorePH.png"),  width=600, height=600)
stat.test <- mydata %>%
  t_test(shannon ~ ScorePH) %>%
  mutate(y.position = 4.6)
bxp <- ggboxplot(mydata, x = "ScorePH", y = "shannon", color = "ScorePH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()





png(filename=paste0(filename,"observed.PHCat.png"), width=380, height=600)
stat.test <- mydata %>%
  t_test(observed ~ PHCat) %>%
  mutate(y.position = 250)
bxp <- ggboxplot(mydata, x = "PHCat", y = "observed", color = "PHCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"observed.RashCat.png"), width=380, height=600)
stat.test <- mydata %>%
  t_test(observed ~ RashCat) %>%
  mutate(y.position = 250)
bxp <- ggboxplot(mydata, x = "RashCat", y = "observed", color = "RashCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"observed.ScoreCat.png"), width=380, height=600)
stat.test <- mydata %>%
  t_test(observed ~ ScoreCat) %>%
  mutate(y.position = 250)
bxp <- ggboxplot(mydata, x = "ScoreCat", y = "observed", color = "ScoreCat")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"observed.ScorePH.png"), width=600, height=600)
stat.test <- mydata %>%
  t_test(observed ~ ScorePH) %>%
  mutate(y.position = 250)
bxp <- ggboxplot(mydata, x = "ScorePH", y = "observed", color = "ScorePH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"observed.RashPH.png"), width=800, height=600)
stat.test <- mydata %>%
  t_test(observed ~ RashPH) %>%
  mutate(y.position = 250)
bxp <- ggboxplot(mydata, x = "RashPH", y = "observed", color = "RashPH")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

#########################################################


png(filename=paste0(filename,"shannon.country.png"))
stat.test2 <- mydata %>%
  t_test(shannon ~ Country) %>%
  mutate(y.position = 4.6)
bxp <- ggboxplot(mydata, x = "Country", y = "shannon", color = "Country")
bxp + stat_pvalue_manual(stat.test2,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter()
dev.off()

png(filename=paste0(filename,"shannon.SkinPH-RashCat-Country.png"),width=1000, height = 380)
bxp <- ggboxplot(mydata, x = "RashCat", y = "SkinPH", color = "RashCat", facet.by = "Country")
bxp + geom_jitter()
dev.off()

png(filename=paste0(filename,"shannon.SkinPH-ScoreCat-Country.png"),width=600, height = 380)
bxp <- ggboxplot(mydata, x = "ScoreCat", y = "SkinPH", color = "ScoreCat", facet.by = "Country")
bxp + geom_jitter()
dev.off()

png(filename=paste0(filename,"shannon.RashScoreN-PHCat-Country.png"),width=1000, height = 380)
bxp <- ggboxplot(mydata, x = "PHCat", y = "RashScoreN", color = "PHCat", facet.by = "Country")
bxp + geom_jitter()
dev.off()


KP_PHCat=kruskal.test(shannon ~ PHCat, data = mydata)$p.value
KP_ScoreCat=kruskal.test(shannon ~ ScoreCat, data = mydata)$p.value
KP_RashCat=kruskal.test(shannon ~ RashCat, data = mydata)$p.value
KP_RashPH=kruskal.test(shannon ~ RashPH, data = mydata)$p.value
KP_ScorePH=kruskal.test(shannon ~ ScorePH, data = mydata)$p.value
KP_Country=kruskal.test(shannon ~ Country, data = mydata)$p.value

print(paste0("KP_PHCat=", KP_PHCat, "; KP_Country=", KP_Country, "; KP_RashCat=", KP_RashCat,
             "; KP_ScoreCat=", KP_ScoreCat, "; KP_RashPH=", KP_RashPH, "; KP_ScorePH=", KP_ScorePH
     ))

########################################
x<-test_d
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- cbind(as.data.frame(mds$points), mydata)
mds_data$SampleID <- rownames(mds_data)

#install.packages('devtools')
library(devtools)
#install_github('fawda123/ggord')
library(ggord)
library(ggplot2)
# https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
 #png(filename=paste0(filename,".beta.PHCat.png"),  width=600, height=380)
 #p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = PHCat)) +geom_point()+theme_bw()+geom_text(aes(label=Site),hjust=0.5, vjust=0)
 #print(p1)
 #dev.off() 
 
 png(filename=paste0(filename,".beta.PHCat.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = PHCat)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=PHCat, fill=PHCat), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 
 
 png(filename=paste0(filename,".beta.RashCat.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = RashCat)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=RashCat, fill=RashCat), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 
 
 png(filename=paste0(filename,".beta.ScoreCat.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = ScoreCat)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=ScoreCat, fill=ScoreCat), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 
 
 png(filename=paste0(filename,".beta.RashPH.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = RashPH)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=RashPH, fill=RashPH), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 
 
 png(filename=paste0(filename,".beta.ScorePH.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = ScorePH)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=ScorePH, fill=ScorePH), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 
 
 png(filename=paste0(filename,".beta.Country.png"),  width=600, height=380)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Country)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=Country, fill=Country), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 
 
 #png(filename=paste0(filename,".beta.SID.png"),  width=600, height=380)
 #p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = SID)) +geom_point(size=6)+theme_bw()+stat_chull(aes(color=SID, fill=SID), alpha=0.1, geom="polygon")+ theme(legend.position = "none")
 #print(p1)
 #dev.off() 
 
############################################################
print("all samples")
# calculate Bray-Curtis distance among samples
bc.dist <- vegdist(test_d, method = "bray")
# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram
png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=2000, height=800)
plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
dev.off()
# Taxonomic (Bray-Curtis) dissimilarity explained
print("adonis result bc.dist ~ ScoreCat")
adonis2(bc.dist ~ ScoreCat)
adonis2(bc.dist ~ PHCat)
adonis2(bc.dist ~ RashCat)
adonis2(bc.dist ~ RashPH)
adonis2(bc.dist ~ ScorePH)
adonis2(bc.dist ~ Country)
####P=0.097 ScoreCat|| Sex 0.024 ||Score 1: 0.26| Score 2 0.088 | ScoreN 0.044 ######

####without 3408 ScoreCat:0.099 ||Sex 0.025 ||Score 1 0.286 | Score 2 0.104 |score N 0.037
print("adonis  result bc.dist ~ Sex")
adonis2(bc.dist ~ Sex)


print("adonis  result bc.dist ~ RashScore")
adonis2(bc.dist ~ RashScoreN)

######################################
jaccard.dist <- vegdist(x, method = "jaccard")
mds2 <- metaMDS(jaccard.dist)
mds_data2 <- as.data.frame(cbind(mds2$points, mydata))
mds_data2$SampleID <- rownames(mds_data2)

png(filename=paste0(filename,".beta.PHCat.jaccard.png"),  width=380, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = PHCat)) +geom_point(size=6)+theme_bw()+
  stat_chull(aes(color=PHCat, fill=PHCat), alpha=0.1, geom="polygon") 
print(p1)
dev.off()

png(filename=paste0(filename,".beta.PHCat.RashCat.jaccard.png"),  width=800, height=380)
p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = PHCat)) +geom_point(size=6)+theme_bw()+
  stat_chull(aes(color=PHCat, fill=PHCat), alpha=0.1, geom="polygon") + facet_grid(. ~ RashCat)
print(p1)
dev.off()

# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram
png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=2000, height=800)
plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
dev.off()


test.adonis <- adonis2(bc.dist ~ RashPH)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis
#p=0.001
###################################################################
#PAIRWISE PERMANOVA
cbn <- combn(x=unique(RashPH), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(RashPH %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(RashPH %in% cbn[,i]), ]
  
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$RashPH)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
# Output the results to a file
write.table(p.table, "beta_RashPH.txt", sep = "\t", quote = FALSE)


cbn <- combn(x=unique(ScorePH), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(ScorePH %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(ScorePH %in% cbn[,i]), ]
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$ScorePH)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
write.table(p.table, "beta_ScorePH.txt", sep = "\t", quote = FALSE)


cbn <- combn(x=unique(Country), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- x [c(Country %in% cbn[,i]), ]
  metadata_sub <- mds_data[c(Country %in% cbn[,i]), ]
  permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$Country)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
write.table(p.table, ".Country.txt", sep = "\t", quote = FALSE)
###########################Jaccard test
print("adonis result jaccard.dist ~ RashPH")
adonis2(jaccard.dist ~ RashPH)
test.adonis <- adonis2(jaccard.dist ~ RashPH)
test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

test.adonis <- adonis2(jaccard.dist ~ ScorePH)
#test.adonis <- as.data.frame(test.adonis$aov.tab)
test.adonis

#



