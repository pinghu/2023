##########################
rm(list=ls())
library(stringr)
library("plotrix")
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggplot2)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se= std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
#filename="GSS3001.metaphlan4.relab10.7"
filename="GSS3001_Sequence_reads_Analysis.txt"

library(ggplot2)
library("ggpubr")
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
truefc<-function(VVV){
	XXX=VVV
	if(VVV==0){
	    XXX=NA
   	}else if(VVV<1){
	    XXX=-1/VVV
    	}
	return(XXX)
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=A[1:d[1], 2:d[2]]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ

Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
Collect=rep("NA", Clen)
Sid=rep("NA", Clen)
BodySite=rep("NA", Clen)
id2=rep("NA", Clen)
splitname<-strsplit(Cname, "_")
for(mm in  1:Clen ){
  Collect[mm]=splitname[[mm]][5]
  Sid[mm]=splitname[[mm]][2]
  BodySite[mm]=splitname[[mm]][3]
  id2[mm]=splitname[[mm]][4]
}
Collect[35]="Swab1"
Collect[36]="Tape1"
for( i in c(1,2,3,10)) {Collect[i]="Control"}
for( i in 11:26) {BodySite[i]=id2[i]}
for( i in 33:35) {BodySite[i]=id2[i]}
BodySite[36]="SkinMC"
BodySite[1]="SkinMC"
BodySite[2]=id2[2]
BodySite[10]=id2[2]
for (i in c(1,2,3,10,35,36)){Sid[i]=BodySite[i]}
NewID=paste0(Collect, ".", BodySite, ".", Sid)
for (i in 4:9){NewID[i]=paste0(NewID[i],id2[i])}
for (i in 27:32){NewID[i]=paste0(NewID[i],id2[i])}
SiteSID=paste0(BodySite, ".", Sid)
for (i in 4:9){SiteSID[i]=paste0(SiteSID[i],id2[i])}
for (i in 27:32){SiteSID[i]=paste0(SiteSID[i],id2[i])}
SiteCollect=paste0(BodySite, Collect)
SiteCollectSid=paste0(BodySite, Collect, Sid)
SidCollect=paste0(Sid, Collect)

print(paste("Taxon","ShortName","percentZeroNA",
            "Pwilcox_Swab_Tape","Pwilcox_Swab1_Tape1_noSkinMC","Pwilcox_Swab2_Tape2_babyOnly","KP_Collect",
            "Pwilcox_Butt_BChest", "Pwilcox_Butt_Genital", "Pwilcox_BChest_Genital","KP_BodySite",
            "KP_Sid", "Pwilcox_Butt_BLANK", "Pwilcox_BChest_BLANK", "Pwilcox_Genital_BLANK", "Pwilcox_Butt_Reagent", "Pwilcox_BChest_Reagent", "Pwilcox_Genital_Reagent",
            "MSwab", "MTape", "MSwab1", "MTape1","MSwab2", "MTape2","MSkinMC","MReagents", "MZymo", "MArm", "MBChest", "MBlank", "MButt", "MGenital", 
            "truefc_Swab_Tape", "truefc_Swab1_Tape1", "truefc_Swab2_Tape2",sep=","))

for (i in 1:d[1]){
    genename=A[i, 1]
    if(is.na(genename)){next;}
    gene<-as.numeric(C[i,])
    splitG<-strsplit(as.character(genename), "[.|]")
    LLL=length(splitG[[1]])
    genus=splitG[[1]][LLL]
    genus1=str_replace(genus, "[kposgfc]__", "")
    genus2=str_replace(genus1, "_", " ")
    mydata=data.frame(gene,log(gene), BodySite, SiteSID, Collect, id2, Sid, SiteCollect, SiteCollectSid, SidCollect)
    mydata2=mydata[BodySite == "BChest" | BodySite == "Genital" | BodySite == "Butt",]
    mydata1=mydata[BodySite == "BChest" | BodySite == "Genital" | BodySite == "Butt" | BodySite== "Arm",]
    KP_Collect=kruskal.test(gene ~ Collect, data = mydata1)$p.value
    KP_BodySite=kruskal.test(gene ~ BodySite, data = mydata1)$p.value
    KP_Sid=kruskal.test(gene ~ Sid, data = mydata1)$p.value
    Swab<-mydata$gene[mydata$Collect=="Swab1"]
    Tape<-mydata$gene[mydata$Collect=="Tape1"]
    MSwab=mean(Swab);
    MTape=mean(Tape);
    FC_Swab_Tape<-MSwab/MTape
    truefc_Swab_Tape=truefc(FC_Swab_Tape)
    Pwilcox_Swab_Tape<-my.wilcox.p.value(Swab, Tape, na.rm=TRUE,paired = TRUE, alternative = "two.sided")
    Swab1<-mydata1$gene[mydata1$Collect=="Swab1"]
    Tape1<-mydata1$gene[mydata1$Collect=="Tape1"]
    
    
    
    MSwab1=mean(Swab1);
    MTape1=mean(Tape1);
    
    Pwilcox_Swab1_Tape1<-my.wilcox.p.value(Swab1, Tape1, na.rm=TRUE,paired = TRUE, alternative = "two.sided")
    FC_Swab1_Tape1<-MSwab1/MTape1
    truefc_Swab1_Tape1=truefc(FC_Swab1_Tape1)
    Swab2<-mydata2$gene[mydata2$Collect=="Swab1"]
    Tape2<-mydata2$gene[mydata2$Collect=="Tape1"]
    
  
    MSwab2=mean(Swab2);
    MTape2=mean(Tape2);
    
    Pwilcox_Swab2_Tape2<-my.wilcox.p.value(Swab2, Tape2, na.rm=TRUE,paired = TRUE, alternative = "two.sided")
    FC_Swab2_Tape2<-MSwab2/MTape2
    truefc_Swab2_Tape2=truefc(FC_Swab2_Tape2)   
  
    percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
    
    Genital=gene[BodySite=="Genital"]; MGenital=mean(Genital)
    BLANK=gene[BodySite=="Blank"]; MBLANK=mean(BLANK)
    Butt=gene[BodySite=="Butt"]; MButt=mean(Butt)
    BChest=gene[BodySite=="BChest"]; MBChest=mean(BChest)
    Reagent=gene[BodySite=="reagents"]; MReagent=mean(Reagent)
    SkinMC=gene[BodySite=="SkinMC"]; MSkinMC=mean(SkinMC)
    Zymo=gene[BodySite=="Zymo"]; MZymo=mean(Zymo)
    Arm=gene[BodySite=="Arm"]; MArm=mean(Arm)
    
    Pwilcox_Butt_Genital<-my.wilcox.p.value(Butt, Genital, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_Butt_BChest<-my.wilcox.p.value(Butt, BChest, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_BChest_Genital<-my.wilcox.p.value(BChest, Genital, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    
    Pwilcox_Butt_BLANK<-my.wilcox.p.value(Butt, BLANK, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_BChest_BLANK<-my.wilcox.p.value(BChest,BLANK, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_Genital_BLANK<-my.wilcox.p.value(Genital,BLANK, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    
    Pwilcox_Butt_Reagent<-my.wilcox.p.value(Butt, Reagent, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_BChest_Reagent<-my.wilcox.p.value(BChest,Reagent, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_Genital_Reagent<-my.wilcox.p.value(Genital, Reagent, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    
    if(percentZeroNA<0.75){
      
      png(filename=paste0(genus2,".Collect.png"),width=280, height = 480)
      stat.test <- mydata1 %>%
        wilcox_test(gene ~ Collect) %>%
        mutate(y.position = max(gene))
      p0 <- ggboxplot(mydata1, x = "Collect", y = "gene", color = "Collect") + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
      print(p0)
      dev.off()
      
      
      png(filename=paste0(genus2,"BodySite.png"), width=480, height = 480)
      #stat.test <- mydata %>%
       # wilcox_test(gene ~ BodySite) %>%
        #mutate(y.position = max(gene))
      p1 <- ggboxplot(mydata1, x = "BodySite", y = "gene", color = "BodySite")+ stat_compare_means()
      print(p1)   
      #+ stat_pvalue_manual(stat.test2,   label = "p", tip.length = 0.01, step.increase = 0.08)
      #print(p1)
      dev.off()
      
      png(filename=paste0(genus2,".Collect.Site.png"), height = 600, width=600)
      p2 <- ggboxplot(mydata1, x = "Collect", y = "gene", color = "Collect", facet.by = "BodySite")+ stat_compare_means(comparisons = list(c("Swab1", "Tape1")), method = "wilcox")
      print(p2)   
      dev.off()
      
      jpeg(file=paste0(genus2, ".collect-line.jpeg"))
      p2=ggplot(mydata1, aes(Collect, gene)) +
        geom_boxplot(width=0.3, size=1.5, fatten=1.5, colour="grey70") +
        geom_point(colour="red", size=2, alpha=0.5) +
        geom_line(aes(group=SiteSID), colour="red", linetype="11") +
        theme_classic()
      print(p2)
      dev.off()
      
      jpeg(file=paste0(genus2, ".collect-line.site.jpeg"))
      p2=ggplot(mydata1, aes(SiteCollect, gene)) +
        geom_boxplot(width=0.3, size=1.5, fatten=1.5, colour="grey70") +
        geom_point(colour="red", size=2, alpha=0.5) +
        geom_line(aes(group=SiteSID), colour="red", linetype="11") +
        theme_classic()
      print(p2)
      dev.off()
      
      
      
    }
    print(paste(genename,genus2,percentZeroNA, Pwilcox_Swab_Tape,Pwilcox_Swab1_Tape1,Pwilcox_Swab2_Tape2,KP_Collect, 
                Pwilcox_Butt_BChest, Pwilcox_Butt_Genital, Pwilcox_BChest_Genital, KP_BodySite, 
                KP_Sid, Pwilcox_Butt_BLANK, Pwilcox_BChest_BLANK, Pwilcox_Genital_BLANK, Pwilcox_Butt_Reagent, Pwilcox_BChest_Reagent, Pwilcox_Genital_Reagent,
                MSwab, MTape, MSwab1, MTape1, MSwab2, MTape2,MSkinMC, MReagent, MZymo, MArm, MBChest, MBLANK, MButt, MGenital, truefc_Swab_Tape, truefc_Swab1_Tape1, truefc_Swab2_Tape2,sep=","))
    
    }
