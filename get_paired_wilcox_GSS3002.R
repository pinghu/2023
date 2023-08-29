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
#filename="GSS3002.metaphlan4.relab10.7"
filename="GSS3002_Sequence_reads_Analysis.txt"
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
  id2[mm]=splitname[[mm]][6]
  Visit[mm]=splitname[[mm]][3]
}
Collect[Collect=="Swab1"]="Expert"
Collect[Collect=="Swab3"]="Mom"
NewID=paste0(Collect, ".", BodySite, ".", Sid, ".", Visit)
BodySiteSid=paste0(BodySite, ".", Sid, ".", Visit)

SiteCollect=paste0(BodySite, Collect)
SiteCollectSid=paste0(BodySite, Collect, Sid)
SidCollect=paste0(Sid, Collect)


print(paste("Taxon","ShortName","percentZeroNA","Pwilcox_Expert1_Mom1","KP_Collect", "Pwilcox_Butt_BChest", "Pwilcox_Butt_Genital", "Pwilcox_BChest_Genital",
            "KP_BodySite", "KP_Sid", "Pwilcox_Butt_BLANK", "Pwilcox_BChest_BLANK", "Pwilcox_Genital_BLANK",
             "MExpert1", "MMom1", "MButt", "MBChest", "MGenital",  "truefc_Expert1_Mom1", "FC_Expert1_Mom1",sep=","))

for (i in 1:d[1]){
    genename=A[i, 1]
    if(is.na(genename)){next;}
    gene<-as.numeric(C[i,])
    splitG<-strsplit(as.character(genename), "[.|]")
    LLL=length(splitG[[1]])
    genus=splitG[[1]][LLL]
    genus1=str_replace(genus, "[kposgfc]__", "")
    genus2=str_replace(genus1, "_", " ")
    mydata=data.frame(gene,log(gene), Collect, Sid, BodySite, NewID, Visit, id2, BodySiteSid)
    mydata1=mydata[BodySite !="Blank",]
    mydata2=mydata1[BodySite != "Blank" & BodySite  != "Genital",]                   
    
    
    KP_Collect=kruskal.test(gene ~ Collect, data = mydata1)$p.value
    KP_BodySite=kruskal.test(gene ~ BodySite, data = mydata1)$p.value
    KP_Sid=kruskal.test(gene ~ Sid, data = mydata1)$p.value
    Expert<-gene[Collect=="Expert"]
    Expert1<-gene[Collect=="Expert"&BodySite != "BLANK" & BodySite !="Genital"]
    Mom1<-gene[Collect=="Mom"& BodySite != "BLANK"& BodySite !="Genital"]
    MExpert1=mean(Expert1);
    MMom1=mean(Mom1);
    
    Pwilcox_Expert1_Mom1<-my.wilcox.p.value(Expert1, Mom1, na.rm=TRUE,paired = TRUE, alternative = "two.sided")
    FC_Expert1_Mom1<-MExpert1/MMom1
    truefc_Expert1_Mom1=truefc(FC_Expert1_Mom1)
    Genital=gene[BodySite=="Genital"]; MGenital=mean(Genital)
    BLANK=gene[BodySite=="BLANK"]; MBLANK=mean(BLANK)
    Butt=gene[BodySite=="Butt"]; MButt=mean(Butt)
    BChest=gene[BodySite=="BChest"]; MBChest=mean(BChest)
    Pwilcox_Butt_Genital<-my.wilcox.p.value(Butt, Genital, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_Butt_BChest<-my.wilcox.p.value(Butt, BChest, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_BChest_Genital<-my.wilcox.p.value(BChest, Genital, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    
    Pwilcox_Butt_BLANK<-my.wilcox.p.value(Butt, BLANK, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_BChest_BLANK<-my.wilcox.p.value(BChest,BLANK, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_Genital_BLANK<-my.wilcox.p.value(Genital,BLANK, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    
    percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
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
      p2=ggplot(mydata2, aes(Collect, gene)) +
        geom_boxplot(width=0.3, size=1.5, fatten=1.5, colour="grey70") +
        geom_point(colour="red", size=2, alpha=0.5) +
        geom_line(aes(group=BodySiteSid), colour="red", linetype="11") +
        theme_classic()
      print(p2)
      dev.off()
      
      jpeg(file=paste0(genus2, ".collect-line.site.jpeg"))
      p2=ggplot(mydata1, aes(SiteCollect, gene)) +
        geom_boxplot(width=0.3, size=1.5, fatten=1.5, colour="grey70") +
        geom_point(colour="red", size=2, alpha=0.5) +
        geom_line(aes(group=BodySiteSid), colour="red", linetype="11") +
        theme_classic()
      print(p2)
      dev.off()
      
    }
    print(paste(genename,genus2,percentZeroNA, Pwilcox_Expert1_Mom1,KP_Collect, Pwilcox_Butt_BChest, Pwilcox_Butt_Genital, Pwilcox_BChest_Genital,
                KP_BodySite, KP_Sid, Pwilcox_Butt_BLANK, Pwilcox_BChest_BLANK, Pwilcox_Genital_BLANK,
                MExpert1, MMom1, MButt, MBChest, MGenital,  truefc_Expert1_Mom1, FC_Expert1_Mom1,sep=","))
}
