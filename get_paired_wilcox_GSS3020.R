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
#filename="metaphlan.relab10.7"
#filename="pathabundance.community.relab10"
filename="GSS3020_read_count_data.txt"
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

HeightInch=rep("NA", Clen)
Weightlbs=rep("NA", Clen)
BMI=rep("NA", Clen)
SID=rep("NA", Clen)
BodySite=rep("NA", Clen)
Ethinicity=rep("NA", Clen)
Age=rep("NA", Clen)
Disease=rep("NA", Clen)
for(mm in  1:Clen ){
  HeightInch[mm]=splitname[[mm]][2]
  Weightlbs[mm]=splitname[[mm]][3]
  BMI[mm]=splitname[[mm]][4]
}
splitname<-strsplit(Cname, "[.]")
for(mm in  1:Clen ){
  BodySite[mm]=splitname[[mm]][1]
  Disease[mm]=splitname[[mm]][2]
  Ethinicity[mm]=splitname[[mm]][3]
  SID[mm]=splitname[[mm]][4]
  Age[mm]=splitname[[mm]][5]
}
DiseaseEthBodySite=paste0(Disease, "_", Ethinicity, "_", BodySite)
DiseaseEthBodySiteSID=paste0(Disease, "_", Ethinicity, "_", BodySite, "_", SID)
DiseaseBodySite=paste0(Disease, "_",  BodySite)
EthBodySite=paste0(Ethinicity, "_", BodySite)
EthDisease=paste0(Ethinicity, "_", Disease)

print(paste("genename","genus2","percentZeroNA","KP_Disease","KP_DiseaseChinese","Pwilcox_PostAcne_Healthy","Pwilcox_ChinChineseH_ChinChineseP","Pwilcox_CheekChineseH_CheekChineseP",
            "Pwilcox_CheekChineseH_CheekAAP","Pwilcox_CheekChineseH_CheekCaucasianP","Pwilcox_ChinChineseH_ChinAAP","Pwilcox_ChinChineseH_ChinCaucasianP", 
            "Pwilcox_ChineseH_AAP","Pwilcox_ChineseH_CaucasianP","Pwilcox_ChineseH_ChineseP",
            "truefc_PostAcne_Healthy","truefc_ChineseP_ChineseH","turefc_ChinChineseP_ChinChineseH","truefc_CheekChineseP_CheekChineseH", 
            "MChineseH","MChineseP","MHealthy","MPostAcne","MChinChineseH","MChinChineseP","MChinCaucasianP","MChinAAP","MCheekChineseH","MCheekChineseP","MCheekCaucasianP","MCheekAAP",
            "KP_Ethinicity","KP_DiseaseEthinicity","KP_EthinicityDisease", 
            "Pwilcox_CheekChineseP_CheekAAP","Pwilcox_CheekChineseP_CheekCaucasianP","Pwilcox_CheekAAP_CheekCaucasianP",
            "Pwilcox_ChinChineseP_ChinAAP","Pwilcox_ChinChineseP_ChinCaucasianP","Pwilcox_ChinAAP_ChinCaucasianP",
            "Pwilcox_AA_Caucasian","Pwilcox_AA_Chinese","Pwilcox_Caucasian_Chinese","Pwilcox_AAP_CaucasianP","Pwilcox_CheekAAP_CheekCaucasianP","Pwilcox_ChineseP_AAP","Pwilcox_ChineseP_CaucasianP",
            "MChinese","MAA","MCaucasian",
            "KP_BodySite","KP_BodySiteChinese","Pwilcox_Cheek_Chin","Pwilcox_CheekP_ChinP","Pwilcox_CheekChinese_ChinChinese",
            "Pwilcox_ChinCaucasianP_CheekCaucasianP","Pwilcox_ChinAAp_CheekAAP","Pwilcox_ChinChineseP_CheekChineseP","Pwilcox_ChinChineseH_CheekChineseH",
            "MChin","MCheek","MChinChinese","MCheekChinese","MChinP","MCheekP",
            "truefc_Cheek_Chin","truefc_CheekChinese_ChinChinese","truefc_CheekP_ChinP",
            "KP_Sid","KP_BodySiteEthinicityDisease","KP_DiseaseBodysiteChinese","KP_BodySiteDisease","KP_BodySiteEthinicity","KP_DiseaseBodysite","KP_DiseaseBodysiteEthinicity",
            sep=","))

for (i in 1:d[1]){
    genename=A[i, 1]
    if(is.na(genename)){next;}
    gene<-as.numeric(C[i,])
    splitG<-strsplit(as.character(genename), "[.|]")
    LLL=length(splitG[[1]])
    genus=splitG[[1]][LLL]
    genus1=str_replace(genus, "[kposgfc]__", "")
    genus2=str_replace(genus1, "_", " ")
    mydata=data.frame(gene,log(gene), BodySite, Disease, SID, Ethinicity, Age, BMI, Weightlbs, HeightInch, DiseaseEthBodySite, DiseaseEthBodySiteSID, DiseaseBodySite, EthBodySite, EthDisease)
    mydataDisease=mydata[Disease != "HLY",]
    mydataChinese=mydata[Ethinicity =="Chinese",]
    KP_DiseaseChinese=kruskal.test(gene ~ Disease, data = mydataChinese)$p.value
    KP_Disease=kruskal.test(gene ~ Disease, data = mydata)$p.value
    KP_BodySite=kruskal.test(gene ~BodySite, data = mydata)$p.value
    KP_BodySiteChinese =kruskal.test(gene ~BodySite, data = mydataChinese)$p.value
    KP_BodySiteDisease=kruskal.test(gene ~BodySite, data = mydataDisease)$p.value
    KP_Ethinicity=kruskal.test(gene ~ Ethinicity, data = mydata)$p.value
    KP_EthinicityDisease=kruskal.test(gene ~ Ethinicity, data = mydataDisease)$p.value
    KP_Sid=kruskal.test(gene ~SID, data = mydata)$p.value
    KP_DiseaseBodysite=kruskal.test(gene ~ DiseaseBodySite, data = mydata)$p.value
    KP_DiseaseBodysiteChinese=kruskal.test(gene ~ DiseaseBodySite, data = mydataChinese)$p.value
    KP_DiseaseEthinicity=kruskal.test(gene ~ EthDisease, data = mydata)$p.value
    KP_BodySiteEthinicity=kruskal.test(gene ~ EthBodySite, data = mydata)$p.value
    KP_BodySiteEthinicityDisease=kruskal.test(gene ~ EthBodySite, data = mydataDisease)$p.value
    KP_DiseaseBodysiteEthinicity=kruskal.test(gene ~ DiseaseEthBodySite, data = mydata)$p.value

    PostAcne<-gene[Disease=="PIH"]
    Healthy<-gene[Disease=="HLY"]
    MPostAcne=mean(PostAcne);
    MHealthy=mean(Healthy);
    FC_PostAcne_Healthy<-MPostAcne/MHealthy
    truefc_PostAcne_Healthy=truefc(FC_PostAcne_Healthy)
    Pwilcox_PostAcne_Healthy<-my.wilcox.p.value(PostAcne, Healthy, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
############################################################################
    Cheek<-gene[BodySite=="LeftCheek"]
    Chin<-gene[BodySite=="Chin"]
    MCheek=mean(Cheek);
    MChin=mean(Chin);
    FC_Cheek_Chin<-MCheek/MChin
    truefc_Cheek_Chin=truefc(FC_Cheek_Chin)
    Pwilcox_Cheek_Chin<-my.wilcox.p.value(Cheek, Chin, na.rm=TRUE,paired = TRUE, alternative = "two.sided")
##############################################################################
    Chinese <-gene[Ethinicity=="Chinese"]; MChinese=mean(Chinese)
    AA <-gene[Ethinicity=="AA"]; MAA=mean(AA)
    Caucasian <-gene[Ethinicity=="Caucasian"]; MCaucasian=mean(Caucasian)
    Pwilcox_AA_Chinese<-my.wilcox.p.value(AA, Chinese, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_AA_Caucasian<-my.wilcox.p.value(AA, Caucasian, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_Caucasian_Chinese<-my.wilcox.p.value(Caucasian, Chinese, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
#################################################################################
    ChineseH<-gene[EthDisease=="Chinese_HLY"]; MChineseH=mean(ChineseH);
    ChineseP<-gene[EthDisease=="Chinese_PIH"]; MChineseP=mean(ChineseP);
    AAP<-gene[EthDisease=="AA_PIH"]; MAAP=mean(AAP);
    CaucasianP <-gene[EthDisease=="Caucasian_PIH"]; MCaucasianP=mean(CaucasianP)
    Pwilcox_ChineseH_ChineseP<-my.wilcox.p.value(ChineseH, ChineseP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChineseH_CaucasianP<-my.wilcox.p.value(ChineseH, CaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChineseH_AAP<-my.wilcox.p.value(ChineseH, AAP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChineseP_CaucasianP<-my.wilcox.p.value(ChineseP, CaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChineseP_AAP<-my.wilcox.p.value(ChineseP, AAP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_AAP_CaucasianP<-my.wilcox.p.value(AAP, CaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
###########################################################################################################################
    ChinChineseH<-gene[DiseaseEthBodySite=="HLY_Chinese_Chin"]; MChinChineseH=mean(ChinChineseH)
    ChinChineseP<-gene[DiseaseEthBodySite =="PIH_Chinese_Chin"]; MChinChineseP = mean(ChinChineseP)
    ChinAAP<-gene[DiseaseEthBodySite=="PIH_AA_Chin"]; MChinAAP=mean(ChinAAP);
    ChinCaucasianP <-gene[DiseaseEthBodySite == "PIH_Caucasian_Chin"]; MChinCaucasianP= mean(ChinCaucasianP)
    CheekChineseH<-gene[DiseaseEthBodySite=="HLY_Chinese_LeftCheek"]; MCheekChineseH =mean(CheekChineseH)
    CheekChineseP<-gene[DiseaseEthBodySite =="PIH_Chinese_LeftCheek"]; MCheekChineseP =mean(CheekChineseP)
    CheekAAP<-gene[DiseaseEthBodySite=="PIH_AA_LeftCheek"]; MCheekAAP=mean(CheekAAP)
    CheekCaucasianP <-gene[DiseaseEthBodySite == "PIH_Caucasian_LeftCheek"]; MCheekCaucasianP =mean(CheekCaucasianP)
    
    Pwilcox_ChinChineseH_ChinChineseP<-my.wilcox.p.value(ChinChineseH, ChinChineseP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinChineseH_ChinCaucasianP<-my.wilcox.p.value(ChinChineseH, ChinCaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinChineseH_ChinAAP<-my.wilcox.p.value(ChinChineseH, ChinAAP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinChineseP_ChinCaucasianP<-my.wilcox.p.value(ChinChineseP, ChinCaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinChineseP_ChinAAP<-my.wilcox.p.value(ChinChineseP, ChinAAP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinAAP_ChinCaucasianP<-my.wilcox.p.value(ChinAAP, ChinCaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
 
    Pwilcox_CheekChineseH_CheekChineseP<-my.wilcox.p.value(CheekChineseH, CheekChineseP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_CheekChineseH_CheekCaucasianP<-my.wilcox.p.value(CheekChineseH, CheekCaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_CheekChineseH_CheekAAP<-my.wilcox.p.value(CheekChineseH, CheekAAP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_CheekChineseP_CheekCaucasianP<-my.wilcox.p.value(CheekChineseP, CheekCaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_CheekChineseP_CheekAAP<-my.wilcox.p.value(CheekChineseP, CheekAAP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_CheekAAP_CheekCaucasianP<-my.wilcox.p.value(CheekAAP, CheekCaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    
    Pwilcox_ChinChineseH_CheekChineseH<-my.wilcox.p.value(ChinChineseH, CheekChineseH, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinChineseP_CheekChineseP<-my.wilcox.p.value(ChinChineseP, CheekChineseP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinAAp_CheekAAP<-my.wilcox.p.value(ChinAAP, CheekAAP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    Pwilcox_ChinCaucasianP_CheekCaucasianP<-my.wilcox.p.value(ChinCaucasianP, CheekCaucasianP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
   
    percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
    turefc_ChinChineseP_ChinChineseH = truefc(MChinChineseP/MChinChineseH)
    truefc_CheekChineseP_CheekChineseH =truefc(MCheekChineseP/MCheekChineseH)
    truefc_ChineseP_ChineseH=truefc(MChineseP/MChinChineseH)
    
    CheekChinese<-gene[EthBodySite=="Chinese_LeftCheek"]; MCheekChinese=mean(CheekChinese)
    ChinChinese<-gene[EthBodySite=="Chinese_Chin"]; MChinChinese=mean(ChinChinese)
    truefc_CheekChinese_ChinChinese=truefc(MCheekChinese/MChinChinese)
    Pwilcox_CheekChinese_ChinChinese<-my.wilcox.p.value(CheekChinese, ChinChinese, na.rm=TRUE,paired = FALSE, alternative = "two.sided")
    
    CheekP<-gene[DiseaseBodySite=="PIH_Chin"]; MCheekP=mean(CheekP)
    ChinP<-gene[DiseaseBodySite=="PIH_LeftCheek"]; MChinP=mean(ChinP)
    truefc_CheekP_ChinP=truefc(MCheekP/MChinP)
    Pwilcox_CheekP_ChinP<-my.wilcox.p.value(CheekP, ChinP, na.rm=TRUE,paired = FALSE, alternative = "two.sided")

    if(percentZeroNA<0.75){
      pvalues <- c(Pwilcox_Cheek_Chin,  KP_EthinicityDisease,  Pwilcox_CheekChineseH_CheekChineseP, Pwilcox_ChinChineseH_ChinChineseP)
      formatted_pvalues <- sprintf("%.2f", pvalues)
      
      df2 <- data_summary(mydata, varname="gene", 
                          groupnames=c("DiseaseEthBodySite", "DiseaseBodySite"))
      p <- ggplot(df2, aes(x = gene, y = DiseaseEthBodySite, group = DiseaseBodySite, color = DiseaseBodySite)) +
        geom_point(size = 5) +
        geom_errorbar(aes(xmin = gene - se, xmax = gene + se), width = .2,
                      position = position_dodge(0.05)) +
        labs(x = genename, y = "") +
        theme_classic() +
        theme(legend.position = "none") +
        ggtitle(paste(genename,"Psite:", formatted_pvalues[1],"Peth:",formatted_pvalues[2],"PcheekDH:",formatted_pvalues[3],"PChinDH:",formatted_pvalues[4] )) +
        theme(plot.margin = margin(0.5, 0.5, 0.5, 1.5, "cm"),plot.title = element_text(size = 6))  # Adjust the margins as needed
      
      ggsave(paste0(genename, ".png"), p, width = 6, height = 4)
    }

    print(paste(genename,genus2,percentZeroNA,KP_Disease,KP_DiseaseChinese,Pwilcox_PostAcne_Healthy,Pwilcox_ChinChineseH_ChinChineseP,Pwilcox_CheekChineseH_CheekChineseP,
                Pwilcox_CheekChineseH_CheekAAP,Pwilcox_CheekChineseH_CheekCaucasianP,Pwilcox_ChinChineseH_ChinAAP,Pwilcox_ChinChineseH_ChinCaucasianP, 
                Pwilcox_ChineseH_AAP,Pwilcox_ChineseH_CaucasianP,Pwilcox_ChineseH_ChineseP,
                truefc_PostAcne_Healthy,truefc_ChineseP_ChineseH,turefc_ChinChineseP_ChinChineseH,truefc_CheekChineseP_CheekChineseH, 
                MChineseH,MChineseP,MHealthy,MPostAcne,MChinChineseH,MChinChineseP,MChinCaucasianP,MChinAAP,MCheekChineseH,MCheekChineseP,MCheekCaucasianP,MCheekAAP,
                KP_Ethinicity,KP_DiseaseEthinicity,KP_EthinicityDisease, 
                Pwilcox_CheekChineseP_CheekAAP,Pwilcox_CheekChineseP_CheekCaucasianP,Pwilcox_CheekAAP_CheekCaucasianP,
                Pwilcox_ChinChineseP_ChinAAP,Pwilcox_ChinChineseP_ChinCaucasianP,Pwilcox_ChinAAP_ChinCaucasianP,
                Pwilcox_AA_Caucasian, Pwilcox_AA_Chinese, Pwilcox_Caucasian_Chinese,Pwilcox_AAP_CaucasianP,Pwilcox_CheekAAP_CheekCaucasianP,Pwilcox_ChineseP_AAP, Pwilcox_ChineseP_CaucasianP,
                MChinese, MAA, MCaucasian,
                KP_BodySite,KP_BodySiteChinese,Pwilcox_Cheek_Chin,Pwilcox_CheekP_ChinP,Pwilcox_CheekChinese_ChinChinese,
                Pwilcox_ChinCaucasianP_CheekCaucasianP,Pwilcox_ChinAAp_CheekAAP, Pwilcox_ChinChineseP_CheekChineseP, Pwilcox_ChinChineseH_CheekChineseH,
                MChin,MCheek,MChinChinese,MCheekChinese,MChinP,MCheekP,
                truefc_Cheek_Chin,truefc_CheekChinese_ChinChinese,truefc_CheekP_ChinP,
                KP_Sid,KP_BodySiteEthinicityDisease,KP_DiseaseBodysiteChinese,KP_BodySiteDisease,KP_BodySiteEthinicity, KP_DiseaseBodysite,KP_DiseaseBodysiteEthinicity,
                sep=","))

    
    
    
}
