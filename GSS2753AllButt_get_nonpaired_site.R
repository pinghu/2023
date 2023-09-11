##########################
####Truefc staill have problem
####################################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
#filename="metaphlan.relab10.7.txt"
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.kruskal.p.value <- function(...) {
  obj<-try(kruskal.test(...), silent=TRUE)
  if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}



truefc<-function(VVV){
  #print(VVV)
  if (is.finite(VVV )){
	  XXX=VVV
	  if(VVV==0){
	      XXX=NA
   	  }else if(VVV<1){
	      XXX=-1/VVV
    	}
	  return(XXX)
  }else{
    return("NA")
  }
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=A[1:d[1], 2:d[2]]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ

Cname=colnames(A)[2:d[2]]
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
result=paste("Name", "ShortName","percentZeroNA","mean_all", 
             "mean_NoRash",  "mean_Rash","mean_MildRash_0.5_1",
             "mean_HighRash_1.5_2", "TFC_Rash_NoRash","TFC_HighRash_NoRash",
             "KP_RashCat","KP_ScoreCat",  
             "Pwilcox_Rash_NoRash", "Pwilcox_HighRash_NoRash", 
             "mean_LowPH","mean_HighPH",   "TFC_HighPH_LowPH",
             "KP_PHCat",  "Pwilcox_HighPH_LowPH", 
             "mean_NoRashLowPH","mean_NoRashHighPH",   
             "mean_RashLowPH","mean_RashHighPH", 
             "mean_MildRashLowPH","mean_MildRashHighPH", 
             "mean_HighRashHighPH", 
             "KP_RashPH","KP_ScorePH",  
             "Pwilcox_HighRashHighPH_NoRashHighPH", 
             "Pwilcox_RashHighPH_NoRashHighPH", 
             "Pwilcox_RashLowPH_NoRashLowPH", 
             "Pwilcox_RashHighPH_RashLowPH", 
             "Pwilcox_NoRashHighPH_NoRashLowPH", 
             "TFC_HighRashHighPH_NoRashHighPH",
             "TFC_RashHighPH_NoRashHighPH", 
             "TFC_RashLowPH_NoRashLowPH", 
             "TFC_RashHighPH_RashLowPH", 
             "TFC_NoRashHighPH_NoRashLowPH", 
             "KP_Country", "KP_Sex", 
             sep=",")
print(result)
for (i in 1:d[1]){
  
  genename=A[i,1]
  gene<-as.numeric(C[i,])
  loggene<-log(gene)
  splitG<-strsplit(as.character(genename), "[.|]")
  LLL=length(splitG[[1]])
  genus=splitG[[1]][LLL]
  
  mydata=data.frame(gene,loggene,RashCat, RashPH,PHCat, RashScore, RashScoreN, SkinPH, ScoreCat, ScorePH,  Site, Country,
                    Sex, Age,SID, TEWL)
  #https://github.com/kassambara/rstatix
  
    KP_RashCat=my.kruskal.p.value(gene ~ RashCat, data = mydata, na.action=na.omit)
    KP_PHCat=my.kruskal.p.value(gene ~ PHCat, data = mydata, na.action=na.omit)
    KP_ScoreCat=my.kruskal.p.value(gene ~ ScoreCat, data = mydata, na.action=na.omit)
    KP_RashPH=my.kruskal.p.value(gene ~ RashPH, data = mydata, na.action=na.omit)
    KP_ScorePH=my.kruskal.p.value(gene ~ ScorePH, data = mydata, na.action=na.omit)
    KP_Country = my.kruskal.p.value(gene ~ Country, data = mydata, na.action=na.omit)
    KP_Sex = my.kruskal.p.value(gene ~ Sex, data = mydata, na.action=na.omit)
   
    #print(paste0("KP_PHCat=", KP_PHCat, "; KP_Country=", KP_Country, "; KP_RashCat=", KP_RashCat,
     #            "; KP_ScoreCat=", KP_ScoreCat, "; KP_RashPH=", KP_RashPH, "; KP_ScorePH=", KP_ScorePH, "; KP_Sex=", KP_Sex
    #))
  
    Rash<-gene[ScoreCat=="Rash"]; mRash=mean(Rash);
    NoRash<-gene[ScoreCat=="NoRash"]; mNoRash=mean(NoRash);
    LowPH<-gene[PHCat=="LowPH"]; mLowPH=mean(LowPH);
    HighPH <- gene[PHCat =="HighPH"]; mHighPH=mean(HighPH);
    MildRash<-gene[RashCat=="MildRash(0.5,1)"]; mMildRash=mean(MildRash);
    HighRash<-gene[RashCat=="HighRash(>=1.5)"]; mHighRash=mean(HighRash);
    
    Pwilcox_Rash_NoRash=my.wilcox.p.value(Rash, NoRash, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_Rash_NoRash=truefc(mRash/mNoRash)
    Pwilcox_HighPH_LowPH=my.wilcox.p.value(HighPH, LowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_HighPH_LowPH=truefc(mHighPH/mLowPH)
    Pwilcox_HighRash_NoRash=my.wilcox.p.value(HighRash,NoRash, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_HighRash_NoRash=truefc(mHighRash/mNoRash)

   
    NoRashLowPH<-gene[ScorePH=="NoRash.LowPH"]; mNoRashLowPH=mean(NoRashLowPH);
    NoRashHighPH<-gene[ScorePH=="NoRash.HighPH"]; mNoRashHighPH=mean(NoRashHighPH);
    RashLowPH<-gene[ScorePH=="Rash.LowPH"]; mRashLowPH=mean(RashLowPH);
    RashHighPH<-gene[ScorePH=="Rash.HighPH"]; mRashHighPH=mean(RashHighPH);
    HighRashHighPH<-gene[RashPH=="HighRash(>=1.5).HighPH"]; mHighRashHighPH=mean(HighRashHighPH);
    MildRashHighPH<-gene[RashPH=="MildRash(0.5,1).HighPH"]; mMildRashHighPH=mean(MildRashHighPH);
    MildRashLowPH<-gene[RashPH=="MildRash(0.5,1).LowPH"]; mMildRashLowPH=mean(MildRashLowPH);
    Pwilcox_HighRashHighPH_NoRashHighPH=my.wilcox.p.value(HighRashHighPH, NoRashHighPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    Pwilcox_RashHighPH_NoRashHighPH=my.wilcox.p.value(RashHighPH, NoRashHighPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    Pwilcox_RashLowPH_NoRashLowPH=my.wilcox.p.value(RashLowPH, NoRashLowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    Pwilcox_RashHighPH_RashLowPH=my.wilcox.p.value(RashHighPH, RashLowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    Pwilcox_NoRashHighPH_NoRashLowPH=my.wilcox.p.value(NoRashHighPH, NoRashLowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_HighRashHighPH_NoRashHighPH = truefc(mHighRashHighPH/mNoRashHighPH)
    TFC_RashHighPH_NoRashHighPH = truefc(mRashHighPH/mNoRashHighPH)
    TFC_RashLowPH_NoRashLowPH = truefc(mRashLowPH/mNoRashLowPH)
    TFC_RashHighPH_RashLowPH = truefc(mRashHighPH/mRashLowPH)
    TFC_NoRashHighPH_NoRashLowPH = truefc(mNoRashHighPH/mNoRashLowPH)
    percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
    
    result0=paste(genename, genus,"mean_all", mean(gene),
                 "mean_NoRash", mNoRash, "mean_Rash",mRash,"mean_MildRash", mMildRash, 
                 "mean_HighRash", mHighRash, "TFC_Rash_NoRash",TFC_Rash_NoRash,"TFC_HighRash_NoRash",TFC_HighRash_NoRash,
                 "KP_RashCat",KP_RashCat,"KP_ScoreCat",  KP_ScoreCat, 
                 "Pwilcox_Rash_NoRash", Pwilcox_Rash_NoRash, "Pwilcox_HighRash_NoRash", Pwilcox_HighRash_NoRash,
                 "mean_LowPH",mLowPH,"mean_HighPH", mHighPH,  "TFC_HighPH_LowPH",TFC_HighPH_LowPH,
                 "KP_PHCat",  KP_PHCat,"Pwilcox_HighPH_LowPH", Pwilcox_HighPH_LowPH,
                 "mean_NoRashLowPH",mNoRashLowPH,"mean_NoRashHighPH", mNoRashHighPH,  
                 "mean_RashLowPH",mRashLowPH,"mean_RashHighPH", mRashHighPH,
                 "mean_MildRashLowPH",mMildRashLowPH,"mean_MildRashHighPH", mMildRashHighPH,
                 "mean_HighRashHighPH", mHighRashHighPH,
                 "KP_RashPH",KP_RashPH,"KP_ScorePH",  KP_ScorePH, 
                 "Pwilcox_HighRashHighPH_NoRashHighPH", Pwilcox_HighRashHighPH_NoRashHighPH, 
                 "Pwilcox_RashHighPH_NoRashHighPH", Pwilcox_RashHighPH_NoRashHighPH,
                 "Pwilcox_RashLowPH_NoRashLowPH", Pwilcox_RashLowPH_NoRashLowPH,
                 "Pwilcox_RashHighPH_RashLowPH", Pwilcox_RashHighPH_RashLowPH,
                 "Pwilcox_NoRashHighPH_NoRashLowPH", Pwilcox_NoRashHighPH_NoRashLowPH,
                 "TFC_HighRashHighPH_NoRashHighPH",TFC_HighRashHighPH_NoRashHighPH,
                 "TFC_RashHighPH_NoRashHighPH", TFC_RashHighPH_NoRashHighPH,
                 "TFC_RashLowPH_NoRashLowPH", TFC_RashLowPH_NoRashLowPH,
                 "TFC_RashHighPH_RashLowPH", TFC_RashHighPH_RashLowPH,
                 "TFC_NoRashHighPH_NoRashLowPH", TFC_NoRashHighPH_NoRashLowPH,
                 "KP_Country",  KP_Country,"KP_Sex",  KP_Sex,
                 sep=",") 
    result=paste(genename, genus,percentZeroNA, mean(gene),
                  mNoRash, mRash, mMildRash, 
                 mHighRash, TFC_Rash_NoRash,TFC_HighRash_NoRash,
                 KP_RashCat, KP_ScoreCat, 
                  Pwilcox_Rash_NoRash,  Pwilcox_HighRash_NoRash,
               mLowPH,mHighPH,  TFC_HighPH_LowPH,
                 KP_PHCat, Pwilcox_HighPH_LowPH,
                 mNoRashLowPH, mNoRashHighPH,  
                mRashLowPH, mRashHighPH,
                mMildRashLowPH, mMildRashHighPH,
              mHighRashHighPH,
               KP_RashPH, KP_ScorePH, 
                 Pwilcox_HighRashHighPH_NoRashHighPH, 
                  Pwilcox_RashHighPH_NoRashHighPH,
                  Pwilcox_RashLowPH_NoRashLowPH,
                 Pwilcox_RashHighPH_RashLowPH,
               Pwilcox_NoRashHighPH_NoRashLowPH,
                 TFC_HighRashHighPH_NoRashHighPH,
                 TFC_RashHighPH_NoRashHighPH,
                TFC_RashLowPH_NoRashLowPH,
                 TFC_RashHighPH_RashLowPH,
                 TFC_NoRashHighPH_NoRashLowPH,
                KP_Country, KP_Sex,
                 sep=",") 
    print(result)
    minP=min(KP_RashCat, KP_ScoreCat,Pwilcox_Rash_NoRash,  Pwilcox_HighRash_NoRash, KP_PHCat, Pwilcox_HighPH_LowPH, KP_RashPH, KP_ScorePH, 
                 Pwilcox_HighRashHighPH_NoRashHighPH, 
                  Pwilcox_RashHighPH_NoRashHighPH,
                  Pwilcox_RashLowPH_NoRashLowPH,
                 Pwilcox_RashHighPH_RashLowPH,
               Pwilcox_NoRashHighPH_NoRashLowPH)
	       
    if(percentZeroNA<0.75){
       if(minP<0.05){
      png(filename=paste0(genus,".PHCat.png"),width=600, height = 380)
      stat.test <- mydata %>%
        wilcox_test(gene ~ PHCat) %>%
        mutate(y.position = max(gene))
      p0 <- ggboxplot(mydata, x = "PHCat", y = "gene", color = "PHCat") + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
      print(p0)
      dev.off()
      
      png(filename=paste0(genus,".ScoreCat.png"),width=600, height = 380)
      stat.test <- mydata %>%
        wilcox_test(gene ~ ScoreCat) %>%
        mutate(y.position = max(gene))
      p0 <- ggboxplot(mydata, x = "ScoreCat", y = "gene", color = "ScoreCat") + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
      print(p0)
      dev.off()     
      
      png(filename=paste0(genus,".ScorePH.png"),width=600, height = 380)
      stat.test <- mydata %>%
        wilcox_test(gene ~ ScorePH) %>%
        mutate(y.position = max(gene))
      p0 <- ggboxplot(mydata, x = "ScorePH", y = "gene", color = "ScorePH") + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
      print(p0)
      dev.off()  
      
      png(filename=paste0(genus,".RashPH.png"),width=800, height = 580)
      stat.test <- mydata %>%
        wilcox_test(gene ~ RashPH) %>%
        mutate(y.position = max(gene))
      p0 <- ggboxplot(mydata, x = "RashPH", y = "gene", color = "RashPH") + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)
      print(p0)
      dev.off()      
      
      png(filename=paste0(genus,"Site.png"), width=480, height = 480)
      p1 <- ggboxplot(mydata, x = "Site", y = "gene", color = "Site")+ stat_compare_means()
      print(p1)   
      dev.off()
      
      png(filename=paste0(genus,"Country.png"), width=480, height = 480)
      p1 <- ggboxplot(mydata, x = "Country", y = "gene", color = "Country")+ stat_compare_means()
      print(p1)   
      dev.off()
      
      png(filename=paste0(genus,"Sex.png"), width=480, height = 480)
      p1 <- ggboxplot(mydata, x = "Sex", y = "gene", color = "Sex")+ stat_compare_means()
      print(p1)   
      dev.off()
      }
    }
     
}
