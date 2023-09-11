args <- commandArgs(trailingOnly = TRUE)
#print(args)
filename <- args[1]

#library(gridExtra)
#pdf("data_output.pdf", height=11, width=8.5)
#grid.table(mtcars)
#dev.off()

#filename="UADeepGeneRA.stat.xls"
A<-read.table(filename,header=TRUE, sep="\t",stringsAsFactors=F)
d=dim(A)
cols=c(11,12,13,14,18,19,27,28,29,30,31,32,33,39,40)
Ulen=length(cols)
df<-data.frame(matrix(ncol=Ulen, nrow=d[1]))
colnames(df)=paste0("fdr.", colnames(A)[cols])
rownames(df)=A[,1]
for(i in 1:Ulen){
 colN=cols[i]
 p=A[,colN]
 df[,i]=p.adjust(p, method="fdr")
}
write.table(df, file=paste0(filename, ".adjustp"), sep = "\t",eol = "\n", na = "NA", row.names = TRUE,col.names = TRUE)
