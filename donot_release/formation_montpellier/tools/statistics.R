#!/usr/bin/Rscript

#usage :  ./statistics.R  resu_gassst.tab 15000
# or      RScript statistics.R  resu_gassst.tab 15000
args=commandArgs(TRUE)
if (length(args)<2) {
cat("Usage:\n")
cat("./statistics.R  disco.cvs fp.txt [figures.png]\n")
cat("2 obligatory arguments :\n   - csv dico file,\n   - fp list file\n")
cat("optional arguments :\n   - output figure file\n")
quit()
}


fic=args[1]
fp.fic=args[2]

output="0"
if(length(args)>2){
  output=args[3]
}

mycols=c("skyblue2","salmon3")

lines=readLines(fic)
tmp=strsplit(lines,"\\|")
id=as.numeric(unlist(lapply(tmp,function(x) {strsplit(x[1],"_")[[1]][4]}),use.names=F))
rank=as.numeric(unlist(lapply(tmp,function(x) {strsplit(x[length(x)],"_")[[1]][2]}),use.names=F))
type=unlist(lapply(tmp,function(x) {sub(">","",strsplit(x[1],"_")[[1]][1])}),use.names=F)
type[grep("P_2",lines,fixed=T)]="CLOSE"

ntype=length(unique(type))
type.names=sort(unique(type))

fp=read.table(fp.fic)$V3
valid=ifelse(is.element(id,fp),"FP","TP")

cat("Statistics on discoSNP++ results, counts :\n")
cat("------------------------------------------\n")
N=length(lines)
nfp=length(fp)
cat(paste("Nb predictions = ",N,"\n",sep=""))
cat(paste("Nb TP          = ",N-nfp," (",round((N-nfp)*100/N)," %)\n",sep=""))
cat("\n")
validation=valid
validation[validation=="TP"]=" TP" ## to be in the good order (TP before FP)
print(table(type,validation))
cat("------------------------------------------\n")
cat("\n")

d=aggregate(rank,by=list(v=valid,type=type),'median')
cat("Median rank of variants\n")
cat("-----------------------\n")
t=data.frame(TP=d$x[d$v=="TP"],FP=d$x[d$v=="FP"])
rownames(t)=unique(d$type)
print(t)
cat("-----------------------\n")
cat("\n")

if(length(grep("unitig",tmp[[1]]))>0){

## min 
unitig=as.numeric(unlist(lapply(tmp,function(x) {a=strsplit(grep("unitig",x,value=T),"_");n=length(a[[1]]);return(min(as.numeric(a[[1]][n]),as.numeric(a[[2]][n])))}),use.names=F))
d=aggregate(unitig,by=list(v=valid,type=type),'median')
cat("Median distance to the first polymorphism or branching (min unitig length in bp)\n")
cat("------------------------------------------------------------------------------- \n")
t=data.frame(TP=d$x[d$v=="TP"],FP=d$x[d$v=="FP"])
rownames(t)=unique(d$type)
print(t)
cat("------------------------------------------------------------------------------- \n")
cat("\n")

#sum
contig=as.numeric(unlist(lapply(tmp,function(x) {a=strsplit(grep("contig",x,value=T),"_");n=length(a[[1]]);return(sum(as.numeric(a[[1]][n]),as.numeric(a[[2]][n])))}),use.names=F))
d=aggregate(contig,by=list(v=valid,type=type),'median')
cat("Median size of the contig including the variants (in bp)\n")
cat("--------------------------------------------------------\n")
t=data.frame(TP=d$x[d$v=="TP"],FP=d$x[d$v=="FP"])
rownames(t)=unique(d$type)
print(t)
cat("--------------------------------------------------------\n")


## Plotting :
if(output!="0"){
  bitmap(output,type="png256",width=12,height=7,res=400)
  par(mfrow=c(2,2),cex=1.2)
  boxplot(rank~validation+type,ylab="rank",col=mycols,main="Phi coefficient (rank)",names=rep(c("TP","FP"),ntype))
  yrange=par("usr")[4]-par("usr")[3]
  text((1:ntype)*2-0.5,par("usr")[3]-yrange/3,labels=type.names,xpd=T)

  boxplot(unitig~validation+type,log="y",ylab="unitig length (bp)",col=mycols,main="Distance to first branching",names=rep(c("TP","FP"),ntype))
  yrange=par("usr")[4]-par("usr")[3]
  text((1:ntype)*2-0.5,exp(par("usr")[3]-yrange/1.3),labels=type.names,xpd=T)
  
  boxplot(contig~validation+type,log="y",ylab="contig length (pb)",col=mycols,main="Contig size including the variant",names=rep(c("TP","FP"),ntype))
  yrange=par("usr")[4]-par("usr")[3]
  text((1:ntype)*2-0.5,exp(par("usr")[3]-yrange/1.3),labels=type.names,xpd=T)
  
  d=dev.off()
  
  cat("\n")
  cat(paste("Figures in file",output,"\n"))
}
}else{
  bitmap(output,type="png256",width=12,height=7,res=400)
  par(cex=1.4)
  boxplot(rank~validation+type,ylab="rank",col=mycols,main="Phi coefficient (rank)",names=rep(c("TP","FP"),ntype))
  yrange=par("usr")[4]-par("usr")[3]
  text((1:ntype)*2-0.5,par("usr")[3]-yrange/8,labels=type.names,xpd=T)
}
