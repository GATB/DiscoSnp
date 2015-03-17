args <- commandArgs(TRUE)
print(args[1])


lim=0.2
col=c("salmon3","darkseagreen4","skyblue4","pink4")

bitmap(args[1],type="png256",width=8,height=6.5,res=400)


#dev.off()
par(mar=c(4, 4, 0.5, 0.5) + 0.1,cex=2.5)
plot(c(0,80),c(30,100),type="n",xlab="Recall (%)",ylab="Precision (%)",main="",bty='l',axes=F)


axis(side = 1, at=c(0,50,70,80,90), labels=c(0,50,70,80,90), lwd = 2,tcl=0.5,tck=0.02)
axis(side = 2, at=c(50,70,80,90,100), labels=c(50,70,80,90,100), lwd = 2,tcl=0.5,tck=0.02)

legend(9.6,40,legend=c(expression(discoSnp++ rank >= 0.2),expression(discoSnp++ rank < 0.2)),col="black",lty=c(1,3),lwd=3, bty='n',yjust=0,text.width=10)



snp=read.table("roc_SNP")
df <- try(read.table("roc_INDEL")) 
if (inherits(df, 'try-error')){  
       print("No indels") 


       legend(52.3,28,legend=c("SNPs"),col=c(col[1]),lty=0,lwd=3,pch=15,bty='n',yjust=0,pt.cex=1.5)

       cool_snp=snp[which(snp$V3 > lim),]
       

       lines(cool_snp$V2~cool_snp$V1,col=col[1],pch=5,type="l",lty=1,lwd=3)
       lines(snp$V2~snp$V1,col=col[1],pch=5,type="l",lty = "dotted",lwd=3)
       } else {
       indel=read.table("roc_INDEL")

       legend(52.3,28,legend=c("SNPs", "Indels" ),col=c(col[1],col[3]),lty=0,lwd=3,pch=15,bty='n',yjust=0,pt.cex=1.5)


cool_snp=snp[which(snp$V3 > lim),]
cool_indel=indel[which(indel$V3 > lim),]

snp=snp[which(snp$V3 <= lim),]
indel=indel[which(indel$V3 <= lim),]

lines(cool_snp$V2~cool_snp$V1,col=col[1],pch=5,type="l",lty=1,lwd=3)
lines(cool_indel$V2~cool_indel$V1,col=col[3],pch=5,type="l",lty=1,lwd=3)
lines(snp$V2~snp$V1,col=col[1],pch=5,type="l",lty = "dotted",lwd=3)
lines(indel$V2~indel$V1,col=col[3],pch=5,type="l",lty = "dotted",lwd=3)
} 

d=dev.off()
