tab=read.table("res_k.tab",sep=';',h=T)
tab
col=c("skyblue4","salmon2")

lim=0.2
col=c("salmon3","darkseagreen4","skyblue4","pink4")
bitmap("roc_human_k.png",type="png256",width=8,height=6.5,res=400)

par(mar=c(4, 4, 0.5, 0.5) + 0.1,cex=2.5)
plot(c(11,43),c(0,100),type="n",xlab="k",ylab=" %",main="prec/recall depending on k",bty='l',axes=F)
axis(side = 1, at=c(13,17,21,25,29,33,37,41), labels=c(13,17,21,25,29,33,37,41), lwd = 2,tcl=0.5,tck=0.02)
axis(side = 2, at=c(10,50,70,80,100), labels=c(10,50,70,80,100), lwd = 2,tcl=0.5,tck=0.02)

lines(tab$rec~tab$k,col=col[1],pch=5,type="l",lty=1,lwd=3)
lines(tab$prec~tab$k,col=col[3],pch=5,type="l",lty=1,lwd=3)
legend(21,18,legend=c("recall", "Precision" ),col=c(col[1],col[3]),lty=0,lwd=3,pch=15,bty='n',yjust=0,pt.cex=1.5)

d=dev.off()


tab=read.table("res_c.tab",h=T)

bitmap("roc_human_c.png",type="png256",width=8,height=6.5,res=400)

par(mar=c(4, 4, 0.5, 0.5) + 0.1,cex=2.5)
plot(c(1,25),c(0,100),type="n",xlab="c",ylab=" %",main="prec/recall depending on k",bty='l',axes=F)
axis(side = 1, at=c(2,4,6,8,10,12,14,16,18,20,22,24), labels=c(2,4,6,8,10,12,14,16,18,20,22,24), lwd = 2,tcl=0.5,tck=0.02)
axis(side = 2, at=c(10,50,70,80,100), labels=c(10,50,70,80,100), lwd = 2,tcl=0.5,tck=0.02)

lines(tab$rec~tab$c,col=col[1],pch=5,type="l",lty=1,lwd=3)
lines(tab$prec~tab$c,col=col[3],pch=5,type="l",lty=1,lwd=3)
legend(12,18,legend=c("recall", "Precision" ),col=c(col[1],col[3]),lty=0,lwd=3,pch=15,bty='n',yjust=0,pt.cex=1.5)

d=dev.off()
