perf=read.table("data_perf.tab",h=T)
perf
## 2 barplot :	
#barplot(perf$memory[perf$type=="bact"],names=perf$name[perf$type=="bact"],col="lightblue")
#library(plotrix)

doubleBarPlot=function(tab){

  tab=tab[order(-tab$memory),]

  # tab$time=tab$time/60
  # my_col=c("slateblue4","aquamarine4","chocolate")
  col=c("gold3","indianred3")
  max1=max(tab$memory)
  max2=max(tab$time)
  n=nrow(tab)

  timeNorm=tab$time*max1/max2
  
  x1=(1:n)-0.35
  x2=(1:n)-0.02
  x3=(1:n)+0.02
  x4=(1:n)+0.35
  
  par(mar=c(4.1, 4, 4, 4.2) + 0.1,cex=2.5,xpd=T,mgp = c(1.3, 0.3, 0))

  plot(c(0.5,n+0.5),c(0,max1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty='l')
  rect(x1,0,x2,tab$memory,col=col[2])
  rect(x3,0,x4,timeNorm,col=col[1])

  axis(1,at=1:n,label=F,tick=F)
  text(1:n, par("usr")[3] - 3, srt = 45, adj = 1,
          labels = tab$name, xpd = TRUE)


  mem.tick=pretty(c(0,tab$memory))
  axis(2,at=mem.tick,label=mem.tick,tcl=0.5,tck=0.02,lwd=2)
  time.tick=pretty(c(0,tab$time))
  axis(4,at=time.tick*max1/max2,label=time.tick,tcl=0.5,tck=0.02,lwd=2)
  mtext("Time (m)", side=4, line=3, cex.lab=1,las=0, cex=2.5,padj=-1)
  mtext("Memory (GB)", side=2, line=3, cex.lab=1,las=0, cex=2.5, padj=1)
  legend(n+0.6,max1*1.55,legend=c("Memory","Time"),fill=c(col[2],col[1]),xjust=1,y.intersp=1.5,bty='n')
}

png(file="perf_human_article.png", width = 1024, height = 850, units = "px", pointsize = 12)
# bitmap("perf_human_article.png",type="png256",width=8,height=7,res=800)
#doubleBarPlot(perf[perf$type=="human",])
tab=perf[perf$type=="human",]
tab=tab[order(-tab$memory),]

  tab$time=tab$time/60
  col=c("gold3","indianred3")
  # col=c("skyblue4","salmon2")
  max1=max(tab$memory)
  max2=max(tab$time)
  n=nrow(tab)

  timeNorm=tab$time*max1/max2
  
  # width rectangles
  x1=(1:n)-0.35
  x2=(1:n)-0.02
  x3=(1:n)+0.02
  x4=(1:n)+0.35
  
  par(mar=c(3.8, 4, 2.5, 4.2) + 0.1,cex=2.5,xpd=T,mgp = c(1.3, 0.3, 0))

  plot(c(0.5,n+0.5),c(0,max1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty='l')
  rect(x1,0,x2,tab$memory,col=col[2])
  rect(x3,0,x4,timeNorm,col=col[1])

  axis(1,at=1:n,label=F,tick=F)
  text(1:n, par("usr")[3] - 6, srt = 35, adj = 1,
          labels = tab$name, xpd = TRUE)


  #mem.tick=pretty(c(0,tab$memory))
mem.tick=c(0,50,100)
  axis(2,at=mem.tick,label=mem.tick,tcl=0.5,tck=0.02,lwd=2)
  #time.tick=pretty(c(0,tab$time))
time.tick=c(0,250,500,750,1000)
  axis(4,at=time.tick*max1/max2,label=time.tick,tcl=0.5,tck=0.02,lwd=2)
  mtext("Time (m)", side=4, line=3, cex.lab=1,las=0, cex=2.5,padj=-1)
  mtext("Memory (GB)", side=2, line=3, cex.lab=1,las=0, cex=2.5, padj=1)
  legend(n+0.6,110,legend=c("Memory","Time"),fill=c(col[2],col[1]),xjust=1,yjust=0,y.intersp=1.5,bty='n',horiz=T)
dev.off()

#png(file="perf_human_article.png", width = 1024, height = 850, units = "px", pointsize = 12)
bitmap("perf_human_article2.png",type="png256",width=7,height=6,res=800)
#doubleBarPlot(perf[perf$type=="human",])
tab=perf[perf$type=="human",]
tab=tab[order(-tab$memory),]

  tab$time=tab$time/60
  #col=c("gold3","indianred3")
  col=c("skyblue4","salmon2")
  max1=max(tab$memory)
  max2=max(tab$time)
  n=nrow(tab)

  timeNorm=tab$time*max1/max2
  
  # width rectangles
  x1=(1:n)-0.35
  x2=(1:n)-0.02
  x3=(1:n)+0.02
  x4=(1:n)+0.35
  
  par(mar=c(3, 4, 2.5, 4.2) + 0.1,cex=3,xpd=T,mgp = c(1.3, 0.3, 0))

  plot(c(0.5,n+0.5),c(0,max1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty='l')
  rect(x1,0,x2,tab$memory,col=col[2])
  rect(x3,0,x4,timeNorm,col=col[1])

  axis(1,at=1:n,label=F,tick=F)
  text(1:n, par("usr")[3] - 6, srt = 35, adj = 1,
          labels = c("Cortex","hybrid","discoSnp"), xpd = TRUE)


  #mem.tick=pretty(c(0,tab$memory))
mem.tick=c(0,50,100)
  axis(2,at=mem.tick,label=mem.tick,tcl=0.5,tck=0.02,lwd=2)
  #time.tick=pretty(c(0,tab$time))
time.tick=c(0,250,500,750,1000)
  axis(4,at=time.tick*max1/max2,label=time.tick,tcl=0.5,tck=0.02,lwd=2)
  mtext("Time (m)", side=4, line=3, cex.lab=1,las=0, cex=2.5,padj=-1)
  mtext("Memory (GB)", side=2, line=3, cex.lab=1,las=0, cex=2.5, padj=1)
  legend(n+0.6,110,legend=c("Memory","Time"),fill=c(col[2],col[1]),xjust=1,yjust=0,y.intersp=1.5,bty='n',horiz=T)
dev.off()

png(file="perf_bact_article.png", width = 1024, height = 850, units = "px", pointsize = 12)
# bitmap("perf_bact_article.png",type="png256",width=10,height=11,res=800)
tab=perf[perf$type=="bact",]
tab=tab[c(3,5,1,4,6,2),]
tab

  tab$time=tab$time/60
    col=c("gold3","indianred3")
  # col=c("skyblue4","salmon2")
  max1=max(tab$memory)
  max2=max(tab$time)
  n=nrow(tab)
  

  timeNorm=tab$time*max1/max2
  
  # width rectangles
s=c(1,2,3,4.5,5.5,6.5)
  x1=s-0.35
  x2=s-0.02
  x3=s+0.02
  x4=s+0.35
  
  par(mar=c(5.1, 4, 3, 4.2) + 0.1,cex=2.5,xpd=T,mgp = c(1.3, 0.3, 0))

  plot(c(0.5,tail(s,1)+0.5),c(0,max1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty='l')
  rect(x1,0,x2,tab$memory,col=col[2])
  rect(x3,0,x4,timeNorm,col=col[1])

  axis(1,at=s,label=F,tick=F)
  text(s, par("usr")[3] - 0.5, srt = 35, adj = 1, labels = tab$name, xpd = TRUE)

#segments(x4[2]+0.15,0,x4[2]+0.15,10,col="grey50",lty=2,lwd=2)
segments(3.75,0,3.75,30,col="grey50",lty=2,lwd=2)

  mem.tick=pretty(c(0,tab$memory))
  axis(2,at=mem.tick,label=mem.tick,tcl=0.5,tck=0.02,lwd=2)
  time.tick=pretty(c(0,tab$time))
# time.tick=seq(0,250,by=50)
  axis(4,at=time.tick*max1/max2,label=time.tick,tcl=0.5,tck=0.02,lwd=2)
  mtext("Time (m)", side=4, line=3, cex.lab=1,las=0, cex=2.5,padj=-1)
  mtext("Memory (GB)", side=2, line=3, cex.lab=1,las=0, cex=2.5, padj=1)
 # legend(tail(s,1)+1,49,legend=c("Memory","Time"),fill=c(col[2],col[1]),xjust=1,yjust=0,bty='n',horiz=T)
legend(s[3],20,legend=c("Memory","Time"),col=c(col[2],col[1]),xjust=1,yjust=0,bty='n',horiz=F, pch=c(15,15))
dev.off()
