rm(tab,t1,t2,tab2,t3,t4,final)

tab=read.table("res_b_0",sep="=")
t1d=tab[tab$V1=="precision_DEL",]
t2d=tab[tab$V1=="recall_DEL",]

t1s=tab[tab$V1=="precision_SNP",]
t2s=tab[tab$V1=="recall_SNP",]

tab2=read.table("res_b_1",sep="=")

t3d=tab2[tab2$V1=="precision_DEL",]
t4d=tab2[tab2$V1=="recall_DEL",]
t3s=tab2[tab2$V1=="precision_SNP",]
t4s=tab2[tab2$V1=="recall_SNP",]

final=data.frame(p_D_b0=t1d$V2,r_D_b0=t2d$V2, p_S_b0=t1s$V2, r_S_b0=t2s$V2, p_D_b1=t3d$V2,r_D_b1d=t4d$V2, p_S_b1=t3s$V2,r_S_b1s=t4s$V2)


# ONLY B0
#final=data.frame(precision_b0=t1$V2,recall_b0=t2$V2)

# ONLY B1
#final=data.frame(precision_b1=t3$V2,recall_b1=t4$V2)

head(final)
boxplot(final, main="DiscoDel, coli, 1000 SNP, 2x50x, 1% error rate, \n100 insertions, insert size [1;10], k=31, D 10 c 5",xlab = expression(bold("p/r: precision/recall, D/S: Deletion/Snp")))
summary(final)
nrow(final)
