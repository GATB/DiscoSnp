rm(tab,t1,t2,tab2,t3,t4,final)

tab=read.table("res_b_0",sep="=")
t1=tab[tab$V1=="precision",]
t2=tab[tab$V1=="recall",]

tab2=read.table("res_b_1",sep="=")

t3=tab2[tab2$V1=="precision",]
head(t3$V2)
t4=tab2[tab2$V1=="recall",]

final=data.frame(precision_b0=t1$V2,recall_b0=t2$V2, precision_b1=t3$V2,recall_b1=t4$V2)


final=data.frame(precision_b1=t3$V2,recall_b1=t4$V2)

head(final)
boxplot(final, main="DiscoDel, coli, noSNP, 50del, insert size [1;10], k=31, D 10 c 5",ylabel="%")

