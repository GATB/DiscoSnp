
final=NULL
nb=2:30
seuils=c("0_0","0_2","0_5")




#    _____  _                                     _ _       
#   |  __ \(_)                                   | | |      
#   | |  | |_ ___  ___ ___    _ __ ___  ___ _   _| | |_ ___ 
#   | |  | | / __|/ __/ _ \  | '__/ _ \/ __| | | | | __/ __|
#   | |__| | \__ \ (_| (_) | | | |  __/\__ \ |_| | | |_\__ \
#   |_____/|_|___/\___\___/  |_|  \___||___/\__,_|_|\__|___/
#                                                           
#                                                           



for (seuil in seuils){
  temp=read.table(paste("res_discosnp++/res_b1_D_10_P_4_c_4_d_1/res_temp_afac_t_",seuil,"_pos_SNP_precision.txt",sep=""))
  a=data.frame(seuil=rep(seuil,length(nb)),nb=nb,disco_precision_SNP=temp[,1])
  a
  
  temp=read.table(paste("res_discosnp++/res_b1_D_10_P_4_c_4_d_1/res_temp_afac_t_",seuil,"_pos_SNP_recall.txt",sep=""))
  a=cbind(a,data.frame(disco_recall_SNP=temp[,1]))
  
  
  temp=read.table(paste("res_discosnp++/res_b1_D_10_P_4_c_4_d_1/res_temp_afac_t_",seuil,"_pos_INDEL_precision.txt",sep=""))
  a=cbind(a,data.frame(disco_precision_INDEL=temp[,1]))
  
  
  temp=read.table(paste("res_discosnp++/res_b1_D_10_P_4_c_4_d_1/res_temp_afac_t_",seuil,"_pos_INDEL_recall.txt",sep=""))
  a=cbind(a,data.frame(disco_recall_INDEL=temp[,1]))
  
  final=rbind(final,a)
}




#     _____           _                                  _ _       
#    / ____|         | |                                | | |      
#   | |     ___  _ __| |_ _____  __  _ __ ___  ___ _   _| | |_ ___ 
#   | |    / _ \| '__| __/ _ \ \/ / | '__/ _ \/ __| | | | | __/ __|
#   | |___| (_) | |  | ||  __/>  <  | | |  __/\__ \ |_| | | |_\__ \
#    \_____\___/|_|   \__\___/_/\_\ |_|  \___||___/\__,_|_|\__|___/
#                                                                  
#                                                                  


seuil="0_0"
temp=read.table(paste("res_cortex_sabre/res_temp_afac_t_",seuil,"_pos_SNP_precision.txt",sep=""))
a=data.frame(cortex_precision_SNP=temp[,1])
a

temp=read.table(paste("res_cortex_sabre/res_temp_afac_t_",seuil,"_pos_SNP_recall.txt",sep=""))
a=cbind(a,data.frame(cortex_recall_SNP=temp[,1]))


temp=read.table(paste("res_cortex_sabre/res_temp_afac_t_",seuil,"_pos_INDEL_precision.txt",sep=""))
a=cbind(a,data.frame(cortex_precision_INDEL=temp[,1]))


temp=read.table(paste("res_cortex_sabre/res_temp_afac_t_",seuil,"_pos_INDEL_recall.txt",sep=""))
a=cbind(a,data.frame(cortex_recall_INDEL=temp[,1]))

final=cbind(final,a)




head(final)
#                                        _   _                          _ _       
#                                       | | | |                        | | |      
#    ___  ___   __ _ _ __     __ _  __ _| |_| | __  _ __ ___  ___ _   _| | |_ ___ 
#   / __|/ _ \ / _` | '_ \   / _` |/ _` | __| |/ / | '__/ _ \/ __| | | | | __/ __|
#   \__ \ (_) | (_| | |_) | | (_| | (_| | |_|   <  | | |  __/\__ \ |_| | | |_\__ \
#   |___/\___/ \__,_| .__/   \__, |\__,_|\__|_|\_\ |_|  \___||___/\__,_|_|\__|___/
#                   | |______ __/ |                                               
#                   |_|______|___/                                                


seuil="0_0"
temp=read.table(paste("res_soap_gatk_genocluster/res_temp_afac_t_",seuil,"_pos_SNP_precision.txt",sep=""))
a=data.frame(soap_gatk_precision_SNP=temp[,1])
a

temp=read.table(paste("res_soap_gatk_genocluster/res_temp_afac_t_",seuil,"_pos_SNP_recall.txt",sep=""))
a=cbind(a,data.frame(soap_gatk_recall_SNP=temp[,1]))


temp=read.table(paste("res_soap_gatk_genocluster/res_temp_afac_t_",seuil,"_pos_INDEL_precision.txt",sep=""))
a=cbind(a,data.frame(soap_gatk_precision_INDEL=temp[,1]))


temp=read.table(paste("res_soap_gatk_genocluster/res_temp_afac_t_",seuil,"_pos_INDEL_recall.txt",sep=""))
a=cbind(a,data.frame(soap_gatk_recall_INDEL=temp[,1]))

final=cbind(final,a)




#          _       _   
#         | |     | |  
#    _ __ | | ___ | |_ 
#   | '_ \| |/ _ \| __|
#   | |_) | | (_) | |_ 
#   | .__/|_|\___/ \__|
#   | |                
#   |_|                

my_col=c("skyblue4","burlywood4","salmon2")

###############
## INDELS
###############
png("n_coli.png",width=825,height=646,pointsize=12,bg="white")
par(mfrow=c(1,2))
my_col=c("skyblue4","burlywood4","salmon2")

plot(c(2,30),c(10,100),type="n",xlab="nb read sets",ylab="Precision and recall (%)",main="INDEL precision/recall")
id_seuil=1 # 0_0 for cortex and soap_gatk
# CORTEX
lines(final$nb[final$seuil==seuils[id_seuil]],final$cortex_precision_INDEL[final$seuil==seuils[id_seuil]],col=my_col[1],pch=5,type="b",lwd=2,cex=1)
lines(final$nb[final$seuil==seuils[id_seuil]],final$cortex_recall_INDEL[final$seuil==seuils[id_seuil]],col=my_col[1],pch=23,bg=my_col[1],type="b",lwd=2,cex=1)
# SOAP GATK
lines(final$nb[final$seuil==seuils[id_seuil]],final$soap_gatk_precision_INDEL[final$seuil==seuils[id_seuil]],col=my_col[2],pch=1,type="b",lwd=2,cex=1)
lines(final$nb[final$seuil==seuils[id_seuil]],final$soap_gatk_recall_INDEL[final$seuil==seuils[id_seuil]],col=my_col[2],pch=16,type="b",lwd=2,cex=1)
# DISCO
lines(final$nb[final$seuil==seuils[id_seuil]],final$disco_precision_INDEL[final$seuil==seuils[id_seuil]],col=my_col[3],pch=2,type="b",lwd=2,cex=1)
lines(final$nb[final$seuil==seuils[id_seuil]],final$disco_recall_INDEL[final$seuil==seuils[id_seuil]],col=my_col[3],pch=17,type="b",lwd=2,cex=1)
# legends
legend(15,60,legend=c("Precision","Recall"),col=c("grey","grey"),lty=c(0,0),lwd=2,pch=c(0,15),bg="white",bty='n',pt.cex=1.5)
legend(20,40,legend=c("Cortex", "Soap+Gatk", "DiscoSnp++"),bg=my_col[1],col=c(my_col[1],my_col[2],my_col[3]),lty=1,lwd=1,pch=c(18,16,17),bty='n', pt.cex=c(2,1.5,1.5))



###############
## SNP 
###############

plot(c(2,30),c(10,100),type="n",xlab="nb read sets",ylab="Precision and recall (%)",main="SNP precision/recall")
id_seuil=1 # 0_0 for cortex and soap_gatk
# CORTEX
lines(final$nb[final$seuil==seuils[id_seuil]],final$cortex_precision_SNP[final$seuil==seuils[id_seuil]],col=my_col[1],pch=5,type="b",lwd=2,cex=1)
lines(final$nb[final$seuil==seuils[id_seuil]],final$cortex_recall_SNP[final$seuil==seuils[id_seuil]],col=my_col[1],pch=23,bg=my_col[1],type="b",lwd=2,cex=1)
# SOAP GATK
lines(final$nb[final$seuil==seuils[id_seuil]],final$soap_gatk_precision_SNP[final$seuil==seuils[id_seuil]],col=my_col[2],pch=1,type="b",lwd=2,cex=1)
lines(final$nb[final$seuil==seuils[id_seuil]],final$soap_gatk_recall_SNP[final$seuil==seuils[id_seuil]],col=my_col[2],pch=16,type="b",lwd=2,cex=1)
# DISCO
lines(final$nb[final$seuil==seuils[id_seuil]],final$disco_precision_SNP[final$seuil==seuils[id_seuil]],col=my_col[3],pch=2,type="b",lwd=2,cex=1)
lines(final$nb[final$seuil==seuils[id_seuil]],final$disco_recall_SNP[final$seuil==seuils[id_seuil]],col=my_col[3],pch=17,type="b",lwd=2,cex=1)
# legends
legend(15,60,legend=c("Precision","Recall"),col=c("grey","grey"),lty=c(0,0),lwd=2,pch=c(0,15),bg="white",bty='n',pt.cex=1.5)
legend(20,40,legend=c("Cortex", "Soap+Gatk", "DiscoSnp++"),bg=my_col[1],col=c(my_col[1],my_col[2],my_col[3]),lty=1,lwd=1,pch=c(18,16,17),bty='n', pt.cex=c(2,1.5,1.5))


