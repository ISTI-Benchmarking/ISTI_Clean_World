moo=scan("HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess015CLS_FEB2015.txt",what="list",sep="\n")
st1=strsplit(moo[1],split=" ")
st1<-st1[[1]][st1[[1]] != ""]
st2=strsplit(moo[2],split=" ")
st2<-st2[[1]][st2[[1]] != ""]
st3=strsplit(moo[3],split=" ")
st3<-st3[[1]][st3[[1]] != ""]
st4=strsplit(moo[4],split=" ")
st4<-st4[[1]][st4[[1]] != ""]
st5=strsplit(moo[5],split=" ")
st5<-st5[[1]][st5[[1]] != ""]
st6=strsplit(moo[6],split=" ")
st6<-st6[[1]][st6[[1]] != ""]
st7=strsplit(moo[7],split=" ")
st7<-st7[[1]][st7[[1]] != ""]
st8=strsplit(moo[8],split=" ")
st8<-st8[[1]][st8[[1]] != ""]
st9=strsplit(moo[9],split=" ")
st9<-st9[[1]][st9[[1]] != ""]
st10=strsplit(moo[10],split=" ")
st10<-st10[[1]][st10[[1]] != ""]

days<-replicate(1908,1)
mons<-replicate(1908/12,seq(12))
years<-replicate(12,seq(1908/12)+1859)
mydates<-paste(days,"/",array(mons,c(1,1908)),"/",array(t(years),c(1,1908)),sep="")
datees<-as.Date(mydates,format="%d/%m/%Y")

setEPS()
postscript("../../IMAGES/HadGEM2-ES_WigglyLoess015_FEB2015.eps")
plot(datees,st1[2:length(st1)],ylab='temperature anomaly (degrees C)',main='HadGEM2-ES Loess 0.15
(Wiggly)',ylim=c(-2,2),cex=0.5)
points(datees,st2[2:length(st2)],col='darkgoldenrod',cex=0.5)
points(datees,st3[2:length(st3)],col='darkgoldenrod1',cex=0.5)
points(datees,st4[2:length(st4)],col='darkgoldenrod4',cex=0.5)

points(datees,st5[2:length(st5)],col='cornflowerblue',cex=0.5)
points(datees,st6[2:length(st6)],col='blue',cex=0.5)
points(datees,st7[2:length(st7)],col='cyan',cex=0.5)
points(datees,st8[2:length(st8)],col='darkblue',cex=0.5)

points(datees,st9[2:length(st9)],col='red',cex=0.5)
points(datees,st10[2:length(st10)],col='brown4',cex=0.5)
dev.off()

moo=scan("HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess021CLS_FEB2015.txt",what="list",sep="\n")
st1=strsplit(moo[1],split=" ")
st1<-st1[[1]][st1[[1]] != ""]
st2=strsplit(moo[2],split=" ")
st2<-st2[[1]][st2[[1]] != ""]
st3=strsplit(moo[3],split=" ")
st3<-st3[[1]][st3[[1]] != ""]
st4=strsplit(moo[4],split=" ")
st4<-st4[[1]][st4[[1]] != ""]
st5=strsplit(moo[5],split=" ")
st5<-st5[[1]][st5[[1]] != ""]
st6=strsplit(moo[6],split=" ")
st6<-st6[[1]][st6[[1]] != ""]
st7=strsplit(moo[7],split=" ")
st7<-st7[[1]][st7[[1]] != ""]
st8=strsplit(moo[8],split=" ")
st8<-st8[[1]][st8[[1]] != ""]
st9=strsplit(moo[9],split=" ")
st9<-st9[[1]][st9[[1]] != ""]
st10=strsplit(moo[10],split=" ")
st10<-st10[[1]][st10[[1]] != ""]

days<-replicate(1908,1)
mons<-replicate(1908/12,seq(12))
years<-replicate(12,seq(1908/12)+1859)
mydates<-paste(days,"/",array(mons,c(1,1908)),"/",array(t(years),c(1,1908)),sep="")
datees<-as.Date(mydates,format="%d/%m/%Y")

setEPS()
postscript("../../IMAGES/HadGEM2-ES_WavyLoess021_FEB2015.eps")
plot(datees,st1[2:length(st1)],ylab='temperature anomaly (degrees C)',main='HadGEM2-ES Loess 0.21
(Wavy)',ylim=c(-2,2),cex=0.5)
points(datees,st2[2:length(st2)],col='darkgoldenrod',cex=0.5)
points(datees,st3[2:length(st3)],col='darkgoldenrod1',cex=0.5)
points(datees,st4[2:length(st4)],col='darkgoldenrod4',cex=0.5)

points(datees,st5[2:length(st5)],col='cornflowerblue',cex=0.5)
points(datees,st6[2:length(st6)],col='blue',cex=0.5)
points(datees,st7[2:length(st7)],col='cyan',cex=0.5)
points(datees,st8[2:length(st8)],col='darkblue',cex=0.5)

points(datees,st9[2:length(st9)],col='red',cex=0.5)
points(datees,st10[2:length(st10)],col='brown4',cex=0.5)
dev.off()

moo=scan("HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess031CLS_FEB2015.txt",what="list",sep="\n")
st1=strsplit(moo[1],split=" ")
st1<-st1[[1]][st1[[1]] != ""]
st2=strsplit(moo[2],split=" ")
st2<-st2[[1]][st2[[1]] != ""]
st3=strsplit(moo[3],split=" ")
st3<-st3[[1]][st3[[1]] != ""]
st4=strsplit(moo[4],split=" ")
st4<-st4[[1]][st4[[1]] != ""]
st5=strsplit(moo[5],split=" ")
st5<-st5[[1]][st5[[1]] != ""]
st6=strsplit(moo[6],split=" ")
st6<-st6[[1]][st6[[1]] != ""]
st7=strsplit(moo[7],split=" ")
st7<-st7[[1]][st7[[1]] != ""]
st8=strsplit(moo[8],split=" ")
st8<-st8[[1]][st8[[1]] != ""]
st9=strsplit(moo[9],split=" ")
st9<-st9[[1]][st9[[1]] != ""]
st10=strsplit(moo[10],split=" ")
st10<-st10[[1]][st10[[1]] != ""]

days<-replicate(1908,1)
mons<-replicate(1908/12,seq(12))
years<-replicate(12,seq(1908/12)+1859)
mydates<-paste(days,"/",array(mons,c(1,1908)),"/",array(t(years),c(1,1908)),sep="")
datees<-as.Date(mydates,format="%d/%m/%Y")

setEPS()
postscript("../../IMAGES/HadGEM2-ES_NormalLoess031_FEB2015.eps")
plot(datees,st1[2:length(st1)],ylab='temperature anomaly (degrees C)',main='HadGEM2-ES Loess 0.31
(Normal)',ylim=c(-2,2),cex=0.5)
points(datees,st2[2:length(st2)],col='darkgoldenrod',cex=0.5)
points(datees,st3[2:length(st3)],col='darkgoldenrod1',cex=0.5)
points(datees,st4[2:length(st4)],col='darkgoldenrod4',cex=0.5)

points(datees,st5[2:length(st5)],col='cornflowerblue',cex=0.5)
points(datees,st6[2:length(st6)],col='blue',cex=0.5)
points(datees,st7[2:length(st7)],col='cyan',cex=0.5)
points(datees,st8[2:length(st8)],col='darkblue',cex=0.5)

points(datees,st9[2:length(st9)],col='red',cex=0.5)
points(datees,st10[2:length(st10)],col='brown4',cex=0.5)
dev.off()

moo=scan("HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess04CLS_FEB2015.txt",what="list",sep="\n")
st1=strsplit(moo[1],split=" ")
st1<-st1[[1]][st1[[1]] != ""]
st2=strsplit(moo[2],split=" ")
st2<-st2[[1]][st2[[1]] != ""]
st3=strsplit(moo[3],split=" ")
st3<-st3[[1]][st3[[1]] != ""]
st4=strsplit(moo[4],split=" ")
st4<-st4[[1]][st4[[1]] != ""]
st5=strsplit(moo[5],split=" ")
st5<-st5[[1]][st5[[1]] != ""]
st6=strsplit(moo[6],split=" ")
st6<-st6[[1]][st6[[1]] != ""]
st7=strsplit(moo[7],split=" ")
st7<-st7[[1]][st7[[1]] != ""]
st8=strsplit(moo[8],split=" ")
st8<-st8[[1]][st8[[1]] != ""]
st9=strsplit(moo[9],split=" ")
st9<-st9[[1]][st9[[1]] != ""]
st10=strsplit(moo[10],split=" ")
st10<-st10[[1]][st10[[1]] != ""]

days<-replicate(1908,1)
mons<-replicate(1908/12,seq(12))
years<-replicate(12,seq(1908/12)+1859)
mydates<-paste(days,"/",array(mons,c(1,1908)),"/",array(t(years),c(1,1908)),sep="")
datees<-as.Date(mydates,format="%d/%m/%Y")

setEPS()
postscript("../../IMAGES/HadGEM2-ES_FatLoess04_FEB2015.eps")
plot(datees,st1[2:length(st1)],ylab='temperature anomaly (degrees C)',main='HadGEM2-ES Loess 0.4
(Fat)',ylim=c(-2,2),cex=0.5)
points(datees,st2[2:length(st2)],col='darkgoldenrod',cex=0.5)
points(datees,st3[2:length(st3)],col='darkgoldenrod1',cex=0.5)
points(datees,st4[2:length(st4)],col='darkgoldenrod4',cex=0.5)

points(datees,st5[2:length(st5)],col='cornflowerblue',cex=0.5)
points(datees,st6[2:length(st6)],col='blue',cex=0.5)
points(datees,st7[2:length(st7)],col='cyan',cex=0.5)
points(datees,st8[2:length(st8)],col='darkblue',cex=0.5)

points(datees,st9[2:length(st9)],col='red',cex=0.5)
points(datees,st10[2:length(st10)],col='brown4',cex=0.5)
dev.off()
