######################
####VEGAN AN√ÅLISIS####
######################
rm(list = ls())
graphics.off()
library(reshape2)
library(vegan)
library(pvclust)
library(ggplot2)
library(car)
??vegan
vegangraphics.off()
rm(list = ls())
data(BCI)
BCI
str(BCI)
View(BCI)
diversity(BCI, index = "simpson")
diversity(BCI, index = "shannon")
?simpson
H<-diversity(BCI)
H
##ÕNDICES DDE DIVERSIDAD Y EQUITATIVIDAD SEG⁄N KAZUYA NAOKI##
setwd("C:/Users/Gabriel Torrico/Documents/BiologÌa/PROYECTO COFFEE/Bases de datos")
matabun<-read.csv("matabun2.csv",sep = ";", dec = ",", header = T)
row.names(matabun)<-as.character(matabun$Parcela)
matabun$Parcela<-NULL
str(matabun)
specnumber(matabun) #riqueza de especies por fila#
specnumber(matabun, groups = "pooled") #riqueza de especies para el juego de datos completo#
S<-specnumber(matabun)
#Ìndices de diversidad/heterogeneidad#
shannon<-diversity(matabun, index = "shannon")
invsimp<-diversity(matabun, index = "invsimpson")
eqshannon<-shannon/log(S) #equitatitividad de shannon#
eqinvsimpson<-invsimp/S #equitatividad de inverso de simpson#
#n˙mero de especies#
Abum<-rowSums(matabun)
#todo como dataframe#
index<-data.frame(S, Abum, shannon, invsimp, eqshannon, eqinvsimpson)
export(index, "Ìndices_FINAL.xlsx", row.names=T)
 # # # # # # # # # # # # # ## # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # preparar la base de datos # # # # ## # # # # # # matriz de abundancias absolutas #
# # # # # # # # # # # # # ## # # # # # # # # # # # # # # matriz de presencia y ausencia #
 # # # # # # # # # # # # # ## # # # # # # # # # # # # #
setwd("C:/Users/Gabriel Torrico/Documents/BiologÌa/PROYECTO COFFEE/Bases de datos")
getwd()
basecero<-read.csv("basecero.csv",sep = ";", dec = ",", header = T) #<- LA BASE IMPORTANTE
head(basecero)
sppabu<-dcast(basecero, Parcela~Nombre.com˙n, value.var = "DAP")
sppabun<-as.data.frame(sppabu)
sppabun[is.na(sppabun)] <- 0
row.names(sppabun)<-sppabun$Parcela ##This makes the row names the same as plot ID, especially useful for plot IDs for ordinations#
sppabun<-sppabun[-1] #dropping first column (we already made it the row names)
sppabun
library(rio)
export(sppabun, "mat.abun.xlsx", row.names=T)
# a l t e r n a t i v o #
sppabu<-dcast(basecero, Parcela~species, value.var = "DAP")
sppabun<-as.data.frame(sppabu)
sppabun[is.na(sppabun)] <- 0
row.names(sppabun)<-sppabun$Parcela ##This makes the row names the same as plot ID, especially useful for plot IDs for ordinations#
sppabun<-sppabun[-1] #dropping first column (we already made it the row names)
sppabun
rowSums(sppabun)

library(rio)
export(sppabun, "mat.abun2.xlsx", row.names=T)
all.equal(basecero$Nombre.com˙n,basecero$common.name)

#Presence absence matrix#
sppabun[sppabun > 0] <- 1 #converts from abundance to P/A
sppabun
export(sppabun, "mat.01.xlsx", row.names=T)
#otra base#
basecont<-read.csv("baseabu.csv", sep = ";", dec = ",", header = T)
basecont
row.names(basecont)<-basecont$Propietario
basecont$Propietario<-NULL
basecont22 <- basecont/rowSums(basecont)
basecont<-prop.table(basecont)
head(basecont22)
rowSums(basecont22)
baseabu<-read.csv("baseabu.csv", sep = ";", dec = ",", header = T)
row.names(baseabu)<-baseabu$Propietario
baseabu$Propietario<-NULL
baseabu
####ÕNDICES DE DIVERSIDAD####
#Ìndice de shannon-wiener#
H<-diversity(sppabun)
H
print(H, digits = 3)
#Õndice de simpson#
s<-diversity(sppabun, index = "simpson")
s
print(s, digits = 2)
#Õndice de equitabilidad de Pielou#
J <- H/log(specnumber(sppabun))
J
print(J, digits = 2)
#n˙mero de especies#
No<-specnumber(sppabun)
Ab<-rowSums(sppabun)
#todo como dataframe#
index<-data.frame(H, s, J, No, Ab)
export(index, "Ìndicesfinal.xlsx", row.names=T)
####curva de acumulaciÛnn####
curve<-specaccum(sppabun)
plot(curve)
plot(curve, ci.type="polygon", ci.col="yellow")
#number of unseen species#
specpool(basecont)
#####ordenaci√≥n####
modd<-decorana(baseabu)
plot(modd)
##### CA ####
ca.basecont<-cca(baseabu)
summary(ca.basecont)
plot(ca.basecont)
biplot(ca.basecont)
#Visualice los objetos que contiene el resultado del CA#
str(ca.basecont)
str(ca.basecont$CA)
#ELABORAR LAS FIGURAS PARA LA PUBLICACI√É¬ìN#
#juego de datos con las puntuaciones de los primeros tres ejes CA1,2,3, de sitios
datos.ca.sitios<-data.frame(ca.basecont$CA$u[,1:3])
#juego de datos con las puntuaciones de los primeros tres ejes CA1,2,3 de especies
datos.ca.especies<-data.frame(ca.basecont$CA$v[,1:3])
###abreviar los nombres de las especies
datos.ca.especies$esp<-row.names(datos.ca.especies)
datos.ca.especies$esp2<-abbreviate(datos.ca.especies$esp, minlength = 6, strict = FALSE)
datos.ca.especies
#####GR√ÅFICO CA EN GGPLOT2####
library(ggplot2)
ggplot()+
  geom_text(data = datos.ca.sitios, aes(CA1, CA2, label=row.names(datos.ca.sitios)),
            col="#FF0000", size=4)+
  geom_text(data = datos.ca.especies, aes(CA1, CA2, label=esp2),
            col="#006699", size=3)+
  theme(axis.title.x = element_text(face = "bold", size = 14, color = "black"),
        axis.title.y = element_text(face = "bold", size = 14),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"))+
  labs(title = "CA de los √°rboles acompa√±antes en las parcelas de estudio",
       x="CA1(21%)",
       y="CA2(17%)")+
  theme_test()
#ordiplot 3d#
ordiplot3d(ca.basecont, type = "h", display = "sites")
?ordiplot

#################
#### I V I ######
#################

##### IVI 1 ####
library(BiodiversityR)
library(tidyverse)
library(readxl)
#import data#
tutorial<-"C:/Users/Gabriel Torrico/Documents/BiologÌa/PROYECTO COFFEE/Bases de datos/base 0000.xlsx"
excel_sheets(tutorial)
data<-read_excel(tutorial, "basecero")
head(data)
?importancevalue
?mutate
tabla.ivi<-data%>%as.data.frame()%>%
  mutate(count=rep(1, each=nrow(.)),
         ba=0.7854*(DAP/100)^2)%>%
  importancevalue(site="Parcela", species="Nombre cientÌfico", 
                                                  count="count", basal="ba", 
                                                  factor="", level="") 
tabla.ivi  
View(tabla.ivi)
tabla.ivi2<-print(tabla.ivi, digits = 2)
export(tabla.ivi, "IVI finales.xlsx", row.names=T)
##### IVI 2 ####
tutor<-"C:/Users/Gabriel Torrico/Documents/BiologÌa/PROYECTO COFFEE/Bases de datos/tutor.xlsx"
excel_sheets(tutor)
data<-read_excel(tutor, "Hoja1")
head(data)
?importancevalue
?mutate
tabla.ivi<-data%>%as.data.frame()%>%
  mutate(count=rep(1, each=nrow(.)),
         ba=0.7854*(dbh/100)^2)%>%
  importancevalue(site="plot", species="common_name", 
                  count="count", basal="ba", 
                  factor="", level="") 
tabla.ivi  
View(tabla.ivi)
tabla.ivi2<-print(tabla.ivi, digits = 2)
##### exportar datos a excel #####
library(rio)
export(tabla.ivi2, "IVI FINAL.xlsx", row.names=T)

write.csv(tabla.ivi, file = "tabla_ivi_final.csv", row.names = T, col.names = T, sep = ",")
basecont

#clusters#
clua <- hclust(baseabu, "average")
plot(clua) 
view(baseabu)
?histogram

########################################
######### G I G A P L O T S  ###########
########################################
#####plot histograma altura####
#√°rea de trabajo#

base<-read.csv("base0.csv", dec = ",", sep = ";")
str(base)
head(base)
summary(base$Ht)
basecero<-read.csv("basecero.csv", dec = ",", sep = ";")
head(basecero)
str(basecero$Nombre.com˙n)
ggplot(basecero, aes(x=reorder(Familia, Familia, function(x)-length(x)))) +
  geom_bar(fill='red') +  labs(x='Familia')+geom_bar(stat = "count") + 
  stat_count(geom = "text", colour = "white", size = 3.5,
             aes(label = ..count..),position=position_stack(vjust=0.5))
###################################################################################
######### G I G A H I S T O G R A M A S - histograma clases de altura #############
###################################################################################
 
library(viridis)
library(hrbrthemes)
hist(base$Ht, main = "Companion trees height vs. frequency", xlab = "Height classes (m)", ylab = "Absolute frequency", freq = F, las=1)
range(base$Ht)
brkvec<-seq(1, 33, by=2)
length(brkvec)
hist<-hist(base$Ht, breaks = brkvec, col="cyan3", border="green",xaxp=c(1,33, 16),main = "Altura de ·rboles acompaÒantes vs. frecuencia", 
     xlab = "CategorÌas de altura (m)", ylab = "Frecuencia absoluta", freq = T, las=1, cex.lab=1.4)
lines(c(min(hist$breaks),hist$mids,max(hist$breaks)),c(0,hist$counts,0),type = "l",col="darkblue")
text(hist$mids, hist$counts, labels = hist$counts)
hist<-hist(base$Ht, breaks = brkvec, col="khaki4", xaxp=c(1,33, 16),
           xlab = "CategorÌas de altura (m)", ylab = "Frecuencia absoluta", freq = T, las=1)
lines(c(min(hist$breaks),hist$mids,max(hist$breaks)),c(0,hist$counts,0),type = "l",col="blue")
text(hist$mids, hist$counts, labels = hist$counts)

midpt=c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)
freq=c(9,107,153,149,138,136,58,31,10,6,2,4,2,0,1,1)
data=rep(midpt, freq);data
brk=c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33)
hhistt<-hist(data, breaks=brk, col=c("yellow","blue"),lwd=2,main = "Altura de ·rboles acompaÒantes vs. frecuencia",
             xlab = "CategorÌas de altura (m)", ylab = "Frecuencia absoluta", freq = T, las=1)
lines(c(min(hist$breaks),hist$mids,max(hist$breaks)),c(0,hist$counts,0),type = "l",col="darkblue")
text(hist$mids, hist$counts, labels = hist$counts)

#########################
#NUEVOS GIGA HISTOGRAMAS#
#########################
head(basecero)
hist(basecero$DAP)
range(basecero$DAP)
range(basecero$Ht)
range(basecero$Abasal)
range(basecero$Hc)
par (mfrow = c(2, 1)) #definir diposicÌon c(filas, columas)#
par(mar=c(5,5,0.4,0.1)+0.1) #definir m·rgenes: base, izquierda, arriba, derecha + yapa numÈrica#
brkvec<-seq(1, 33, by=2)
length(brkvec)
hist1<-hist(basecero$Ht, breaks = brkvec, col="gray90", border="gray70",xaxp=c(1,33, 16),main = "", 
           xlab = "Height categories (m)", ylab = "Absolute frequency", freq = T, las=1, cex.lab=1.15)
lines(c(min(hist1$breaks),hist1$mids,max(hist1$breaks)),c(0,hist1$counts,0),type = "l",col="black")
text(hist1$mids, hist1$counts, labels = hist1$counts, cex=1)
legend("topright", "a", bg="black", text.col="white", adj=0.7)
?lines
#DAP#
range(basecero$DAP)
brkvec<-seq(0, 130, by=5)
length(brkvec)
hist2<-hist(basecero$DAP, breaks = brkvec, col="gray90", border="gray70",xaxp=c(0,130, 26),main = "", 
           xlab = "DBH categories (cm)", ylab = "Absolute frequency", freq = T, las=1, cex.lab=1.15)
lines(c(min(hist2$breaks),hist2$mids,max(hist2$breaks)),c(0,hist2$counts,0),type = "l",col="black")
text(hist2$mids, hist2$counts, labels = hist2$counts, cex = 1)
legend("topright", "b", bg="black", text.col="white", adj=0.7)

midpt=c(3,9,15,21,27,33,39,45,51,57,63,69,75,81,87,93,99,105,111,117,123)
freq=c(11,540,430,217,128,62,32,21,11,9,4,3,3,0,0,0,0,1,0,0,1)
data=rep(midpt, freq);data
brk=c(0,6,12,18,24,30,36,42,48,54,60,66,72,78,84,90,96,102,108,114,120,126)
hhistt<-hist(data, breaks=brk, col="coral", border="dodgerblue4",lwd=2,xaxp=c(0,126,21),main = "DAP de ·rboles acompaÒantes vs. frecuencia",
             xlab = "CategorÌas de DAP (cm)", ylab = "Frecuencia absoluta", freq = T, las=1, cex.lab=1.4)
lines(c(min(hist$breaks),hist$mids,max(hist$breaks)),c(0,hist$counts,0),type = "l",col="darkblue")
text(hist$mids, hist$counts, labels = hist$counts)
?xaxp
#AREA BASAL#
range(basecero$AB)
brkvec<-seq(0, 1.2, by=0.12)
length(brkvec)
hist3<-hist(basecero$AB, breaks = brkvec, col="gray90", border="gray70",xaxp=c(0,1.2, 10),main = "", 
           xlab = "BA categories (m2 Ha-1)", ylab = "Absolute frequency", freq = T, las=1, cex.lab=1.9)
lines(c(min(hist3$breaks),hist3$mids,max(hist3$breaks)),c(0,hist3$counts,0),type = "l",col="black")
text(hist3$mids, hist3$counts, labels = hist3$counts, cex = 1.55)
?text
par (mfrow = c(1, 3))

#cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5#

# # # # # # # # # # # # # # #
# histogrmas index SIMPSON RIQUEZA ABUNDANCIA 
# # # # # # # # # # # # # # #
rm(list = ls())
graphics.off()
par (mfrow = c(4, 1))
par(mar=c(4.1,4,0.6,1)+0.1)
setwd("C:/Users/Gabriel Torrico/Documents/BiologÌa/PROYECTO COFFEE/Bases de datos")
getwd()
index<-read.csv("index_hist.csv", dec = ",", sep = ";")
range(index$S)
brkvec<-seq(5, 25, by=5)
length(brkvec)
hist3<-hist(index$S, breaks = brkvec, col="gray90", border="gray70",xaxp=c(0,25, 5),main = "", 
            xlab = "Species richness", ylab = "Absolute frequency", freq = T, las=1.3, cex.lab=1.3, cex.axis=1.3)
text(hist3$mids, hist3$counts, labels = hist3$counts, cex = 1.3)
legend("topright", "a", bg="black", text.col="white", adj=0.7, cex = 1.3)
range(index$Abum)
brkvec<-seq(10, 110, by=20)
length(brkvec)
hist3<-hist(index$Abum, breaks = brkvec, col="gray90", border="gray70",xaxp=c(10,110, 5),main = "", 
            xlab = "Abundance of trees", ylab = "Absolute frequency", freq = T, las=1.3, cex.lab=1.3, cex.axis=1.3)
text(hist3$mids, hist3$counts, labels = hist3$counts, cex = 1.3)
legend("topright", "b", bg="black", text.col="white", adj=0.7, cex = 1.3)
range(index$invsimp)
brkvec<-seq(3, 15, by=3)
length(brkvec)
hist3<-hist(index$invsimp, breaks = brkvec, col="gray90", border="gray70",xaxp=c(3,15, 4),main = "", 
            xlab = "Inverse Simpson's diversity index", ylab = "Absolute frequency", freq = T, las=1.3, cex.lab=1.3, cex.axis=1.3)
text(hist3$mids, hist3$counts, labels = hist3$counts, cex = 1.3)
legend("topright", "c", bg="black", text.col="white", adj=0.7, cex = 1.3)
range(index$eqinvsimpson)
brkvec<-seq(0.10, 0.70, by=0.10)
length(brkvec)
hist3<-hist(index$eqinvsimpson, breaks = brkvec, col="gray90", border="gray70",xaxp=c(0.10,0.70, 6),main = "", 
            xlab = "Shannon's evenness index", ylab = "Absolute frequency", freq = T, las=1.3, cex.lab=1.3, cex.axis=1.3)
text(hist3$mids, hist3$counts, labels = hist3$counts, cex = 1.3)
legend("topright", "d", bg="black", text.col="white", adj=0.7, cex = 1.3)

?hist
??cex
par (mfrow = c(1, 1))
range(index$disim)
brkvec<-seq(0, 1, by=0.2)
length(brkvec)
hist3<-hist(index$disim, breaks = brkvec, col="gray90", border="gray70",xaxp=c(0,1, 5),main = "", 
            xlab = "Morisita-Horn dissimilarity index", ylab = "Absolute frequency", freq = T, las=0.5, cex.lab=1)
text(hist3$mids, hist3$counts, labels = hist3$counts, cex = 1)
range(index$Abum)
####       ####
#### C C A ####
####       ####
###############
library(vegan)
library(ggplot2)
library(vegan3d)
graphics.off()
rm(list = ls())
setwd("C:/Users/Gabriel Torrico/Documents/Biolog√≠a/PROYECTO COFFEE/Bases de datos")
getwd()
baseprop<-read.csv("baseprop.csv", sep = ";", dec = ",", header = T)
head(baseprop)
rowSums(baseprop)
baseenv<-read.csv("baseenv.csv", sep = ";", dec = ",", header = T)
baseenv
row.names(baseprop)<-baseprop$Cafi
baseprop$Cafi<-NULL
baseprop
row.names(baseenv)<-baseenv$Cafi
baseenv$Cafi<-NULL
baseenv
all.equal(row.names(baseprop),row.names(baseenv))
# P1 CCA# #sin punttuaci√≥n y sin distancia#
cca1<-cca(baseprop~age+size+masl+slope+shade+cafes+punt+manag+varie+dist, data = baseenv)
vif.cca(cca1)
cca2<-cca(baseprop~age+size+masl+slope+shade+cafes+manag+varie+dist, data = baseenv)
vif.cca(cca2)
cca3<-cca(baseprop~age+size+masl+slope+shade+cafes+manag+varie, data = baseenv)
vif.cca(cca3)
summary(cca3)
plot(cca3)
anova(cca3, permutations = 1000)
ordistep(cca3)
anova(cca3, step=1000, by="terms")
cca4<-cca(baseprop~masl, data = baseenv)
summary(cca4)
anova(cca4, permutations = 1000)
RsquareAdj(cca4)
plot(cca4)
###############
####       ####
####  C A  ####
####       ####
###############
baseprop<-read.csv("baseprop.csv", sep = ";", dec = ",", header = T)
row.names(baseprop)<-baseprop$Cafi
baseprop$Cafi<-NULL
baseprop
ca.prop<-cca(baseprop)
summary(ca.prop)
plot(ca.prop)
ca.prop$CA$u
ca.prop$CA$v
datos.ca.parcelas<-data.frame(ca.prop$CA$u[,1:3])
datos.ca.spps<-data.frame(ca.prop$CA$v[,1:3])
str(datos.ca.parcelas)
head(datos.ca.parcelas)
str(datos.ca.spps)
head(datos.ca.spps)
ggplot()+
  theme_bw()+
  labs(x="CA1 (20%)", y="CA2 (15%)")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18, face = "bold"))+
  geom_text(data = datos.ca.parcelas,
            aes(CA1, CA2, label=row.names(datos.ca.parcelas)),
            size=5, fontface="bold")+
  geom_text(data = datos.ca.spps,
            aes(CA1, CA2, label=row.names(datos.ca.spps)),
            col="red", size=4)

########################################
# # # # # # # # # # # # # # # # # # #  #
#####                             ######
#####  G I G A   C L U S T E R S  ######
#####                             ######
# # # # # # # # # # # # # # # # # # #  #
########################################
library(pvclust)
##############
setwd("C:/Users/Gabriel Torrico/Documents/BiologÌa/PROYECTO COFFEE/Bases de datos")
matabun<-read.csv("matabun2.csv", dec = ",", sep = ";", header = T, row.names = 1)
mata<-as.matrix(matabun)
matadist<-dist(mata, method = "euclidean")
?dist
comm.clust <- hclust(matadist, method = "average")
cutree(comm.clust, k=2:7)
comm3 <- cutree(comm.clust, k=3)
#Make plot labelled with cluster numbers:
plot(comm.clust, labels = as.character(comm3))
#Add cluster names to plot:
comm.lab <- paste(as.character(comm3), names(comm3), sep=":")
plot(comm.clust, labels=comm.lab)
#Highlight the clusters:
rect.hclust(comm.clust, k=3)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
View(sppabun)
dim(sppabun)
matabun<-read.csv("matabun2.csv", dec = ",", sep = ";", header = T, row.names = 1)
head(matabun)
dim(matabun)
datosmat<-as.matrix(matabun)#convertir dataframe en matriz
datosmatdist<-dist(datosmat, method = "canberra")#elaborar la matriz de distancia
datosclustabu<-hclust(datosmatdist, method = "average")#elaborar el cluster con el m√É¬©todo de vinculaci√É¬≥n#
plot(datosclustabu,hang = -1)#gr√É¬°fico del dendrograma#


fviz_dend(x = datosclustabu, cex = 0.95, lwd = 0.75,
          k_colors = c("grey5", "grey40", "grey70", "grey50","grey25","grey85"),
          rect = TRUE, 
          rect_border = "grey35", 
          rect_fill = F, horiz = T)

?pvclust
transabu<-t(matabun)

sorensen<-function(x){
  x<-as.matrix(x)
  y<-t(x)
  mat.dist<-vegdist(y, method = "bray")
  return(mat.dist)
}
pvclustsormed<-pvclust(transabu, method.dist = sorensen,
                       method.hclust = "average")
plot(pvclustsormed, hang=-1, cex = 1)
pvrect(pvclustsormed)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#cluster con datos de presencia ausencia#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
murcis<-read.csv("matriz0y1.csv",dec = ",", sep = ";", header = T, row.names = 1)
head(murcis)
dim(murcis)
head(murcis)
#2.1. Dendrograma con hclust#
datosmurcis<-as.matrix(murcis) #convertir dataframe en matriz
datosmurcisdist<-dist(datosmurcis, method = "binary") #elaborar la matriz de distancia
datosmurcisdist
?vegdist
datosmurcisdist2<-vegdist(datosmurcis, method = "jaccard")
datosmurcisdist2
datosmurcisdist3<-vegdist(datosmurcis, method = "canberra") #<- no funciona para los dastos
datosmurcisdist3
#elaborar el cluster con los 3 m√©todos de vinculaci√≥n#
datosclustmur1<-hclust(datosmurcisdist, method = "single")
datosclustmur2<-hclust(datosmurcisdist, method = "average")
datosclustmur3<-hclust(datosmurcisdist, method = "complete")
#gr√°fico del dendrograma#
plot(datosclustmur1,hang = -1)
plot(datosclustmur2,hang = -1)
library(factoextra)
library(dendextend)
fviz_dend(x = datosclustmur2, cex = 0.95, lwd = 0.75, k = 3,
          k_colors = c("grey5", "grey40", "grey70"),
          rect = TRUE, 
          rect_border = "grey35", 
          rect_fill = F, horiz = T)
plot(datosclustmur3,hang = -1)
#2.2. Dendrograma con pvclust#
### trasponer los datos originales
transmur<-t(murcis)
head(transmur)
#elaborar el cluster con pvclust#
pvclustmur<-pvclust(transmur, method.dist = "binary",
                    method.hclust = "average", nboot = 1000)  #an√°lisis de dendrograma#
plot(pvclustmur, hang=-1)
pvrect(pvclustmur)
#Elaborar cluster con el √≠ndice de similitud de sorensen#
#elaborar funci√≥n para √≠ndices de (dis)similitud#
dist.pvclust<-function(x){
  x<-as.matrix(x)
  y<-t(x)
  mat.dist<-vegdist(y, method = "bray")
  return(mat.dist)
}
#m√©todo de vinculaci√≥nn "vecino m√°s cercano"#
pvclustsorcer<-pvclust(transmur, method.dist = dist.pvclust,
                       method.hclust = "single")
plot(pvclustsorcer, hang=-1)
pvrect(pvclustsorcer)
#m√©todo de vinculaci√≥n average "media"#
pvclustsormed<-pvclust(transmur, method.dist = dist.pvclust,
                       method.hclust = "average")
plot(pvclustsormed, hang=-1)
pvrect(pvclustsormed)
ggsave("ceroyuno.jpg", device = "jpg", height = 12.5, width = 10, units = "cm", dpi = 300)
#m√©todo de vinculaci√≥nn "vecino m√°s lejano"#
pvclustsorlej<-pvclust(transmur, method.dist = dist.pvclust,
                       method.hclust = "complete")
plot(pvclustsorlej, hang=-1)
pvrect(pvclustsorlej)