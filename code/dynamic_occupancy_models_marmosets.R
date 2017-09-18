###########################################################
#
# Models of dynamic occupancy for Marmosets
#
# Marcio Morais Jr <moraisjr@gmail.com>
# Bernardo Niebuhr <bernardo_brandaum@yahoo.com.br>
#
# Dec. 2016
###########################################################

# Load packages
if(!require(unmarked)) { install.packages("unmarked", dep = T); library(unmarked)}
if(!require(MuMIn)) { install.packages("MuMIn", dep = T); library(MuMIn)}

# Data folder
datadir <- "/home/leecb/Github/ms_landscape_ecology_invase_primate_MoraisJr_etal/data"

# Change to data folder
setwd(datadir)

# Read data
ocdyn <- read.table("occudyn.csv", header=T, sep=",")

# Explore detection data
detection.dat <- ocdyn[,1:26]
det <- cbind(detection.dat, det = rowSums(detection.dat[,2:26], na.rm = T))
table(det$det) # number of detections
table(det$det)/nrow(det) # percentage of detections

# Standardizing continuous variables
ocdyn_std <- cbind(ocdyn[,1:26], apply(X=ocdyn[,27:ncol(ocdyn)], MARGIN=2, FUN=function(x){ (x-mean(x))/sd(x) }))

# Detection data (y - response variable)
y.marmoset <- as.matrix(ocdyn_std[,2:26])
years <- as.character(2002:2006)
years <- matrix(years, nrow(ocdyn_std), 5, byrow=TRUE)

# Formating data for unmarked
umf <- unmarkedMultFrame(y=y.marmoset, yearlySiteCovs=list(year=years), siteCovs=ocdyn_std[,27:32], numPrimary=5)

####################################################
# 1) Checking for effect of covariates on detection
####################################################

# Fitting models with covariates for detection
det1 <- colext(~1, ~1, ~1, ~area+nneigh+circle+nurb+nroad+urbarea, umf)
det2 <- colext(~1, ~1, ~1, ~area, umf)
det3 <- colext(~1, ~1, ~1, ~nneigh, umf)
det4 <- colext(~1, ~1, ~1, ~circle, umf)
det5 <- colext(~1, ~1, ~1, ~nurb, umf)
det6 <- colext(~1, ~1, ~1, ~nroad, umf)
det7 <- colext(~1, ~1, ~1, ~urbarea, umf)
det8 <- colext(~1, ~1, ~1, ~year, umf)
det9 <- colext(~1, ~1, ~1, ~circle+nurb+nroad+urbarea, umf)
det10 <- colext(~1, ~1, ~1, ~area+circle+nurb+nroad+urbarea, umf)
det11 <- colext(~1, ~1, ~1, ~nneigh+circle+nurb+nroad+urbarea, umf)
det12 <- colext(~1, ~1, ~1, ~area+circle+nurb+nroad+urbarea+year, umf)
det13 <- colext(~1, ~1, ~1, ~circle+nurb+nroad+urbarea+year, umf)
det14 <- colext(~1, ~1, ~1, ~area+nneigh+circle+nurb+nroad+urbarea+year, umf)
det15 <- colext(~1, ~1, ~1, ~nurb+nroad+urbarea+year, umf)
det16 <- colext(~1, ~1, ~1, ~area+nneigh+circle+year, umf)
det17 <- colext(~1, ~1, ~1, ~circle+nurb+nroad+urbarea+year, umf)
det18 <- colext(~1, ~1, ~1, ~nroad+urbarea+year, umf)

# Comparing models for detection
options(na.action = "na.fail")

det_list <- fitList(det1, det2, det3, det4, det5, det6, det7, det8, det9, det10, det11, det12, det13, det14, det15, det16, det17, det18)
det_sel <- modSel(det_list)

# Models fit hypothesis 1 - anthropogenic facilitation for invasion
mantrop_oc1 <- colext(~1, ~1, ~1, ~nurb+nroad+urbarea, umf)
mantrop_oc2 <- colext(~1, ~1, ~1, ~nurb+nroad, umf)
mantrop_oc3 <- colext(~1, ~1, ~1, ~nurb, umf)
mantrop_oc4 <- colext(~1, ~1, ~1, ~nroad, umf)
mantrop_oc5 <- colext(~1, ~1, ~1, ~urbarea, umf)
mantrop_oc6 <- colext(~1, ~1, ~1, ~nroad+urbarea, umf)
mantrop_oc7 <- colext(~1, ~1, ~1, ~nurb+urbarea, umf)
mantrop_oc8 <- colext(~1, ~1, ~1, ~1, umf) 

# Models fit hypothesis 2 - habitat degradation
mhabitat_oc1 <- colext(~1, ~1, ~1, ~area+nneigh+circle, umf)
mhabitat_oc2 <- colext(~1, ~1, ~1, ~area+nneigh, umf)
mhabitat_oc3 <- colext(~1, ~1, ~1, ~area, umf)
mhabitat_oc4 <- colext(~1, ~1, ~1, ~nneigh, umf)
mhabitat_oc5 <- colext(~1, ~1, ~1, ~ circle, umf)
mhabitat_oc6 <- colext(~1, ~1, ~1, ~nneigh+circle, umf)
mhabitat_oc7 <- colext(~1, ~1, ~1, ~area+circle, umf)

# Comparing models for detection according to hypothesis 1 and 2
mdet_list <- fitList(mantrop_oc1, mantrop_oc2, mantrop_oc3, mantrop_oc4, mantrop_oc5, mantrop_oc6, mantrop_oc7, mantrop_oc8, mhabitat_oc1, mhabitat_oc2, mhabitat_oc3, mhabitat_oc4, mhabitat_oc5, mhabitat_oc6, mhabitat_oc7)
mdet_sel <- modSel(mdet_list)

# adding year
mantrop_y0 <- colext(~1, ~1, ~1, ~year, umf)
mantrop_y1 <- colext(~1, ~1, ~1, ~nurb+nroad+urbarea+year, umf)
mantrop_y4 <- colext(~1, ~1, ~1, ~nroad+year, umf)
mantrop_y5 <- colext(~1, ~1, ~1, ~urbarea+year, umf)
mantrop_y6 <- colext(~1, ~1, ~1, ~nroad+urbarea+year, umf)

mdet_list <- fitList(mantrop_y0, mantrop_y1, mantrop_y6, mantrop_y4, mantrop_y5)
mdet_sel <- modSel(mdet_list)





# Modelo efeito antropico na ocorrencia e Year na detecção

mantrop_oc1<-colext(~nurb+nroad+urbarea, ~1, ~1, ~year, umf)
mantrop_oc2<-colext(~nurb+nroad, ~1, ~1, ~year, umf)
mantrop_oc3<-colext(~nurb, ~1, ~1, ~year, umf)
mantrop_oc4<-colext(~nroad, ~1, ~1, ~year, umf)
mantrop_oc5<-colext(~urbarea, ~1, ~1, ~year, umf)
mantrop_oc6<-colext(~nroad+urbarea, ~1, ~1, ~year, umf)
mantrop_oc7<-colext(~nurb+urbarea, ~1, ~1, ~year, umf)
mantrop_oc8<-colext(~1, ~1, ~1, ~year, umf)


#Seleção dos modelos
mantrop_list<-fitList(mantrop_oc1,mantrop_oc2,mantrop_oc3,mantrop_oc4,mantrop_oc5,mantrop_oc6,mantrop_oc7,mantrop_oc8)
mantrop_sel<-modSel(mantrop_list)




# Modelo efeito habitat na ocorrencia
mhabitat_oc1<-colext(~area+nneigh+circle, ~1, ~1, ~year, umf)
mhabitat_oc2<-colext(~area+nneigh, ~1, ~1, ~year, umf)
mhabitat_oc3<-colext(~area, ~1, ~1, ~year, umf)
mhabitat_oc4<-colext(~nneigh, ~1, ~1, ~year, umf)
mhabitat_oc5<-colext(~circle, ~1, ~1, ~year, umf)
mhabitat_oc6<-colext(~nneigh+circle, ~1, ~1, ~year, umf)
mhabitat_oc7<-colext(~area+circle, ~1, ~1, ~year, umf)
mhabitat_oc8<-colext(~1, ~1, ~1, ~year, umf)

mhabitat_list<-fitList(mhabitat_oc1,mhabitat_oc2,mhabitat_oc3,mhabitat_oc4,mhabitat_oc5,mhabitat_oc6,mhabitat_oc7,mhabitat_oc8)
mhabitat_sel<-modSel(mhabitat_list)

nnPars    AIC delta AICwt cumltvWt
mhabitat_oc3     9 292.35 0.000 0.273     0.27
mhabitat_oc2    10 292.37 0.024 0.270     0.54
mhabitat_oc8     8 294.05 1.706 0.116     0.66
mhabitat_oc7    10 294.19 1.844 0.109     0.77
mhabitat_oc1    11 294.37 2.023 0.099     0.87
mhabitat_oc4     9 295.23 2.884 0.065     0.93
mhabitat_oc5     9 296.04 3.691 0.043     0.97
mhabitat_oc6    10 297.06 4.708 0.026     1.00



#Comparação dos modelos mais plausiveis para ocorrencia

mod_list<-fitList(mantrop_oc2, mantrop_oc1, mhabitat_oc3, mhabitat_oc2, mhabitat_oc8, mhabitat_oc7)
mod_sel<-modSel(mod_list)

nPars    	AIC		 delta AIC	wt 	cumltvWt
mantrop_oc2     10 288.87 	 0.00 		0.411     0.41
mantrop_oc1     11 289.11  	0.25 		0.364     0.77
mhabitat_oc3     9 292.14 	 3.27 		0.080     0.85
mhabitat_oc2    10 292.17  	3.30 		0.079     0.93
mhabitat_oc8     8 293.85  	4.98 		0.034     0.97
mhabitat_oc7    10 293.99 	 5.12 		0.032     1.00

#Estimativa do coeficiente médio

modavg_oc<-model.avg(mantrop_oc2, mantrop_oc1, mhabitat_oc3, mhabitat_oc2, mhabitat_oc8, mhabitat_oc7, revised.var = TRUE,rank="AIC")
summary(modavg_oc)

Relative variable importance: 
  psi(year) 	p(nroad) p(nurb) 	p(urbarea)	 p(area) 	p(nneigh)		p(circle)
Importance:          1.00      0.77     0.77   		 0.36       	0.19    		0.08      		0.03     

#Model Fit

chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}
(pb <- parboot(mantrop_oc2, statistic=chisq, nsim=100))
Parametric Bootstrap Statistics:
  t0 mean(t0 - t_B) StdDev(t0 - t_B) Pr(t_B > t0)
1 288          -13.2             26.1        0.683
We fail to reject the null hypothesis, and conclude that the model fit is adequate.



#Predições com variáveis não padronizadas

#Ler tabela de dados:#
ocdyn<-read.table("occudyn.csv",header=T,sep=",")

#dados detecção#

y.marmoset_pred <- as.matrix(ocdyn[,2:26])

years_pred <- as.character(2002:2006)
years_pred <- matrix(years_pred, nrow(ocdyn), 5, byrow=TRUE)

#Formatar dados unmarked

umf_pred<- unmarkedMultFrame(y=y.marmoset_pred, yearlySiteCovs=list(year=years_pred), siteCovs=ocdyn[,27:32], numPrimary=5)

#Ajustar o melhor modelo (mantrop_oc2)

mbest<-colext(~nurb+nroad, ~1, ~1, ~year, umf_pred)

#Graficos de predição

op <- par(mfrow=c(1,2), mai=c(0.8,0.8,0.1,0.1))
ndr <- data.frame(nroad=seq(0, 1000, length=100), nurb=0)
E.psir <- predict(mbest, type="psi", newdata=ndr, appendData=TRUE)
with(E.psir, {
  plot(nroad, Predicted, ylim=c(0,1), type="l",
       xlab="Distance from road",
       ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
  lines(nroad, Predicted+1.96*SE, col=gray(0.7))
  lines(nroad, Predicted-1.96*SE, col=gray(0.7))
})
ndu<- data.frame(nroad=0, nurb=seq(0, 10000, length=100))
E.psiu<- predict(mbest, type="psi", newdata=ndu, appendData=TRUE)
with(E.psiu, {
  plot(nurb, Predicted, ylim=c(0,1), type="l",
       xlab="Distance from urban area",
       ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
  lines(nurb, Predicted+1.96*SE, col=gray(0.7))
  lines(nurb, Predicted-1.96*SE, col=gray(0.7))
})
par(op)

op <- par(mfrow=c(1,1), mai=c(0.8,0.8,0.1,0.1))
nd <- data.frame(year=c('2002','2003','2004','2005','2006'))
E.det <- predict(mbest, type='det', newdata=nd)
with(E.det, { # Plot for detection probability: note 5 years
  plot(1:5, Predicted, pch=1, xaxt='n', xlab='Year',
       ylab=expression(paste('Detection probability ( ', p, ' )')),
       ylim=c(0,1), col=4)
  axis(1, at=1:5, labels=nd$year)
  arrows(1:5, lower, 1:5, upper, code=3, angle=90, length=0.03, col=4)
  points((1:5)-0.1, p, col=1, lwd = 1, pch=16)
  legend(7.5, 1, c('Parameter','Estimate'), col=c(1,4), pch=c(16, 1),
         cex=0.8)
})
par(op)


#Ocupação nos 5 anos
mbest <- nonparboot(mbest, B = 100)
cbind(psi=”psi”, smoothed=smoothed(mbest)[2,], SE=mbest@smoothed.mean.bsse[2,])
psi   	smoothed            	SE                  
1 "psi" "0.69230774625669"  "0.108775672652512" 
2 "psi" "0.693662247201891" "0.0942109090734682"
3 "psi" "0.847182488267696" "0.0696385231022667"
4 "psi" "0.884212381380721" "0.0745510610059246"
5 "psi" "0.885148657710739" "0.0638126215389604"


#predição da probabilidade de ocorrencia nos demais fragmentos

fragdist<-read.table("frag_dist.csv",header=T,sep=",")
attach(fragdist)
mpredict<- predict(mbest, type="psi", newdata=fragdist, appendData=TRUE)
write.xlsx(mpredict, "predicted.xlsx")
Análise de ocupação - MLD

mldoccu<-read.table("OccuData.csv", header=T,sep=",")

mld<-unmarkedFrameOccu(y=as.matrix(mldoccu[,c("R1","R2","R3","R4")]),siteCovs=mldoccu[,c("MEANELEV", "MEDDN")])

#Ajustar modelo

M_mld<- occu(~1~1, knownOcc=c(1,2,3,4,5,6,7,8,9,10,11,12,13,141,5,16,17,18,58,59,60), mld)

#Predição para conjunto com todos quadrados (GridTotal.csv)

beta <- coef(mld6, type="state")

gridT<-read.table("GridTotal.csv", header=T,sep="\t")

logit.psi <- beta[1] + beta[2]*gridT$MEANELEV + beta[3]*gridT$MEDDN

psi <- exp(logit.psi) / (1 + exp(logit.psi)) #probabilidade de ocupação

#plot(psi, col=terrain.colors(100))

#print(spplot(psi, col.regions=terrain.colors(100)))

grid<-unmarkedFrameOccu(y=as.matrix(gridT[,c("XMIN","XMAX","YMIN","YMAX")]),siteCovs=gridT[,c("MEANELEV","MEDDN")])

E.psi <- predict(mld6, type="state", newdata=grid) #probabilidade de ocupação e SE, Limite de confiança

write.csv(E.psi, file = "OccuPredict.csv")

write.csv(gridT, file = "Ref.csv")

#probabilidade de detecção para os locais amostrados

betadet <- coef(mld6, type="det")
logit.p<- betadet[1] + betadet[2]*gridT$MEANELEV
p<- exp(logit.p) / (1 + exp(logit.p))
E.p<- predict(mld6, type="det", newdata=grid) #Probabilidade de detecção para todas observações
write.csv(p, file = "DetectPredict.csv")


#Limite para determinar presença ou ausencia de mld
pa<-read.table("pamld.csv",header=T,sep=",")
library(PresenceAbsence)

auc.roc.plot(pa)

#Influencia do sagui na ocorrencia do MLD

mldoccu<-read.table("OccuData.csv", header=T,sep=",")

mld<-unmarkedFrameOccu(y=as.matrix(mldoccu[,c("R1","R2","R3","R4")]),siteCovs=mldoccu[,c("MEANELEV", "MEDDN","ME")])

ME_mld<- occu(~1~ME, knownOcc=c(1,2,3,4,5,6,7,8,9,10,11,12,13,141,5,16,17,18,58,59,60), mld)



# Occupancy model with imperfect detection - binomial-binomial model for
# presence and detection processes

# Análise de ocupação Sagui – single season

# Read data
meoccu <- read.table("occume.csv", header=T,sep=",")
head(meoccu)

# Explore detection data
detection.dat <- meoccu[,2:6]
det <- cbind(meoccu[,1], detection.dat, det = rowSums(detection.dat, na.rm = T))
table(det$det) # number of detections
table(det$det)/nrow(det) # percentage of detections

# Standardize continuous variables
meoccu_std <- cbind(meoccu[,1:7], apply(X=meoccu[,8:ncol(meoccu)], MARGIN=2, FUN=function(x){(x-mean(x))/sd(x)}))

# Transform it into an unmarkedDataFrame 
umf<-unmarkedFrameOccu(y = as.matrix(meoccu_std[,c("det1","det2","det3","det4","det5")]),
                       siteCovs = meoccu_std[,c("area", "nneigh","circle", "nurb","nroad", "urbarea", "year")])
summary(umf)

occu(~1 ~1, umf)
backTransform(occu(~1 ~1, umf), type = 'state')
backTransform(occu(~1 ~1, umf), type = 'det')

#modelar detectabilidade



mdet <- occu(~ area+nneigh+circle+nroad+nurb+urbarea+year~1, umf)
View(mdet)

#Seleção de modelo (mais parcimonioso) entre todas as combinações possíveis#
library (MuMIn)
options(na.action = "na.fail")

ddet<-dredge(mdet, beta = FALSE, eval = TRUE, rank = "AICc")
View(ddet)

#modelar ocupação com melhor combinação de variaveis 

moccu <- occu(~1 ~area+nneigh+circle+nroad+nurb+urbarea, umf)

#Seleção de modelo (mais parcimonioso) de ocupação entre todas as combinações possíveis#

doccu<-dredge(moccu, beta = FALSE, eval = TRUE, rank = "AICc")

#Ajuste do melhor modelo#

mbest<- occu(~ circle+nroad+nurb+urbarea+area+nneigh ~circle+nroad+nurb+urbarea+area+nneigh, umf)

dbest<-dredge(mbest, beta = FALSE, eval = TRUE, rank = "AICc")
View(dbest)

mantrop.occu <- occu(~nroad+nurb ~nroad+nurb, umf)
mhabitat.occu <- occu(~area+nneigh ~area+nneigh, umf)
mmix.occu <- occu(~nroad+nurb ~area+nneigh, umf)
mmix2.occu <- occu(~area+nneigh ~nroad+nurb, umf)

compare.occu <- model.sel(mantrop.occu, mhabitat.occu, mmix.occu, mmix2.occu)
