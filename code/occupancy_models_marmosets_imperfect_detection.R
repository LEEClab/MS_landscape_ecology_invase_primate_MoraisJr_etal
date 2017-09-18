###########################################################
#
# Binomial-binomial models for occupancy of invasive marmosets
#
# Models consider marmoset presence and detection as
# two binomial variables with covariates
#
# Marcio Morais Jr <moraisjr@gmail.com>
# Bernardo Niebuhr <bernardo_brandaum@yahoo.com.br>
#
# Sept. 2017
###########################################################

# Load packages
if(!require(unmarked)) { install.packages("unmarked", dep = T); library(unmarked)}
if(!require(MuMIn)) { install.packages("MuMIn", dep = T); library(MuMIn)}

# Data folder
datadir <- "/home/leecb/Github/ms_landscape_ecology_invase_primate_MoraisJr_etal/data"

# Change to data folder
setwd(datadir)

# Read data - prepared to consider colonization and extinction - we're not doing that!
ocdyn <- read.table("occudyn.csv", header=T, sep=",")

#----------------
# 1) Explore detection data to understand it
detection.dat <- ocdyn[,1:26]
det <- cbind(detection.dat, det = rowSums(detection.dat[,2:26], na.rm = T))
table(det$det) # number of sampling units in which marmosets were detected in each of 0-5 years
table(det$det)/nrow(det) # percentage of sampling units in which marmosets were detected in each of 0-5 years

# From all places where detection was different from 1 or 0, marmosets appeared 
# in the last 2 or 3years

#----------------
# 2) Occupancy model with imperfect detection - binomial-binomial model for
# presence and detection processes
# If colonization happened, it is included in the detection process

# Single season analysis of occupancy - years are treated as independent

# Read data
meoccu <- read.table("occume.csv", header=T,sep=",")
head(meoccu)

# Explore detection data
detection.dat <- meoccu[,2:6]
det2 <- cbind(meoccu[,1], detection.dat, det = rowSums(detection.dat, na.rm = T))
table(det2$det) # number of detections
table(det2$det)/nrow(det2) # percentage of detections
colSums(detection.dat, na.rm = T) # when individuals were detected
colSums(detection.dat, na.rm = T)/sum(detection.dat, na.rm = T)

# Standardize continuous variables
meoccu_std <- cbind(meoccu[,1:7], apply(X=meoccu[,8:ncol(meoccu)], MARGIN=2, FUN=function(x){(x-mean(x))/sd(x)}))

# Transform it into an unmarkedDataFrame 
umf<-unmarkedFrameOccu(y = as.matrix(meoccu_std[,c("det1","det2","det3","det4","det5")]),
                       siteCovs = meoccu_std[,c("area", "nneigh","circle", "nurb","nroad", "urbarea", "year")])
summary(umf)

# Exploring the joint effects of presence and detection
m0 <- occu(~1 ~1, umf)
backTransform(m0, type = 'state')
backTransform(m0, type = 'det')

options(na.action = "na.fail")

# Modeling the occupancy considering the variables in the most parcimonious models
# when we did not consider imperfect detection
mantrop.occu.0 <- occu(~1 ~nroad+nurb, umf)
mantrop.occu <- occu(~nroad+nurb ~nroad+nurb, umf)
mhabitat.occu.0 <- occu(~1 ~area+nneigh, umf)
mhabitat.occu <- occu(~area+nneigh ~area+nneigh, umf)
mmix.occu <- occu(~nroad+nurb ~area+nneigh, umf)
mmix2.occu <- occu(~area+nneigh ~nroad+nurb, umf)

(compare.occu <- model.sel(mantrop.occu, mhabitat.occu, mmix.occu, mmix2.occu, 
                          mantrop.occu.0, mhabitat.occu.0))

# Although some models that consider detection depends on covariates (eg mmix.occu) were 
# between the most parcimonuous, the estimates for detection parameters are not significant

#----------------
# 3) Occupancy models including colonization

# Standardizing continuous variables
ocdyn_std <- cbind(ocdyn[,1:26], apply(X=ocdyn[,27:ncol(ocdyn)], MARGIN=2, FUN=function(x){ (x-mean(x))/sd(x) }))

# Detection data (y - response variable)
y.marmoset <- as.matrix(ocdyn_std[,2:26])
years <- as.character(2002:2006)
years <- matrix(years, nrow(ocdyn_std), 5, byrow=TRUE)

# Formating data for unmarked
umf.dyn <- unmarkedMultFrame(y=y.marmoset, yearlySiteCovs=list(year=years), 
                             siteCovs=ocdyn_std[,27:32], numPrimary=5)
summary(umf.dyn)

mnull <- colext(~1, ~1, ~1, ~1, umf.dyn) 
myear_det1 <- colext(~1, ~1, ~1, ~year, umf.dyn)
mantrop_det1 <- colext(~1, ~1, ~1, ~nurb+nroad, umf.dyn)
mhabitat_det1 <- colext(~1, ~1, ~1, ~area+nneigh, umf.dyn)

mantrop_oc1 <- colext(~nurb+nroad, ~1, ~1, ~1, umf.dyn)
mhabitat_oc1 <- colext(~area+nneigh, ~1, ~1, ~1, umf.dyn)
mantrop_oc2 <- colext(~nurb+nroad, ~1, ~1, ~year-1, umf.dyn)
mhabitat_oc2 <- colext(~area+nneigh, ~1, ~1, ~year, umf.dyn)
mantrop_oc3 <- colext(~nurb+nroad, ~1, ~1, ~nurb+nroad, umf.dyn)
mhabitat_oc3 <- colext(~area+nneigh, ~1, ~1, ~area+nneigh, umf.dyn)
mantrop_oc4 <- colext(~nurb+nroad, ~1, ~1, ~area+nneigh, umf.dyn)
mhabitat_oc4 <- colext(~area+nneigh, ~1, ~1, ~nurb+nroad, umf.dyn)

mlist <- fitList(mnull, myear_det1, mantrop_det1, mhabitat_det1,
                 mantrop_oc1, mhabitat_oc1,
                 mantrop_oc2, mhabitat_oc2,
                 mantrop_oc3, mhabitat_oc3,
                 mantrop_oc4, mhabitat_oc4)
(mlist_sel <- modSel(mlist))

backTransform(mantrop_oc2, type = "col")
backTransform(mantrop_oc2, type = "ext")
sapply(mantrop_oc2@estimates@estimates$det@estimates, function(X) exp...)
backTransform(mantrop_oc2, type = "det") # arrumar
backTransform(mantrop_oc2, type = "state") # arrumar
