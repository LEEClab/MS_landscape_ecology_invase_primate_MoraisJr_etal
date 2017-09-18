###########################################################
#
# Binomial models for occupancy of invasive marmosets
#
# Marcio Morais Jr <moraisjr@gmail.com>
# Bernardo Niebuhr <bernardo_brandaum@yahoo.com.br>
#
# Sept. 2017
###########################################################

# Load packages
if(!require(lme4)) { install.packages("lme4", dep = T); library(lme4)}
if(!require(MuMIn)) { install.packages("MuMIn", dep = T); library(MuMIn)}
if(!require(rms)) { install.packages("rms", dep = T); library(rms)}
if(!require(visreg)) { install.packages("visreg", dep = T); library(visreg)}
if(!require(ncf)) { install.packages("ncf", dep = T); library(ncf)}
if(!require(languageR)) { install.packages("languageR", dep = T); library(languageR)}
if(!require(hier.part)) { install.packages("hier.part", dep = T); library(hier.part)}

# Data folder
datadir <- "/home/leecb/Github/ms_landscape_ecology_invase_primate_MoraisJr_etal/data"

# Output folder - for figs etc
outputdir <- "/home/leecb/Github/ms_landscape_ecology_invase_primate_MoraisJr_etal/output"

# Change to data folder
setwd(datadir)

# Read data
presme <- read.table("presme.csv", header=T, sep=",")
head(presme)

# Standardizing continuous variables
presme_std <- cbind(presme[,1:5], apply(X=presme[,6:ncol(presme)], MARGIN=2, FUN=function(x){(x-mean(x))/sd(x)}))

#----------------
# 1) Calculating  correlation between continuous explanatory variables and excluding 
# those correlated variables that explain less the variation in the response variable
cor(presme_std[,6:13], method="spearman")

# The variables “area”, “contig”, “Prox”, and “nneigh” are correlated
# We fit univariate models for each and exclude those with higher AICc
marea <- glmer(presence ~ area + (1|frag), family=binomial, data=presme_std)
mcontig <- glmer(presence ~ contig + (1|frag), family=binomial, data=presme_std)
mprox <- glmer(presence ~ Prox + (1|frag), family=binomial, data=presme_std)
mnneigh <- glmer(presence ~ nneigh + (1|frag), family=binomial, data=presme_std)

model.sel(marea, mcontig)
# Area and Contig models present very similar AICc values. 
# We will keep Area because it has a more direct relation with theories of optimal habitat occupancy (Puliam,1998).
model.sel(mprox, mnneigh)
# Prox e Nneigh models present very similar AICc values
# We will keep Nneigh because it is a direct measure of distance and has a smaller AICc value

# Calculating teh factor of variance inflation to check the remaining correlation
glmVIF <- glm(presence~area+circle+nneigh+nurb+urbarea+nroad, family=binomial,data=presme)
vif(glmVIF)
# No value is too high (>10), so this seems to be ok

#----------------
# 2) Calculating univriate models and a priori models and comparing them

# Calculating univariate models for the 6 variables
marea <- glmer(presence ~ area+(1|frag), family=binomial, data=presme_std)
mcircle <- glmer(presence ~ circle+(1|frag), family=binomial, data=presme_std)
mnneigh <- glmer(presence ~ nneigh+(1|frag), family=binomial, data=presme_std)
mnurb <- glmer(presence ~ nurb+(1|frag), family=binomial, data=presme_std)
mnroad <- glmer(presence ~ nroad+(1|frag), family=binomial, data=presme_std)
murbarea <- glmer(presence ~ urbarea+(1|frag), family=binomial, data=presme_std)

# Calculating models for the 2 hypotheses (anthropogenic facilitation and disturbed habitat)
mantrop <- glmer(presence ~ nurb+nroad+urbarea+(1|frag), family=binomial, data=presme_std)
mhabitat <- glmer(presence ~ area+circle+nneigh+(1|frag), family=binomial, data=presme_std)

# Comparing models
options(na.action = "na.fail")

(antrop <- dredge(mantrop, beta = FALSE, eval = TRUE, rank = "AICc"))

(habitat <- dredge(mhabitat, beta = FALSE, eval = TRUE, rank = "AICc"))

# Removing “urbarea” from mantrop model
mantrop2 <- glmer(presence~nurb+nroad+(1|frag),family=binomial,data=presme_std)

# Removing “patch shape - circle” from mantrop model
mhabitat2 <- glmer(presence ~ area+nneigh+(1|frag), family=binomial, data=presme_std)

# Comparing the models based on a priori hypothesis, to select the most parcimonious one
(compare <- model.sel(mantrop2, mhabitat2, marea, mcircle, mnneigh, mnurb, mnroad, murbarea))
# (all univariate models fit poor in comparison to a priori models)
# we really keep the a priori models here
(compare <- model.sel(mantrop2, mhabitat2))

# Calculating mean coefficients for all models selected
mav_compare <- model.avg(mhabitat2, mantrop2)
# mav_compare <- model.avg(habitat)
# mav_compare <- model.avg(antrop)
importance(mav_compare)

# Plot 3D to look at the effect of variables
visreg2d(mantrop2, "nroad", "nurb", scale="response", xlab="Distance from Roads (m)", ylab="Distance from Urban Areas (m)", plot.type="image")
visreg2d(mantrop2, "nroad","nurb", scale="response", xlab="Distance from Roads (m)", ylab="Distance from Urban Areas (m)", plot.type="persp")
visreg2d(mantrop2,"nroad","nurb", scale="response", xlab="Distance from Roads (m)", 
         ylab="Distance from Urban Areas (m)", main="Probability of occurrence", 
         plot.type="image", col=rainbow(40, s = 1, v = .9, start = 1/25, end = 2.7/6, alpha = 1))

#----------------
# 3) Comparing these results with the model combining all variables from a priori models
# and with the full model

# Model 3 - combination of a priori models mantrop2 and mhabitat2
comb.model <- glmer(presence ~ area + nneigh + nurb + nroad + (1|frag), family=binomial, data=presme_std)

# Model selection including this combination model
(comb.mod <- dredge(comb.model, beta = FALSE, eval = TRUE, rank = "AICc"))

(compare <- model.sel(mantrop2, mhabitat2, comb.model))

# Calculating mean coeffients and importance of each, including the combination model
mav_compare.comb <- model.avg(comb.model, mhabitat2, mantrop2)
importance(mav_compare.comb)

# Plotting the effect of variables for the combination model
par(mfrow = c(2,2))
plotLMER.fnc(comb.model, fun=plogis, addlines=TRUE)
par(mfrow = c(1,1))

# Some conclusions:
# From all these analyses, if we look only at the a priori models, we have a
# stronger effect of antropic (invasion facilition) variables
# If we include the combined model, again the antropic variables seem to present
# an stronger effect (see the model comparison tables and the figures of variable effects)
# Patch area may be important too, but from the graphics it seems that its effect is not so large.

#----------------
# 4) Including a full model - does it make sense?

# Including a full model
full <- glmer(presence ~ area + circle + nneigh + nurb + nroad + urbarea + (1|frag), 
            family=binomial, data=presme_std)

# Selecting variables
(full.select <- dredge(full, beta = FALSE, eval = TRUE, rank = "AICc"))
Totalm <- get.models(full.select, cumsum(weight) <= 1.00)

# Mean coefficient for all selected models
full_mav <- model.avg(Totalm)
importance(full_mav)

# Coefficient of determination based on likelihood test
r.squaredLR(full, null = null.fit(full, TRUE))

# For all models, including patch ID as a random intercept in the null model
mnull <- glmer(presence ~ 1+(1|frag), family=binomial, data=presme_std)
r.squaredLR(full, null = mnull)
r.squaredLR(comb.model, null = mnull)
r.squaredLR(mantrop2, null = mnull)
r.squaredLR(mhabitat, null = mnull)

# Hierarchical partition and GOF
env <- subset(presme_std, select=c(6,7,9,10,11))
gof <- all.regs(presme_std$presence, env, fam = "binomial", gof = "Rsqu", print.vars = TRUE)
partition(gof, pcan = 5, var.names = names(subset(presme_std, select=c(6,7,9,10,11))))

# Does it make sense to include a full model?
# Well, at least we've got an interesting result here: If considering all these
# 5 variables, nurb and nroad explain 66% of all variability together
          
#----------------
# 5) Checking spatial autocorrelation and model assumptions

# Analysis of spatial autocorrelation
Correlog_mantrop2 <- spline.correlog(x=presme[,"E"], y=presme[,"N"], xmax=5000, 
                                     z=resid(mantrop2))
plot.spline.correlog(Correlog_mantrop2)
summary(Correlog_mantrop2)

Correlog_mhabitat2 <- spline.correlog(x=presme[,"E"], y=presme[,"N"], xmax=5000, 
                                     z=resid(mhabitat2))
plot.spline.correlog(Correlog_mhabitat2)
summary(Correlog_mhabitat2)

# completar a partir daqui!! ta interessante isso!

# QUANTILE-QUANTILE PLOTS

# Real expected values and residuals
Fitted <- fitted(mantrop2)
Resid <- residuals(mantrop2)
# calculate the ordered residuals for the model
Resids <- sort(Resid)
Resids_0 <- Resids[Resids<0]
Resids_1 <- Resids[Resids>=0]

# Simulated values
# Replicate the process and produce qqplots
# Specify the number of replicates
Reps <- 1000
# Simulate Reps data sets
set.seed(123)
Sims <- matrix(rbinom(n=length(Fitted)*Reps, size=1, p=Fitted), nrow=length(Fitted), ncol=Reps)
# Fit the model to each simulated data set
cont <- 1
Models <- apply(Sims, MARGIN=2, FUN=function(X) { print(cont); cont <<- cont+1; return(list(X, glmer(formula = X ~ presme_std$nurb + presme_std$nroad + (1|presme_std$frag), family = binomial)))})
# Calculate the ordered (interpolated) simulated residuals and the point-wise 95% confidence intervals#
cont <- 1
Resids_Sim <- matrix(unlist(lapply(Models, FUN=function(X) { print(cont); cont <<- cont+1; TempData <- presme_std; TempData[,"presence"] <- X[1][[1]][[1]]; return(sort(resid(X[2][[1]])))} )), ncol=Reps, nrow=nrow(presme_std))
Resids_Sim_0 <- matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X<0],ppoints(Resids_0,a=1))}),ncol=Reps,nrow=length(Resids_0))
Resids_Sim_1 <- matrix(apply(Resids_Sim,MARGIN=2,FUN=function(X){quantile(X[X>=0],ppoints(Resids_1,a=1))}),ncol=Reps,nrow=length(Resids_1))
Resids_Sim_0_Median <- apply(Resids_Sim_0,MARGIN=1,FUN=median)
Resids_Sim_1_Median <- apply(Resids_Sim_1,MARGIN=1,FUN=median)
Resids_Sim_0_Lower <- apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.025)})
Resids_Sim_0_Upper <- apply(Resids_Sim_0,MARGIN=1,FUN=function(X){quantile(X,0.975)})
Resids_Sim_1_Lower <- apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.025)})
Resids_Sim_1_Upper <- apply(Resids_Sim_1,MARGIN=1,FUN=function(X){quantile(X,0.975)})
# plot the qauntile-quantile plot with 95% confidence intervals and 1:1 line
plot(Resids_Sim_0_Median, Resids_0, xlim = c(-1,1), ylim = c(-1,1), pch = 20,
     xlab = "Simulated quantiles",ylab = "Fitted quantiles",
     col = grey(0, alpha = 0.5))
# plot(sort(Resids), sort(Fitted), col = 2)
points(Resids_Sim_1_Median,Resids_1, pch = 20,
       col = grey(0, alpha = 0.5))
lines(Resids_Sim_0_Median,Resids_Sim_0_Lower)
lines(Resids_Sim_0_Median,Resids_Sim_0_Upper)
lines(Resids_Sim_1_Median,Resids_Sim_1_Lower)
lines(Resids_Sim_1_Median, Resids_Sim_1_Upper)
abline(0, 1, lty=3)

# Export figure
setwd(outputdir)
tiff("qqplot_marmosets.tif", width = 15, height = 15, units = "cm", res = 300)
plot(Resids_Sim_0_Median, Resids_0, xlim = c(-1,1), ylim = c(-1,1), pch = 20,
     xlab = "Simulated quantiles",ylab = "Fitted quantiles",
     col = grey(0.3, alpha = 0.5))
# plot(sort(Resids), sort(Fitted), col = 2)
points(Resids_Sim_1_Median,Resids_1, pch = 20,
       col = grey(0.1, alpha = 0.5))
lines(Resids_Sim_0_Median,Resids_Sim_0_Lower)
lines(Resids_Sim_0_Median,Resids_Sim_0_Upper)
lines(Resids_Sim_1_Median,Resids_Sim_1_Lower)
lines(Resids_Sim_1_Median, Resids_Sim_1_Upper)
abline(0, 1, lty=3)
dev.off()

# PARTIAL RESIDUAL PLOTS

# Calculate the partial residuals for each covariate
Part_Res1 <- Resid/Fitted*(1-Fitted)+mantrop2@beta[2]*presme_std [,"nurb"]
Part_Res2 <- Resid/Fitted*(1-Fitted)+mantrop2@beta[3]*presme_std [,"nroad"]

# Plot the partial residuals and the smoothed plot
setwd(outputdir)
tiff("part_resid_marmosets.tif", width = 20, height = 12, units = "cm", res = 300)
split.screen(c(1,2))
screen(1)
plot(presme[,"nurb"], Part_Res1, 
     xlab="Distance to the nearest urban area (m)", ylab="Partial residuals",
     col = grey(0.3, alpha = 0.5), pch = 19)
lines(seq(min(presme[,"nurb"]),max(presme[,"nurb"]),length.out=100),predict(loess(formula=Part_Res1~nurb,data=presme),newdata=data.frame(nurb=seq(min(presme[,"nurb"]),max(presme[,"nurb"]),length.out=100))))
screen(2)
plot(presme[,"nroad"], Part_Res2, 
     xlab="Distance to the nearest road (m)", ylab="Partial residuals",
     col = grey(0.3, alpha = 0.5), pch = 19)
lines(seq(min(presme[,"nroad"]),max(presme[,"nroad"]),length.out=100),predict(loess(formula=Part_Res2~nroad,data=presme),newdata=data.frame(nroad=seq(min(presme[,"nroad"]),max(presme[,"nroad"]),length.out=100))))
close.screen(all = TRUE)
dev.off()

#----------------
# 6) Predicting the probability of occurrence in most of the fragments

# Reading data
setwd(datadir)
fragdist <- read.table("frag_dist.csv", header=T, sep=",")
head(fragdist)

# Predict and export
mdpredict <- predict(mantrop2, newdata = , type = "response", re.form = NA)
write.xlsx(mdpredict, "predict.xlsx")

x <- seq(from = 0, to = max(fragdist$nroad), by = 100)
y <- seq(from = 0, to = max(fragdist$nurb), by = 1000)
z <- predict()


install.packages("akima", dependencies = T)
library(akima)

predict(mc8)

grid <- interp(x = presme, y = fragdist$nurb, z = predict(mantrop2, newdata = fragdist, type = "response", re.form = NA),
               xo = seq(0, max(fragdist$nroad), by = 100), yo = seq(0, max(fragdist$nurb), by = 500),
               duplicate = "mean")


grid <- interp(x = fragdist$nroad, y = fragdist$nurb, z = predict(mantrop2, newdata = fragdist, type = "response", re.form = NA),
               xo = seq(0, max(fragdist$nroad), by = 100), yo = seq(0, max(fragdist$nurb), by = 500),
               duplicate = "mean")

plot(final$HAMOUNT7.5, final$HSPLIT.res7.5)
image(grid, col = rev(heat.colors(15)), xlab = "% HABITAT (7.5 km)",
       ylab = "Res?duos da DH (7.5 km)", main = "Riqueza de esp?cies")
contour(grid, add=TRUE, labcex = 1.1)
points(final$HAMOUNT7.5, final$HSPLIT.res7.5)

HAMOUNT7.5.100 <- 100*final$HAMOUNT7.5
HSPLIT.res7.5.100 <- 100*final$HSPLIT.res7.5
mc8 <- glm(final$S ~ HAMOUNT7.5.100 + HSPLIT.res7.5.100, family = poisson)
predict(mc8)

grid <- interp(x = HAMOUNT7.5.100, y = HSPLIT.res7.5.100, z = predict.glm(mc8, type = "response"),
               xo = seq(0,100, length = 150), yo = seq(min(HSPLIT.res7.5.100), max(HSPLIT.res7.5.100), length = 150))

png("1resultados_S_tudojunto3.png", 500, 500, res = 100)
plot(final$HAMOUNT7.5, final$HSPLIT.res7.5)
image (grid, col = rev(heat.colors(15)), xlab = "% HABITAT (7.5 km)",
       ylab = "Res?duos da DH (%, 7.5 km)", main = "Riqueza de esp?cies")
filled.contour (grid, col = rev(heat.colors(30)), xlab = "% HABITAT (7.5 km)",
                ylab = "Res?duos da DH (%, 7.5 km)", main = "Riqueza de esp?cies")
contour(grid, add=TRUE, labcex = 1.1, )
points(HAMOUNT7.5.100, HSPLIT.res7.5.100, pch = 21, bg = 1)
dev.off()


