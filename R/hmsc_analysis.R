# haiyang learn the Jonint species distribution modelling with R package Hmsc
## appendix S2 from Tikhonov 2019 but again modified by haiyang 20201110

## HMSC use Bayesian interence to fit latent-variable joint species distribution models. Check Ovaskainen concept figure for detials.

# step 1: setting model struture and fitting the model
# step 2: examing MCMC convergence
# step 3: evalating model fit
# step 4: exploring parameter estimates
# step 5: make predictions

## the example dataset is the bird communites in Finland.
rm(list=ls())
setwd("/Users/zhanghaiyang/GitHub/genomesize_transect")
library(Hmsc)


## step 1 read in all the data, species, env, phylo, and traits data

#####################
## PHYLOGENY
phyloTree <- ape::read.tree("Data/phylo_2020.txt")

#####################

#####################
# species data
comm <- read.csv("Data/2020_COMM_data_median_all.csv")
comm <- comm[, -177]
colnames(comm) <- sub('\\.', '_', colnames(comm))
## how to take the dataset to present/absent
comm.data_present <- as.matrix(comm[, -c(1,2)])
rownames(comm.data_present) <- comm$SiteID
comm.data_present[comm.data_present > 0] <- 1 ## here we transfer abundance-based community to absent/present data.
Y1 <- as.matrix(comm.data_present)

library(picante)
combined <- match.phylo.comm(phyloTree, Y1)
###here we have a problem, combined$phylo, because the default name for "combined" was phy
phylo <- combined$phy
Y <- combined$com
#####################

# TRAITS
##################### add gs_species into our dataset
data_gs <- read.csv("Data/2020_COMM_gs_traits.csv")
#data_gs <- na.omit(data_gs)
library(dplyr)
GS_spe <- data_gs %>%
  group_by(new.names) %>%
  dplyr::summarise(GS <- mean(GS_1C_mean, na.rm = T))
GS_spe <- na.omit(GS_spe)
colnames(GS_spe) <- c("Species", "GS")
GS_spe$Species <- sub('\\ ', '_', GS_spe$Species)
TrData = as.data.frame(subset(GS_spe, Species %in% phylo$tip.label))
##########################################

#####################
# locaiton and environmental data
space <- read.csv("Data/2020_spatial_site.csv", header=T)
env <- read.csv("Data/2020_env_site.csv", header = T)
env_res <- env[, c("MAP", "MAT", "Alt","soilN", "soilP")]
XData <- env_res
#####################





## step 2 set up the model and fitting the model
####################################################################################
####################################################################################
####################################################################################

# STUDY DESIGN
studyDesign = matrix(NA,nrow(Y),1)
studyDesign[,1] = rownames(Y)
studyDesign = as.data.frame(studyDesign)
colnames(studyDesign) = c("Site")
studyDesign[,1]=as.factor(studyDesign[,1])
rownames(studyDesign) =  rownames(Y)

# RANDOM EFFECT STRUCTURE, HERE ROUTE AS A SPATIAL LATENT VARIABLE

sRL = space[, c("Longitude", "Latitude")]
colnames(sRL) = c("x","y")
rownames(sRL) = studyDesign[,1]
rL = HmscRandomLevel(sData=sRL, longlat = T)
rL$nfMin = 5
rL$nfMax = 10

XFormula = ~ MAP + MAT + Alt +soilN + soilP
TrFormula = ~ GS

m = Hmsc(Y = Y,
         XData = XData, XFormula = XFormula,
         TrData = TrData, TrFormula = TrFormula,
         phyloTree = phyloTree,
         distr = "probit",
         studyDesign = studyDesign, ranLevels = list(Site=rL))



thin = 10
samples = 1000
nChains = 2

set.seed(1123)
#ptm = proc.time()
m = sampleMcmc(m, samples = samples, thin = thin,
               adaptNf = rep(ceiling(0.4*samples*thin),1),
               transient = ceiling(0.5*samples*thin),
               nChains = nChains, 
               initPar = "fixed effects")

ModelDir = "/Users/zhanghaiyang/GitHub/2019_transect_GS/hmsc_model"
MixingDir = file.path(ModelDir, "mixing")
MFDir = file.path(ModelDir, "model_fit")

filename = file.path(ModelDir, paste("model_", "thin_", ... = as.character(thin),
                                     "_samples_", as.character(samples), ".Rdata", sep = ""))


save(m, file=filename)
#load("/Users/haiyangzhang/haiyangwork/GitHub/2019_transect_GS/hmsc_model/model_1_pa_thin_20_samples_300.Rdata")

## step 3 examing MCMC convergence
####################################################################################
####################################################################################
####################################################################################

# a key step in MCMC-based bayesian analyasis is to assess the performance of the smapling procedure by evaluating chain mixing and 
# convergence. This can be easyily done in Hmsc using tools from the coda package. Within Hmsc, THE FUNCTION convertToCodaObject converts
# the posterior distributions produced by Hmsc into coda-format.

# compute mixing statiscs

thin = 10
samples = 1000
nChains = 2
filename = file.path(ModelDir, paste("model_", "thin_", ... = as.character(thin),
                                         "_samples_", as.character(samples), ".Rdata", sep = ""))
#filename = "/Users/zhanghaiyang/GitHub/2019_transect_GS/hmsc_model/rareremove_model_thin_20_samples_1000.Rdata"
load(filename)
mpost = convertToCodaObject(m)
## we can plot the trace plots of the beta-para or quantitively, we can check the es.beta and ge.v below
## It will be good if effective sampel sizes are high, and potential scale reduction factors are close to one. 

es.beta = effectiveSize(mpost$Beta)
ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
hist(ge.beta)
mean(ge.beta)
median(ge.beta)
quantile(ge.beta, probs = 0.05) 
quantile(ge.beta, probs = 0.95) 

es.gamma = effectiveSize(mpost$Gamma)
ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
hist(ge.gamma)

es.rho = effectiveSize(mpost$Rho)
ge.rho = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
hist(ge.rho)
es.V = effectiveSize(mpost$V)
ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
hist(ge.V)
#es.omega = effectiveSize(mpost$Omega[[1]])
#ge.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
    
mixing = list(es.beta=es.beta, ge.beta=ge.beta,
                  es.gamma=es.gamma, ge.gamma=ge.gamma,
                  es.rho=es.rho, ge.rho=ge.rho,
                  es.V=es.V, ge.V=ge.V#,
                  #es.omega=es.omega, ge.omega=ge.omega
              )
filename = file.path(MixingDir, paste("mixing_", "chains_",as.character(nChains),
                                          "_thin_", as.character(thin),"_samples_",
                                          as.character(samples),
                                          ".Rdata",sep = ""))
save(file=filename, mixing)

    
    
## step 4 evalue the model fit
####################################################################################
####################################################################################
####################################################################################

# code below used to evalue the model's explanatory power.    
set.seed(1123)
predY = computePredictedValues(m, expected=FALSE)
MF = evaluateModelFit(hM=m, predY=predY)
mean(MF$TjurR2) # 41%
mean(MF$AUC) ## 0.91
## we can also evaluate the model's predictive power through two-fold cross validation.
partition = createPartition(m, nfolds = 2)
preds_new = computePredictedValues(m, partition = partition, nParallel = 4)
MF_pred <- evaluateModelFit(hM = m, predY = preds_new)


filename = file.path(MFDir, paste("model_", "_chains_",as.character(nChains),
                                      "_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      ".Rdata",sep = ""))
save(file=filename, MF)


# Predictive power
    
    
## step 4 evalue the model fit
####################################################################################
####################################################################################
####################################################################################
Spec_select <- as.vector(t(TrData[order(TrData$GS),])[1,])

    
### first figure
postBeta = getPostEstimate(m, parName="Beta")
source("hmsc_model/plotBeta_without_intercept.R") ## here I modify the original plotBeta code and generate a new one for removing the intercept.
m$covNames <- c("(Intercept)", "MAP", "MAT", "Altitude","Soil N", "Soil P")
png("Output/fig_401.png", width = 3.5, height = 7, units = 'in', res = 350)
plotBeta_intercept(m, post = postBeta, param = "Support", plotTree = T, supportLevel = 0.95, spNamesNumbers=c(F,F), SpeciesOrder = "Original")
dev.off()

#postGamma = getPostEstimate(m, parName = "Gamma")
#plotGamma(m, post=postGamma, param="Support", supportLevel = 0.5)
#plotGamma(m, post=postGamma, param="Mean", supportLevel = 0.9)
#plotGamma(m, post=postGamma, param="Support", trNamesNumbers=c(TRUE,TRUE))



VP = computeVariancePartitioning(m, group = c(1,1,2,2), groupnames = c("climate","habitat"))
VP_fig_data <- as.data.frame(VP$vals[, Spec_select])
rownames(VP_fig_data) <- c("Climate", "Habitat", "Random")
VP_fig_data_sub <- cbind(VP_fig_data[,(VP_fig_data[3,]) > 0.2], VP_fig_data[, c("Leymus_chinensis", "Agropyron_cristatum")])


##############################################################################
######figure for making variaiton plot for all the species
##############################################################################
VP_fig_data2 <- VP_fig_data
#colnames(VP_fig_data2) <- c(1:169)
cols = heat.colors(3, alpha = 1)

leg = rownames(VP_fig_data2)
means = round(100 * rowMeans(VP_fig_data2), 1)

for (i in 1:3) {
  leg[i] = paste(leg[i], " (", toString(means[i]),
                 "%)", sep = "")
}

png("Output/fig_402_01.png", width = 16, height = 8, units = 'in', res = 300)
par(xpd=T, mar=par()$mar+c(6,0.5,0,0))
barplot(as.matrix(VP_fig_data2), main = "Variation explanation",
        xlab= "", 
        ylab = "Variance proportion", 
        las = 2, col = cols, cex.axis=0.5, cex.names=0.5)
legend(x = "top", legend = rownames(VP_fig_data_sub), fill = cols, bty = "n", ncol = 3, inset = -0.35)
abline(v = 99.5, col = 2,lty=2) # 1.57 Artemisia_scoparia
abline(v = 138, col = 2,lty=2)  # 2.5 Trachydium_tianshanicum
abline(v = 167, col = 2,lty=2)  # 4.03 Hemerocallis_citrina
abline(v = 170.5, col = 2,lty=2) # 5.13 Filifolium_sibiricum
# Now, define a custom axis
dev.off()

#write.csv(TrData, "Data/TrData")
#TrData_sub <- subset(TrData, GS < 5.13)


##############################################################################
######figure for making variaiton plot for 20 species
##############################################################################
png("Output/fig_402_02.png", width = 5, height = 5, units = 'in', res = 250)
par(xpd=T, mar=par()$mar+c(7,0,0,0))
barplot(as.matrix(VP_fig_data_sub), main = "",
        xlab= "", 
        ylab = "Variance proportion", 
        las = 2, col = cols)
legend(x = "top", legend = rownames(VP_fig_data_sub), fill = cols, bty = "n", ncol = 3, inset = -0.35)
dev.off()


#plotVariancePartitioning(m, VP = VP)
library(knitr)
kable(VP$R2T$Beta)
VP$R2T$Y



OmegaCor = computeAssociations(m)
supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel)
          + (OmegaCor[[1]]$support<(1-supportLevel))>0)*OmegaCor[[1]]$mean
#str(toPlot)
toPlot_new <- toPlot[Spec_select, Spec_select]
colnames(toPlot_new) <- c(1:81)
rownames(toPlot_new) <- c(1:81)

library(corrplot)

png("Output/fig_403.png", width = 4.5, height = 3, units = 'in', res = 300)
par(mar=c(1,0,0,0))
corrplot(toPlot_new, method = "color",
         col=colorRampPalette(c("blue","white","red"))(200),
         tl.cex = 0.6,                                   
         tl.srt = 45, tl.col="black", tl.pos='n',
         title="",type="lower", mar=c(0,0,1,0))
dev.off()

## to check the phyloegenetic signal 
summary(mpost$Rho) ## look at the quantiles, whether the interval includes zero or not. If yes, no phylogenetic signal.


####################################### prediction ##########################
## visualizing the results via ploting how the commynity chagnes over some environmental graident of interest.

# this can be done via HMSC by constructing environmental gradients with the funciton constructGradient
# then predicting communities over thoes gradient by the function predict
# Finally using plotGradient to visualize the predicted varations.

Gradient1 = constructGradient(m, focalVariable = "MAP")
Gradient2 = constructGradient(m, focalVariable = "MAT")
Gradient2_2 = constructGradient(m, focalVariable = "Alt")
Gradient3 = constructGradient(m, focalVariable = "soilN")
Gradient4 = constructGradient(m, focalVariable = "soilP")
#Gradient$XDataNew

# as climate is a continuous covariate, the constructed gradient involves a grid of its values ranging from the smallest to the largest values.
# We have set here the value of the non-focal variable habitat to open habitat, and we thus imagine that we sample
# different parts of the climatic gradient solely in open habitats.

predY1 = predict(m, XData = Gradient1$XDataNew, studyDesign = Gradient1$studyDesignNew,
                ranLevels = Gradient1$rLNew, expected = T)

predY2 = predict(m, XData = Gradient2$XDataNew, studyDesign = Gradient2$studyDesignNew,
                 ranLevels = Gradient2$rLNew, expected = T)

predY2_2 = predict(m, XData = Gradient2_2$XDataNew, studyDesign = Gradient2_2$studyDesignNew,
                 ranLevels = Gradient2_2$rLNew, expected = T)

predY3 = predict(m, XData = Gradient3$XDataNew, studyDesign = Gradient3$studyDesignNew,
                 ranLevels = Gradient3$rLNew, expected = T)

predY4 = predict(m, XData = Gradient4$XDataNew, studyDesign = Gradient4$studyDesignNew,
                 ranLevels = Gradient4$rLNew, expected = T)

png("Output/fig_404.png", width = 1.25, height = 5, units = 'in', res = 250)
par(mar=c(2,2,1,1), mfrow=c(5,1))
plotGradient(m, Gradient1, pred = predY1, measure = "T", index = 2, pointsize = 0.2, pointcol="red", showData = T, xlabel = "MAP (mm)", ylabel = "")
plotGradient(m, Gradient2, pred = predY2, measure = "T", index = 2, pointsize = 0.2,pointcol="red",showData = T, xlabel = expression(MAT~(degree*C)), ylabel = "")
plotGradient(m, Gradient2_2, pred = predY2_2, measure = "T", index = 2, pointsize = 0.2,pointcol="red",showData = T, xlabel = expression("Altitude (m)"), ylabel = "")
plotGradient(m, Gradient3, pred = predY3, measure = "T", index = 2, pointsize = 0.2,pointcol="red",showData = T, xlabel = "Soil N (%)", ylabel = "")
plotGradient(m, Gradient4, pred = predY4, measure = "T", index = 2, pointsize = 0.2,pointcol="red",showData = T, xlabel = "Soil P (%)", ylabel = "")

dev.off()


plotGradient(m, Gradient1, pred = predY1, measure = "S", showData = T)
# if model fit on abundance data,measure = "S" to plot the summed abundance over all species, row sum of the predicted communites.
# if present-absent data, measure S would give the expected species richness.

# one can visualize the same prediction for individual species by setting measure = "Y" and by using index to
## select the species to be visualized (as orderd in the matrix m$Y)

plotGradient(m, Gradient, pred = predY, measure = "Y", index= 1, showData = T)

# Finally, by selecting measure = "T", we can visualize how community-weighted mean trait values behave 
# over the environmental gradient. Now the index selects the traits (as ordered in the matrix m$Tr)

plotGradient(m, Gradient1, pred = predY1, measure = "T", index = 2, showData = T)







# one can also construct an environmental gradient over the habitat types
#Gradient = constructGradient(m, focalVariable = "habitat", non.focalVariables = list("climate" = list(1)))
# Let us select the species for which the trait value for habitat use is the highest, and then plot 
# how that species respond to the habitat gradient.
#predY = predict(m, XData = Gradient$XDataNew, studyDesign = Gradient$studyDesignNew, ranLevels = Gradient$rLNew, expected = T)
#plotGradient(m, Gradient, pred = predY, mmeasure = "Y", index= which.max(m$TrData$habitata.use), showData = T, jigger = 0.2)
## i.e. this will show----- the effect of habitat on the abundance of the most habitat-resposne species
## similar, for AMF community, we might want to plot the abundance change from aT to eT for the most temperature-response species
## jigger = 0.2 to randomly move the observed data (the dots) in the horizontal direction to avoid overlap
#plotGradient(m, Gradient, pred = predY, measure = "T", index = 2, showData = T, jigger = 0.2)
## as expected, the averaged value of this trait (weighted by expoonentially transformed species abundances) is higher in forest than in open habitats.


## next step is to make figure for phylogenetic tree

par(fig = c(0, 0.6, 0, 0.8), mar = c(6, 0,2, 0))
library(phylosignal)
library(phylobase)

target_spe <- data.frame(phylo$tip.label)
colnames(target_spe) <- "Species"
library(dplyr)
TrData_new <- left_join(target_spe,TrData,by="Species")



gs_4d <- phylo4d(x=phylo, tip.data=TrData_new$GS)

png("Output/fig_phylo.png", width = 8, height = 16, units = 'in', res = 700)
barplot(gs_4d,center = FALSE,scale = F,trait.bg.col = c("#CEF6CE"), bar.col = "grey35",bar.lwd=2, show.trait = T, trait.labels = "1C DNA content (pg)")
focusTraits(1)
abline(v = 1.57, col = 2,lty=2)
abline(v = 2.5, col = 2,lty=2)
abline(v = 4.02, col = 2,lty=2)
abline(v = 5.13, col = 2,lty=2)
dev.off()


phyloSignal(p4d = gs_4d, method = "all")
