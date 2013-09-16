##### Example routine for calibration model fitting, and prediction of soil properties using from NIR spectra####

# LIBRARIES
library(signal)
library(Cubist)
library(pls)
library(spectroscopy)
#load('correct_steps.RData')


# Set the WORKING DIRECTORY
#setwd("data/")

# import the data in csv format
soil_data<-read.csv('soil_data.txt')
spectra<-read.csv('spectra.txt')

# sort the soil data
soil_data <-  soil_data[order(soil_data[,1]) , ]
nc<-ncol(spectra)
spectra <-spectra[, 2:nc]


wavelength<-seq(350,2500,by =1) # The entire range spectrum wavelengths and resolution

cont_spectra<-log(1/spectra)# convert reflectance to absorbance

#cont_spectra <-correct_step(cont_spectra)

spec_filtered<-filter_sg(cont_spectra, n = 11, p = 2, m = 0)# smooth the spectra

# trim spectra
spec_trim10<-strip_spectra(spec_filtered, wavelength,
                           wavlimits = range(500:2450), which = 10)
wavelength10<-seq(500,2450,by =10)

# Do a BASELINE CORRECTION for absorbance:  Standard normal variate transformation
spec_snv<- snvBLC(spec_trim10)



### Prediction
# select soil variables to predict
soily<- soil_data$Total_Carbon  # change this to the soil variable you want to predict
isrow<-complete.cases(soily)# find which row has complete record (no NA values)
spy<- spec_snv[isrow,] #specify type of spectra to use for prediction
soily<-soily[isrow]
nd<-length(soily)

# Split Calibration: 75% and Validation 25%
set.seed(111)# set a randoms seed (so you don't get different results every time)
ic<-sample(1:nd, round(nd*0.75))# generate a random permutation of data
#calibration set  
spec_c<-spy[ic,]
soil_c<-soily[ic]
#validation set
spec_v<-spy[-ic,]
soil_v<-soily[-ic]


# First we do a principal component analysis
pc_spec<-prcomp(spec_c,center=TRUE,scale=TRUE)
summary(pc_spec)
screeplot(pc_spec, type="lines") # plot variances explained by each component
v<- (pc_spec$sdev)^2 # calculate variances
cumv<-100*cumsum(v)/sum(v)  # percentage of cumulative variances
plot(cumv[1:20], type="b", xlab="PC", ylab="% Cumulative variances")
abline(h=99)
# a biplot
biplot(pc_spec)
# score plot
scoreplot(pc_spec,comps=1:2)
# loading plot
loadingplot(pc_spec,comps=1:2)
loadingplot(pc_spec,comps=1)
loadingplot(pc_spec,comps=2)

# Let's do a simple PC regression
npc<-9# no. of components that explain (roughly) 99% of the data
# Fit a linear model Total C = PC1 + PC2 + ...
sdata<-as.data.frame(pc_spec$x[,1:npc]) # x are the principal components scores
soil_c.pcr_model<-lm(soil_c ~ ., data=sdata) #the 'dot' is for telling that 'all' the varioables will be used in the model
summary(soil_c.pcr_model)
# goodness of fit
soil_c.pcr_predict<-predict(soil_c.pcr_model,sdata)
gfc.pcr_predict<- goof(soil_c,soil_c.pcr_predict, xlab= "Observed", ylab="Predicted", main="PCR Calibration")

# Let's predict to validation data
pc_scores_v<-predict(pc_spec,spec_v) #this gives the pc scores of the v_data ???????
#pc_test <-prcomp(spec_v,center=T,scale=T)
#identical(pc_scores_v,pc_test)

newdata<-as.data.frame(pc_scores_v[,1:npc])
soil_v.pcr_predict<-predict(soil_c.pcr_model,newdata)
gfv.pcr_predict<- goof(soil_v,soil_v.pcr_predict, xlab= "Observed", ylab="Predicted", main="PCR Validation")

gfc.pcr_predict
gfv.pcr_predict

######

