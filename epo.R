# Copyright 2013 Soil security Lab, The University of Sydney
# Example of EPO in removing moisture effect from the spectra

# LIBRARIES
library(signal)
library(Cubist)
library(pls)
library(spectroscopy)

# Set the WORKING DIRECTORY
setwd("data/")

# Import the soil library data 
soil_data<-read.csv('soil_data.txt')
spectra<-read.csv('spectra.txt')
# sort the soil data
soil_data <-  soil_data[order(soil_data[,1]) , ]
nc<-ncol(spectra)
spectra <-spectra[, 2:nc]
spectra<-log(1/spectra)# convert reflectance to absorbance
spec_filtered<-filter_sg(spectra, n = 11, p = 2, m = 0)
wavelength<-seq(350,2500,by =1) # The entire range spectrum wavelengths and resolution
# trim spectra
spec_trim10<-strip_spectra(spec_filtered, wavelength,
                           wavlimits = range(500:2450), which = 10)
# Do a Standard normal variate transformation
spec_snv<- snvBLC(spec_trim10)

# In this example from Minasny et al. (2010), we have 100 soil samples (a subset of the larger dataset) 
# under 3 different moisture conditions. So we have 3 sets of spectra data with different moisture content.
# import the data in csv format
spectra0<-read.csv('moisture_dry.txt',header=FALSE)
spectra1<-read.csv('moisture_wet.txt',header=FALSE)
spectra2<-read.csv('moisture_wet2.txt',header=FALSE)
# import the corresponding soil C data
soilC<-read.csv('C.txt')
# SAVITSKY-GOLAY SMOOTHING FILTER 
abs_filtered0<-filter_sg(spectra0, n = 11, p = 2, m = 0)
abs_filtered1<-filter_sg(spectra1, n = 11, p = 2, m = 0)
abs_filtered2<-filter_sg(spectra2, n = 11, p = 2, m = 0)

# TRIM THE SPECTRA
spectra_trim0<-strip_spectra(abs_filtered0, wavelength,
                             wavlimits = range(500:2450), which = 10)
spectra_trim1<-strip_spectra(abs_filtered1, wavelength,
                             wavlimits = range(500:2450), which = 10)
spectra_trim2<-strip_spectra(abs_filtered2, wavelength,
                             wavlimits = range(500:2450), which = 10)
wavelength10<-seq(500,2450,by =10)

# plot the absorbance spectra
plot(wavelength10,spectra_trim0[1,],type="l",ylim=c(0,1.5))
lines(wavelength10,spectra_trim1[1,],col="blue")
lines(wavelength10,spectra_trim2[1,],col="green")

# perform snv
spec_snvC0<- snvBLC(spectra_trim0)
spec_snvC1<- snvBLC(spectra_trim1)
spec_snvC2<- snvBLC(spectra_trim2)

# plot the absorbance-snv spectra
plot(wavelength10,spec_snvC0[1,],type="l")
lines(wavelength10,spec_snvC1[1,],col="blue")
lines(wavelength10,spec_snvC2[1,],col="green")

####
# Make a Cubist model from original dry spectra
# select soil variables to predict
soilv<- soil_data$Total_Carbon  # change this to the soil variable you want to predict
isrow<-complete.cases(soilv)# find which row has complete record (no NA values)
specv<- spec_snv[isrow,] #specify type of spectra to use for prediction
soilv<-soilv[isrow]

# Generate a Cubist model
soilv.cubist_model<-cubist(x= specv, y=soilv)
gf.cubist_predict<- goof(soilv,predict(soilv.cubist_model,specv))
gf.cubist_predict

# Predict the values from spectra at different moisture content
soilv.cubist_predict.dry<-predict(soilv.cubist_model, spec_snvC0) 
soilv.cubist_predict.wet<-predict(soilv.cubist_model, spec_snvC1) 
soilv.cubist_predict.wet2<-predict(soilv.cubist_model, spec_snvC2) 

gfdry.cubist_predict<- goof(soilC$TotalC,soilv.cubist_predict.dry)
gfwet.cubist_predict<- goof(soilC$TotalC,soilv.cubist_predict.wet)
gfwet2.cubist_predict<- goof(soilC$TotalC,soilv.cubist_predict.wet2)

gfdry.cubist_predict
gfwet.cubist_predict
gfwet2.cubist_predict

plot(soilC$TotalC,soilv.cubist_predict.dry,xlim=c(0,10), ylim=c(0,10),xlab="Observed",ylab="Predicted")
points(soilC$TotalC, soilv.cubist_predict.wet, col="red")
points(soilC$TotalC, soilv.cubist_predict.wet2,col="green")
abline(a = 0, b = 1, col = "brown4")


## Now let's calculate the EPO
# D is the difference matrix (between dry and wet spectra)
D=as.matrix(spec_snvC0-spec_snvC1)
npc<-4  # define no. EPO factors
P<- epo(D,npc)
# P is the projection matrix
source('myImagePlot.R')# we use myImagePlot function to visualize the P matrix
myImagePlot(P,  zlim=c(-0.1,0.1), title=c("Projection matrix"))


# Now project the spectra
Z0 <- as.matrix(spec_snvC0) %*% P    # EPO projected spectra of spec0
Z1 <- as.matrix(spec_snvC1) %*% P    # EPO projected spectra of spec1
Z2 <- as.matrix(spec_snvC2) %*% P    # EPO projected spectra of spec2

# plot  the epo transformed spectra of sample 1
plot(wavelength10,Z0[1,],"l")
lines(wavelength10,Z1[1,],"l",col="blue")
lines(wavelength10,Z2[1,],"l",col="green")


## Now project the library spectra into transformed spectra
specZ <- as.matrix(specv) %*% P
specZ <- as.data.frame(specZ)


# Make a Cubist model from EPO transformed dry spectra
epo.cubist_model<-cubist(x= specZ, y=soilv)
gf.epo_cubist<- goof(soilv,predict(epo.cubist_model,specZ))
gf.epo_cubist

# Predict the values from spectra at different moisture content
epo.cubist_predict.dry<-predict(epo.cubist_model, as.data.frame(Z0))
epo.cubist_predict.wet<-predict(epo.cubist_model, as.data.frame(Z1))
epo.cubist_predict.wet2<-predict(epo.cubist_model, as.data.frame(Z2))

gfdry.epo_cubist<- goof(soilC$TotalC,epo.cubist_predict.dry)
gfwet.epo_cubist<- goof(soilC$TotalC,epo.cubist_predict.wet)
gfwet2.epo_cubist<- goof(soilC$TotalC,epo.cubist_predict.wet2)

gfdry.epo_cubist
gfwet.epo_cubist
gfwet2.epo_cubist

plot(soilC$TotalC,epo.cubist_predict.dry,xlim=c(0,10), ylim=c(0,10),xlab="Observed",ylab="Predicted")
points(soilC$TotalC, epo.cubist_predict.wet, col="red")
points(soilC$TotalC, epo.cubist_predict.wet2,col="green")
abline(a = 0, b = 1, col = "brown4")

# Plot dry vs wet predicted
plot(epo.cubist_predict.dry, epo.cubist_predict.wet)
abline(a = 0, b = 1, col = "brown4")


