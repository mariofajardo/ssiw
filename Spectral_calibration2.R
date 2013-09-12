#####PLS model####
maxc<-30# number of max components
soil_c.pls_model <- plsr(soil_c ~ spec_c, maxc, validation = "CV")
plot(soil_c.pls_model, "val") # RMSEP curves

nc<-14 # no. components to use

str(soil_c.pls_model)
plot(soil_c.pls_model, ncomp = nc) # Plot of cross-validated predictions
plot(soil_c.pls_model, "loadings", comps = 1:3) # The three first loadings
plot(soil_c.pls_model, "coef", ncomp = nc) # Coefficients

# Note: if you want to plot the wavelength on each component's loading you can use:
# plot(wavelength10,soil_c.pls_model$loadings[,1],type="l")

# to predict value based on a given spectra (calibration data)
soil_c.pls_predict <- predict(soil_c.pls_model, ncomp = nc, newdata = spec_c)
# goodness of fit

gfc.pls_predict<- goof(soil_c,soil_c.pls_predict,xlab="Observed", ylab="Predicted", main="PLSR Calibration")
gfc.pls_predict

# to predict value on validation data

soil_v.pls_predict <- predict(soil_c.pls_model, ncomp = nc, newdata = spec_v)
# goodness of fit

gfv.pls_predict<- goof(soil_v,soil_v.pls_predict,xlab="Observed", ylab="Predicted", main="PLSR Validation")

gfc.pls_predict# Prediction statistics
gfv.pls_predict# Validation statistics

#######
# Now let's try bagging PLSR
nbag<-50# no. of bootstrap
maxc<-14# max. number of components
nc<-maxc # no. of PLSR components in model
fit_bag_plsr
bag_plsr<-fit_bag_plsr(soil_c,spec_c,nbag,maxc)

# make the prediction on calibration set
soil_c.bpls_predict<-predict_bag_plsr(bag_plsr$model.bpls,spec_c,nbag,nc)
str(soil_c.bpls_predict)
# goodness of fit
gfc.bpls_predict<- goof(soil_c,soil_c.bpls_predict$pred.ave,xlab="Observed", ylab="Predicted", main="BagPLSR Calibration")

# make the prediction on validation set
soil_v.bpls_predict<-predict_bag_plsr(bag_plsr$model.bpls,spec_v,nbag,nc)
# goodness of fit
gfv.bpls_predict<- goof(soil_v,soil_v.bpls_predict$pred.ave,xlab="Observed", ylab="Predicted", main="BagPLSR Validation")

gfc.bpls_predict# Prediction statistics
gfv.bpls_predict# Validation statistics

hist(soil_c)
hist(sqrt(soil_c))
hist(log(soil_c))

