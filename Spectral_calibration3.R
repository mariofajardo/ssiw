#####LOG_PLS model####
log_soil_c <- log(soil_c)
maxc<-14# number of max components
log_soil_c.pls_model <- plsr(log_soil_c ~ spec_c, maxc, validation = "CV")
plot(log_soil_c.pls_model, "val") # RMSEP curves

nc<-10 # no. components to use

str(log_soil_c.pls_model)
plot(log_soil_c.pls_model, ncomp = nc) # Plot of cross-validated predictions
plot(soil_c.pls_model, ncomp = nc) # Plot of cross-validated predictions


# to predict value based on a given spectra (calibration data)
log_soil_c.pls_predict <- predict(log_soil_c.pls_model, ncomp = nc, newdata = spec_c)
# goodness of fit

gfc.log_pls_predict<- goof(soil_c,exp(log_soil_c.pls_predict),xlab="Observed", ylab="Predicted", main="PLSR Calibration")
gfc.log_pls_predict

# to predict value on validation data

log_soil_v.pls_predict <- predict(log_soil_c.pls_model, ncomp = nc, newdata = spec_v)
# goodness of fit

gfv.log_pls_predict<- goof(soil_v,exp(log_soil_v.pls_predict),xlab="Observed", ylab="Predicted", main="PLSR Validation")

gfc.pls_predict# Prediction statistics
gfv.pls_predict# Validation statistics

gfc.log_pls_predict# Prediction statistics
gfv.log_pls_predict# Validation statistics

#######

# log bagging PLSR
nbag<-50# no. of bootstrap
maxc<-14# max. number of components
nc<-maxc # no. of PLSR components in model
fit_bag_plsr
log_bag_plsr<-fit_bag_plsr(log_soil_c,spec_c,nbag,maxc)

# make the prediction on calibration set
log_soil_c.bpls_predict<-predict_bag_plsr(log_bag_plsr$model.bpls,spec_c,nbag,nc)
str(soil_c.bpls_predict)
# goodness of fit
gfc.log_bpls_predict<- goof(soil_c,exp(log_soil_c.bpls_predict$pred.ave),xlab="Observed", ylab="Predicted", main="BagPLSR Calibration")

# make the prediction on validation set
log_soil_v.bpls_predict<-predict_bag_plsr(log_bag_plsr$model.bpls,spec_v,nbag,nc)
# goodness of fit
gfv.log_bpls_predict<- goof(soil_v,exp(log_soil_v.bpls_predict$pred.ave),xlab="Observed", ylab="Predicted", main="BagPLSR Validation")

gfc.bpls_predict# Prediction statistics
gfv.bpls_predict# Validation statistics

gfc.log_bpls_predict# Prediction statistics
gfv.log_bpls_predict# Validation statistics

###############

# Make a Cubist model
soil_c.cubist_model<-cubist(x= spec_c, y=soil_c)     # fit cubist model
summary(soil_c.cubist_model)# summary of the model

# predict on calibration data
soil_c.cubist_predict<-predict(soil_c.cubist_model, spec_c) 
# goodness of fit
gfc.cpredict<- goof(soil_c,soil_c.cubist_predict,xlab="Observed", ylab="Predicted", main="Cubist Calibration")


# predict on validation data
soil_v.cubist_predict<-predict(soil_c.cubist_model, spec_v) 
gfv.cpredict<- goof(soil_v,soil_v.cubist_predict,xlab="Observed", ylab="Predicted", main="Cubist Validation")
gfv.cpredict
gfv.pls_predict

#the model is overfitting...so I can change some control values

soil_c.cubist_model<-cubist(x= spec_c, y=soil_c,control=cubistControl(rules=4))     # fit cubist model
summary(soil_c.cubist_model)# summary of the model

# predict on calibration data
soil_c.cubist_predict<-predict(soil_c.cubist_model, spec_c) 
# goodness of fit
gfc.cpredict<- goof(soil_c,soil_c.cubist_predict,xlab="Observed", ylab="Predicted", main="Cubist Calibration")


# predict on validation data
soil_v.cubist_predict<-predict(soil_c.cubist_model, spec_v) 
gfv.cpredict<- goof(soil_v,soil_v.cubist_predict,xlab="Observed", ylab="Predicted", main="Cubist Validation")

gfv.cpredict
gfv.pls_predict

gfv.cpredict<- goof(soil_v,soil_v.cubist_predict,xlab="Observed", ylab="Predicted", main="Cubist Validation")


# Make a plot to see the imporatant wavelengths

plot(soil_c.cubist_model$usage[,3], soil_c.cubist_model$usage[,2], type="h", col="plum", xlab="wavelength (nm)", ylab="Percent model usage")# see which variables are important
lines(soil_c.cubist_model$usage[,3], soil_c.cubist_model$usage[,1], type="h", col="blue")# see which variables are important
par(new=T)   
plot(wavelength10, spec_c[13,],axes=F,  
     ylim=c(-2,2), xlab="", ylab="", 
     type="l", main="",xlim=c(500,2450))

 
# Other models: Random Forests
# Just for demonstration
install.packages("randomforest")
library(randomForest)
# Need to rename the "columns" of spectra, randomForest doesnt like to have numerical column names
colnames(spec_c) = paste("x", colnames(spec_c), sep="")
# Build a RF model
soil_c.rf_model <- randomForest(soil_c ~ .,data=spec_c, ntree=500, mtry=10, importance=TRUE)
# Predict on to the calibration data
soil_c.rf_predict<-predict(soil_c.rf_model, spec_c) 
# goodness of fit
gfc.rf_predict<- goof(soil_c,soil_c.rf_predict)

# predict on to the validation data
colnames(spec_v) = paste("x", colnames(spec_v), sep="")
soil_v.rf_predict<-predict(soil_c.rf_model, spec_v) 
# goodness of fit
gfv.rf_predict<- goof(soil_v,soil_v.rf_predict)
# Variable of Importance Plot
varImpPlot(soil_c.rf_model) 
