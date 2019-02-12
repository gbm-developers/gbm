

################################################################################
#############################                 ##################################
#############################   DEMO SPECIES  ##################################
#############################                 ##################################
################################################################################

{
cat("This is a demonstration of a basic exploratory anyalysis of a 
Presence/Absence dataset. We will create tables of accuracy measures, 
make graphs, and examine various methods of optimizing our choice of 
threshold."
,
fill=TRUE,sep="")

if (interactive()) {
	cat("\nType  <Return>\t to continue : ")
	readline()}
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

# Next go to the command window and simply type:

library(PresenceAbsence)# Next go to the command window and simply type:

	library(PresenceAbsence) 	

# To load the dataset for this demo, type:

	data(SPDATA)

# For more details on this dataset see:
#
# Moisen, G.G., Freeman, E.A., Blackard, J.A., Frescino, T.S., 
# Zimmerman N.E., Edwards, T.C. Predicting tree species presence and 
# basal area in Utah: A comparison of stochastic gradient boosting, 
# generalized additive models, and tree-based methods. Ecological 
# Modellng, 199 (2006) 176-187.   
 

################## take a look at the data format ###########################

# Let's take a look at the first few rows of the datasets:

	SPDATA[1:10,]

# First look at the data structure:
#
# Column 1 is the plot ID's. In this case, The ID column is made up the 
# species code of each observation. 
#
# Column 2 is the observed values where '0 = Absent' and '1 = Present'. If 
# the observed values are actual values (for example, basal area or tree 
# counts) any values other than zero will be treated as 'Present'.
#
# The remaining columns are for model predictions. Preferably these will
# be predicted probabilities ranging from zero to one. If all that is 
# available is predicted Presence/Absence values, basic accuracy statistics
# can still be calculated, but you loose the ability to examine the effect
# of threshold choice, and most of the graphs are not meaningful.
#
# In this case, they are the model predictions from three types of models:
# GAM - generalized additive models, See5 - Rulequests software package, and
# SGB - Stochastic gradient boosting.
#
# Note also: the functions rely on column order to identify the Observed and
# Predicted values, not on the column names. Therefore, as long as the columns 
# are in the correct order, column names can be anything you choose. The only
# time the functions use the names is for the default labels on graphs.

##############################################################################

# Lets take a look at the species codes

	unique(SPDATA$SPECIES)
	
# This has generated a vector of the species codes present in SPDATA.
# This vector will be useful later when we want to loop through all
# the species.

# Next, lets check how many observations there are for each species:

	table(SPDATA$SPECIES)
	table(SPDATA$SPECIES,SPDATA$OBSERVED)

# There are 386 observations for each species, but the observed prevalence 
# varies widely from species to species.

###################### Define a few variables ############################

# We will define a few variables to describe this dataset for later use:

	species<-as.character(unique(SPDATA$SPECIES))
	model.names<-as.character(names(SPDATA)[-c(1,2)])

	N.models<-ncol(SPDATA)-2
	N.sp<-length(species)
	N.obs<-386

################### Presence/Absence Histograms ##########################

# Let's get a visual feel for the data. The 'presence.absence.hist' function
# is a bar plot of observed values as a function predicted probability:

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(i in 1:N.sp){
		for(mod in 1:N.models){
			presence.absence.hist(	DATA=SPDATA[SPDATA$SPECIES==species[i],],
							which.model=mod,
							N.bars=20)}
		mtext(species[i],side=3,line=1,cex=1.5,outer=T)}

# Note - use PgUp and PgDn to scroll thru species
#
# For many of the species, the zero bar is so much taller than the 
# other bars that it is difficult to see any other details of the
# species distributions. Setting the argument 'truncate.tallest' = 'TRUE'
# will truncate the tallest bar to fit better in the plot region.
# While we are at it, we will use the 'legend.cex' argument to make the
# legend more readable:

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(i in 1:N.sp){
		for(mod in 1:N.models){
			presence.absence.hist(	DATA=SPDATA[SPDATA$SPECIES==species[i],],
							which.model=mod,
							N.bars=20,
							truncate.tallest=TRUE,
							legend.cex=1.0)}
		mtext(species[i],side=3,line=1,cex=1.5,outer=T)}

# These plots are valuable in two ways:
#
# First, these plots make prevalence dramatically obvious. Many of 
# the traditional accuracy measures (such as percent correctly classified (PCC),
# sensitivity and specificity) are highly dependant on species prevalence.
# Therefore, it is important to check prevalence early in an investigation.
#
# (Note though, that the 'truncate.tallest' option can visually obscure the true 
# prevalence.)
#
#Compare the graphs of 3 species with differing prevalence:
 
#truncate.tallest=FALSE:

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(sp in c("ACGR3","PIEN","POTR5")){
		for(mod in 2){
			presence.absence.hist(	DATA=SPDATA[SPDATA$SPECIES==sp,],
							which.model=mod,
							N.bars=20,
							truncate.tallest=FALSE,
							legend.cex=1.0,
							main=sp)}
		mtext(model.names[mod],side=3,line=1,cex=1.5,outer=T)}

#truncate.tallest=True:

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(sp in c("ACGR3","PIEN","POTR5")){
		for(mod in 2){
			presence.absence.hist(	DATA=SPDATA[SPDATA$SPECIES==sp,],
							which.model=mod,
							N.bars=20,
							truncate.tallest=TRUE,
							legend.cex=1.0,
							main=sp)}
		mtext(model.names[mod],side=3,line=1,cex=1.5,outer=T)}

# Second, these plots are a good visual indicator of model quality. 
#
# Compare the histograms for the See5 models for three species with 
# similar prevalence, but increasing model quality:

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(sp in c("JUSC2","ABCO","PICO")){
		for(mod in 2){
			presence.absence.hist(	DATA=SPDATA[SPDATA$SPECIES==sp,],
							which.model=mod,
							N.bars=20,
							truncate.tallest=TRUE,
							legend.cex=1.0,
							main=sp)}
		mtext(model.names[mod],side=3,line=1,cex=1.5,outer=T)}

# A good model (such as the See5 model for PICO) will have a distinctly 
# double humped histogram. If you choose a threshold in the valley 
# (not necessarily at 0.5!) it will divide the data cleanly into 
# Observed Present and Observed Absent.
#
# A poor model (such as the See5 model for JUSC2) will have considerable  
# overlap between Observed Present and Observed Absent. No threshold exists 
# that will cleanly divide the data into Observe Present and Observed Absent. 
# The challenge here is to find a threshold that offers the best compromise. 
#
# There are several methods available to optimize the choice of threshold.
# Which method is best often depends on the purpose of your model.


#################### Predicted Prevalence Tables ############################

# Let's start by examining prevalence a little more closely...
#
# If it is important that the model predicts the prevalence of a species
# correctly, you can choose a threshold where the predicted threshold
# equals the observed threshold.
#
# 'predicted.prevalence' let's us make a table of the effect of threshold
# choice on the predicted prevalence of our models:

# Lets start with a naive approach, using the standard cutoff of 
# threshold = 0.5, for declaring a species present or absent:

	pred.prev<-predicted.prevalence(DATA=SPDATA[SPDATA$SPECIES==species[1],])
	for(sp in 2:N.sp){
		pred.prev<-rbind(pred.prev,predicted.prevalence(DATA=SPDATA[SPDATA$SPECIES==species[sp],]))}
	pred.prev<-cbind(species=species,pred.prev)
	pred.prev

# We can look even closer at individual species to examine the effect of threshold
# choice on prevalence:

	sp<-1
	print(paste("Species:",species[sp]))
	predicted.prevalence(DATA=SPDATA[SPDATA$SPECIES==species[sp],],threshold=11)

# Setting 'threshold = 11' causes the function to calculate 11 evenly spaced
# thresholds between zero and one. If you are interested in particular 
# thresholds you can, for example, set 'threshold = c(.23,.6,.97)'.
#
# If it is important that the model predicts the prevalence of a species
# correctly, you can choose a threshold where the predicted threshold
# equals the observed threshold.

###################### Accuracy Tables #######################################

# We will make some tables of the accuracy statistics.
#
# First we will look at all three models at the default threshold of 0.5:

	sp<-1
	DATA<-SPDATA[SPDATA$SPECIES==species[sp],]
	print(paste("Species:",species[sp]))
	presence.absence.accuracy(DATA)

# Note: The standard deviations in this table are only appropriate for
# examining each model individually. To do significance testing for differences
# between models it is necessary to take into account correlation. See
# 'help(auc)' for further details.
#
# We can suppress the standard deviations by setting 'st.dev=FALSE'
#
# The default for 'presence.absence.accuracy' is to calculate accuracy measures
# for all models at 'threshold = 0.5'. We can also substitute a different 
# threshold:

	sp<-1
	DATA<-SPDATA[SPDATA$SPECIES==species[sp],]
	print(paste("Species:",species[sp]))
	presence.absence.accuracy(DATA,threshold=0.4,st.dev=FALSE)

# When the table is being produced for multiple models, 'threshold' must be
# either of length 1, or it must be the same length as the number of models
# in the table. In the later case, the first threshold will be used for the 
# first model, the second threshold for the second row, etc...

	sp<-1
	DATA<-SPDATA[SPDATA$SPECIES==species[sp],]
	print(paste("Species:",species[sp]))
	presence.absence.accuracy(DATA, which.model=c(1,3),threshold=c(0.4,.5),st.dev=FALSE)

# 'presence.absence.accuracy' can also be used to produce a table for a single 
# model, illustrating how the accuracy measures change as a function of 
# threshold:

	sp<-1
	DATA<-SPDATA[SPDATA$SPECIES==species[sp],]
	print(paste("Species:",species[sp]))
	presence.absence.accuracy(DATA, which.model=3,threshold=11,st.dev=FALSE)
	presence.absence.accuracy(DATA, which.model=3,threshold=c(.2,.4,.5,.89),st.dev=FALSE)

# You will notice on these tables, that the area under the curve (AUC) is
# constant for all thresholds. The AUC is the only one of these basic accuracy
# statistics that is threshold independent.
#
# Calculating the AUC for large datasets can be very memory intensive. 
# Setting the argument 'find.auc' = 'FLASE' will turn off the AUC calculations:

	sp<-1
	DATA<-SPDATA[SPDATA$SPECIES==species[sp],]
	print(paste("Species:",species[sp]))
	presence.absence.accuracy(DATA,find.auc=FALSE)

################# Graphing the basic Accuracy statistics #####################

# If you need to look up specific values, tables are necessary, but often it
# is easier to spot patterns in a graph. Lets take a graphical look at how
# these error statistics vary with threshold:

# Selected Species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(sp in c("JUSC2","PICO")){
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			error.threshold.plot(DATA, which.model=mod, color=TRUE,legend.cex=1.1)}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)}

# All Species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			error.threshold.plot(DATA, which.model=mod, color=TRUE,legend.cex=1.1)}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev$Obs.Prevalence[i],2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)} 


# These plots always include sensitivity and specificity. Some values of the 
# argument 'opt.methods' will result in adding additional accuracy statistics to  
# the graphs. These lines for the error statistics specified in 'opt.methods'
# will be plotted, even if 'opt.thresholds=FALSE'.
#
# Values for 'opt.methods' that add lines to 'error.threshold.plot' include:
#
# 	"MaxKappa"		adds Kappa
# 	"MaxPCC"		adds PCC (Percent correctly classified)
#	"MinROCdist"	adds the distance from the ROC plot to the upper left corner
#	"MaxSens+Spec"	adds the mean of Sensitivity and Specificity: (Sensitivity+Specificity)/2
# 
# 
# The default is 'opt.methods=1:4' which results in sensitivity, specificity and Kappa.
#
# The argument 'vert.lines' can be used to specify how the optimized thresholds are marked
# on the plot:
#	'vert.lines=TRUE'	a vertical line is drawn at each optimized threshold, and labeled along
#				the top of the plot.
#	'vert.lines=FLASE'the thresholds are marked along the error statistics (the default)
#
# A further example:

# Selected species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(sp in c("JUSC2","PICO")){
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			error.threshold.plot(DATA,which.model=mod,color=TRUE,opt.methods="MaxPCC",vert.lines=TRUE,legend.cex=1.1,opt.legend.cex=1.1)}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)}

# All species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			error.threshold.plot(DATA,which.model=mod,color=TRUE,opt.methods="MaxPCC",vert.lines=TRUE,legend.cex=1.1,opt.legend.cex=1.1)}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev$Obs.Prevalence[i],2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)} 


################### Optimizing Thresholds ##################################

# As mentioned above, the PresenceAbsence library provides multiple methods for optimizing 
# thresholds. The full list includes:
#
# opt.methods	optimization methods:
#	1  "Default"	threshold=0.5
#	2  "Sens=Spec"	threshold where sensitivity equals specificity
#	3  "MaxSens+Spec"	threshold that maximizes the sum of sensitivity and specificity
#	4  "MaxKappa"	threshold that maximizes Kappa
#	5  "MaxPCC"		threshold that maximizes PCC (percent correctly classified)
#	6  "PredPrev=Obs"	threshold where predicted prevalence equals observed prevalence
#	7  "ObsPrev"	threshold set to Observed prevalence
#	8  "MeanProb"	threshold set to mean predicted probability
#	9  "MinROCdist"	threshold where ROC plot makes closest approach to (0,1)
#	10 "ReqSens"	lowest threshold where sensitivity meets user defined requirement
#	11 "ReqSpec"	lowest threshold where specificity meets user defined requirement
#	12 "Cost"		threshold that meets user defined relative costs ratio
#
# The 'opt.methods' argument is specified by a vector of choices from this list. 
# The methods can be given by number ( 'opt.methods=1:10' or 'opt.methods=c(1,2,4)') 
# or by name ('opt.methods=c("Default","Sens=Spec","MaxKappa")').
#
# "Default" - First, the default method of setting 'threshold = 0.5'
#
# "Sens=Spec" - The second method for optimizing threshold choice is by finding the  
# thresholdwhere sensitivity equals specificity. In other words, find the threshold
# where positive observations are just as likely to be wrong as negative
# observations. 
#
# Note: when threshold is optimized by this method it is highly 
# correlated to prevalence, so that rare species are given much lower
# thresholds than widespread species. As a result, rare species may give the
# appearance of inflated distribution, if maps are made with thresholds that
# have been optimized by this method.
#
# "MaxSens+Spec" The third method chooses the threshold that maximizes the sum 
# of sensitivity and specificity. In other words, it is minimizing the mean of 
# the error rate for positive observations and the error rate for negative observations.
# This method is also highly correlated to prevalence. 
#
# "UserAllow" - The third method chooses the threshold that maximizes the sum of  
# sensitivity and specificity. In other words, it is minimizing the mean of the 
# positive error rate and the negative error rate. This method of optimization is   
# also highly correlated to prevalence.
#
# "MaxKappa" - The forth method for optimizing the threshold choice is to find the
# threshold that gives the maximum value of Kappa. Kappa makes full use of the 
# information in the confusion matrix to asses the improvement over chance 
# prediction.
#
# "MaxPCC" - The fifth method is to maximize the total accuracy (PCC - Percent   
# Correctly Classified).
#
# Note: It may seem like maximizing total accuracy would be the obvious goal, 
# however, there are many problems with using PCC to assess model accuracy. 
# For example, with species with very low prevalence, it is possible to maximize 
# PCC simply by declaring the species a absent at all locations -- not a very 
# useful prediction!
#
# "PredPrev=Obs" - The sixth method is to find the threshold where the Predicted  
# prevalence is equal to the Observed prevalence. This is a useful method when 
# preserving prevalence is of prime importance. If you have observed prevalence data 
# from another source, you can use the argument 'obs.prev' to set the observed
# prevalence. Otherwise, 'obs.prev' defaults to the observed prevalence from 'DATA'.
#
# "ObsPrev" - The seventh method is an even simpler variation, where you simply set  
# the threshold to the Observed prevalence. It is nearly as good at preserving 
# prevalence as method 6 and requires no computation. If you have observed prevalence data 
# from another source, you can use the argument 'obs.prev' to set the observed
# prevalence. Otherwise, 'obs.prev' defaults to the observed prevalence from 'DATA'.
#
# "MeanProb" - The eighth method sets the threshold to the mean probability of 
# occurence from the model results.
#
# "MinROCdist" - The ninth  method is to find the threshold that minimizes the 
# distance between the ROC plot and the upper left corner of the unit square.
#
# "ReqSens" - The tenth method allows the user to set a required sensitivity, and then  
# finds the highest threshold (i.e. with the highest possible specificity) that still  
# meets the requiremed sensitivity. In other words, the user can decide that the model  
# must miss no more than, for example 15% of the plots where the species is observed to 
# be present. Therefore they require a sensitivity of at least 0.85. This will then 
# find the threshold that has the highest possible specificiy while still meeting
# the required sensitivity. This may be useful if, for example, the goal is to define 
# a management area for a rare species, and one wants to be certain that the 
# management area doesn't leave populations unprotected.
#
# "ReqSpec" - The eleventh method allows the user to set a required specificity, and  
# then findsthe lowest threshold (i.e. with the highest possible sensitivity) that still  
# meets the requiremed specificity. In other words, the user can decide that the model 
# must miss no more than, for example 15% of the plots where the species is observed 
# to be absent. Therefore they require a specificity of at least 0.85. This will then 
# find the threshold that has the highest possible sensitivity while still meeting
# the required specificity. This may be useful if, for example, the goal is to determine 
# if a species is threatened, and one wants to be certain not to over inflate the 
# population by over declaring true absences as predicted presences.
#
# Note: for "ReqSens" and "ReqSpec", if your model is poor, and your requirement is 
# too strict, it is possible that the only way to meet it will be by declaring every 
# single plot to be Present (for ReqSens) or Absent (for ReqSpec) -- not a very
# useful method of prediction!
#
# "Cost" - The twelth method balances the relative costs of false positive 
# predictions and false negative predictions. A slope is calculated as
# (FPC/FNC)((1 - prevalence)/prevalence). To determine the threshold, a 
# line of this slope is moved from the top left of the ROC plot, till it 
# first touches the ROC curve.

# Note also: 'error.threshold.plot' is a rough and ready applied function. 
# It optimizes thresholds simply by calculating a large number of evenly spaced 
# thresholds and looking for the best ones. This is good enough for graphs,
# but to find the theoretically 'best' thresholds, would require calculating
# every possible unique threshold (not necessarily evenly spaced!).


############################################################################
# The plotting sub-function 'optimal.thresholds' can also be used dirrectly
# to produce tables of optimized thresholds:

	sp<-"QUGA"
	DATA<-SPDATA[SPDATA$SPECIES==sp,]
	print(paste("Species:",sp))
	optimal.thresholds(DATA,opt.methods=1:12,req.sens=0.85,req.spec=0.85,FPC=2,FNC=1)

# The error graphs can be used in several ways to optimize threshold choice, 
# and in fact, the function 'error.threshold.plot' can calculate these 
# optimized thresholds as well as add them to the plots:

# Selected Species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3),mar=c(3,3,2,1))
	for(sp in c("JUSC2","PICO")){
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:N.models){
			error.threshold.plot(	DATA, 
							which.model=mod, 
							color=TRUE,
							opt.thresholds=TRUE,
							opt.methods=c(1,2,10,4),
							req.sens=0.85,
							legend.cex=1.0,
							opt.legend.cex=1.0,
							vert.lines=FALSE)}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)}

# All Species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3),mar=c(3,3,2,1))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:N.models){
			error.threshold.plot(	DATA, 
							which.model=mod, 
							color=TRUE,
							opt.thresholds=TRUE,
							opt.methods=c(1,2,10,4),
							req.sens=0.85,
							legend.cex=1.0,
							opt.legend.cex=1.0,
							vert.lines=FALSE)}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev$Obs.Prevalence[i],2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)} 

# To produce a table of accuracy statistics for the optimized thresholds:

	sp<-"JUSC2"
	which.model=1
	print(paste("Species:",sp))
	print(paste("Model:",model.names[which.model]))
	error.threshold.plot(	DATA=SPDATA[SPDATA$SPECIES==sp,], 
					which.model=which.model,
					opt.thresholds=TRUE,
					opt.methods=c(1,2,10,4),
					req.sens=0.85,
					plot.it=FALSE)

# You can also add the optimized thresholds to the histogram plots:

# Selected species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,3,0),mfrow=c(2,3))
	for(sp in c("JUSC2","PICO")){
		for(mod in 1:N.models){
			error.threshold.plot(	DATA=SPDATA[SPDATA$SPECIES==sp,], 
							which.model=mod, 
							color=TRUE,
							opt.thresholds=TRUE,
							opt.methods=c(1,2,10,4),
							req.sens=0.85,
							legend.cex=.9,
							add.opt.legend=FALSE)}
		for(mod in 1:N.models){
			presence.absence.hist(	DATA=SPDATA[SPDATA$SPECIES==sp,],
							which.model=mod,
							N.bars=20,
							truncate.tallest=TRUE,
							opt.thresholds=TRUE,
							opt.methods=c(1,2,10,4),
							req.sens=0.85,
							legend.cex=.9,
							opt.legend.cex=.9)}
		mtext(sp,side=3,line=0,cex=1.5,outer=T)}

# All species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,3,0),mfrow=c(2,3))
	for(i in 1:N.sp){
		sp<-species[i]
		for(mod in 1:N.models){
			error.threshold.plot(	DATA=SPDATA[SPDATA$SPECIES==sp,], 
							which.model=mod, 
							color=TRUE,
							opt.thresholds=TRUE,
							opt.methods=c(1,2,10,4),
							req.sens=0.85,
							legend.cex=.9,
							add.opt.legend=FALSE)}
		for(mod in 1:N.models){
			presence.absence.hist(	DATA=SPDATA[SPDATA$SPECIES==sp,],
							which.model=mod,
							N.bars=20,
							truncate.tallest=TRUE,
							opt.thresholds=TRUE,
							opt.methods=c(1,2,10,4),
							req.sens=0.85,
							legend.cex=.9,
							opt.legend.cex=.9)}
		mtext(sp,side=3,line=0,cex=1.5,outer=T)}

# In the above plots, the thresholds are optimized by 4 methods:
#
# First, the default method of setting 'threshold = 0.5'
#
# Next, thresholds were optimized by finding the threshold 
# where sensitivity equals specificity. In other words, find the threshold
# where positive observations are just as likely to be wrong as negative
# observations. Notice that when threshold is optimized by this method it is
# highly correlated to prevalence, so that rare species are given much lower
# thresholds than widespread species. As a result, rare species may give the
# appearance of inflated distribution, if maps are made with thresholds that
# have been optimized by this method.
#
# Third, thresholds were optimized by forcing the predictions to meet a user 
# defind required sensitivity. In other words, the user can decide that the 
# model must miss no more than, for example 15% of the plots where the species 
# is observed to be present. Therefore they reqire a sensitivity of at least 0.85. 
# This may be useful if, for example, the goal is to define a management area for 
# an endagered species, and they want to be certain that the management area 
# doesn't leave unprotected too many populations.
#
# Note: The default is for 'RecSens = 0.85'. This should be adjusted to match your
# particular management requirements, but make sure the results are reasonable.  
# If your model is poor, and your requirement is too strict, it is possible that the 
# only way to meet it will be by declaring every single plot to be Present  
# -- not a very useful method of prediction!
#
# This forth method for optimizing the threshold choice is to find the
# threshold that gives the maximum value of Kappa. Kappa makes full use of the 
# information in the confusion matrix to asses the improvement over chance 
# prediction.
#
# Note: It may seem like maximizing total accuracy (PCC - Percent Correctly  
# Classified) would be the obvious goal, however, there are many 
# problems with using PCC to asses model accuracy. For example, with species 
# with very low prevalence, it is possible to maximize PCC simply by declaring 
# the species a absent at all locations -- not a very useful prediction!
#
# Note also: 'error.threshold.plot' is a rough and ready applied function. 
# It optimizes thresholds simply by calculating a large number of evenly spaced 
# thresholds and looking for the best ones. This is good enough for graphs,
# but to find the theoretically 'best' thresholds, would require calculating
# every possible unique threshold (not necessarily evenly spaced!).
#
################################# ROC plots #############################################

# ROC plots, with their associated AUC (Area Under the Curve) provide a
# threshold independent method of assessing model performance. 

# Selected species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0))
	for(sp in c("JUSC2","PICO")){
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		auc.roc.plot(	DATA, 
					color=TRUE,
					legend.cex=.9,
					main=paste(sp,"(",signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence )"))}

# All species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		auc.roc.plot(	DATA, 
					color=TRUE,
					legend.cex=.9,
					main=paste(sp,"(",signif(pred.prev$Obs.Prevalence[i],2),"Prevalence )"))}

# ROC plots are produced by assessing the ratio of true positives to false 
# positives for all possible thresholds. The plot from a good model will rise
# steeply to the upper left corner then level off quickly, resulting in an AUC
# near 1.0. A poor model (i.e. a model that is no better than random assignment)
# will have a ROC plot lying along the diagonal, with an AUC near 0.5.
#
# The AUC is equivalent to the chance that a randomly chosen Present plot will
# have a predicted probability higher than that of a randomly chosen Absent 
# plot.
#
# The 'mark' argument will mark particular thresholds along each models ROC
# plot. Notice that evenly spaced thresholds do not end up even distances apart
# along a ROC plot.

# Single species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0))
	sp="QUGA"
	DATA=SPDATA[SPDATA$SPECIES==sp,]
	auc.roc.plot(	DATA, 
				color=TRUE,
				legend.cex=.9,
				mark=0.5,
				main=paste(sp,"(",signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence )"))

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0))
	sp="QUGA"
	DATA=SPDATA[SPDATA$SPECIES==sp,]
	auc.roc.plot(	DATA, 
				color=TRUE,
				legend.cex=.9,
				mark=5,
				pch=c(1,16),
				main=paste(sp,"(",signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence )"))

# All species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		auc.roc.plot(	DATA, 
					color=TRUE,
					legend.cex=.9,
					mark=0.5,
					main=paste(sp,"(",signif(pred.prev$Obs.Prevalence[i],2),"Prevalence )"))}

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0))
	for(sp in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		auc.roc.plot(	DATA, 
					color=TRUE,
					legend.cex=.9,
					mark=5,
					pch=c(1,16),
					main=paste(sp,"(",signif(pred.prev$Obs.Prevalence[i],2),"Prevalence )"))}


# Even species with the same prevalence can have very different model quality:

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(sp in c("JUSC2","ABCO","PICO")){
	DATA=SPDATA[SPDATA$SPECIES==sp,]
	auc.roc.plot(	DATA,  
				color=TRUE,
				legend.cex=1.1,
				main=paste(	sp,"(",round(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence )"))} 

# The optimal thresholds from 'error.threshold.plot' can be added to the ROC
# plots. This can get hard to read if all the models are on the same plot, so
# we will split the models onto separate plots:

# Selected species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3),mar=c(3,4,2,1))
	for(sp in c("JUSC2","PICO")){
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			auc.roc.plot(	DATA, 
						which.model=mod, 
						color=mod+1,
						opt.thresholds=TRUE,
						opt.methods=c(1,4,6),
						mark.numbers=TRUE,
						main=model.names[mod],
						legend.cex=1.1,
						opt.legend.cex=1.1)}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)}


# All species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3),mar=c(3,4,2,1))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			auc.roc.plot(	DATA, 
						which.model=mod, 
						color=mod+1,
						opt.thresholds=TRUE,
						opt.methods=c(1,4,6),
						mark.numbers=TRUE,
						main=model.names[mod],
						legend.cex=1.1,
						opt.legend.cex=1.1)}
	mtext(sp,side=3,line=1,cex=1.5,outer=T)
	mtext(paste(signif(pred.prev$Obs.Prevalence[i],2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)} 

################################### Calibration Plots #########################################

# Calibration plots provide a goodness-of-fit plot for Presence/Absence models.
#
# The plots are grouped into bins based on their predicted values, and then
# observed bin prevalence (the ratio of plots in each bin with observed values 
# of present verses the total number of plots in each bin) is calculated for each bin. 
# The confidence interval for each bin is also plotted, and the total number of plots 
# is labeled above each the bin.
#
# Unlike a typical goodness-of-fit plot from a linear regression model, with
# Presence/Absence data having all the points lay along the diagonal does not
# necessarily imply a good quality model. The ideal calibration plot for
# Presence/Absence data depends on the intended use of the model.
#
# If the model is to be used to produce probability maps, then it is indeed
# desirable that (for example) 80 percent of plots with predicted probability
# of 0.8 actually do have observed Presence. In this case, having all the
# bins along the diagonal does indicate a good model.
#
# However, if model is to be used simply to predict species presence, then
# all that is required is that some threshold exists(not necessarily 0.5)
# where every plot with a lower predicted probability is observed Absent, and
# every plot with a higher predicted probability is observed Present. In this
# case, a good model will not necessarily (in fact, will rarely) have all the
# bins along the diagonal. (Note: for this purpose 'presence.absence.hist'
# may produce more useful diagnostics.)

# Selected species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(sp in c("JUSC2","PICO")){
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			calibration.plot(	DATA, 
						which.model=mod, 
						color=mod+1,
						main=model.names[mod])}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev[species==sp,]$Obs.Prevalence,2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)}

# All species:
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,3))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:3){
			calibration.plot(	DATA, 
						which.model=mod, 
						color=mod+1,
						main=model.names[mod])}
		mtext(sp,side=3,line=1,cex=1.5,outer=T)
		mtext(paste(signif(pred.prev$Obs.Prevalence[i],2),"Prevalence"),side=3,line=-2,cex=1.2,outer=T)} 

######################## Summary Plots ######################################

# The function 'presence.absence.summary' produces a useful summary page
# of the four types of presence/absence plots, along with optimal thresholds
# and AUC value.

# For a selected species/model:

	mod<-1
	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mar=c(3,10,4,10))
	for(sp in c("JUSC2","PICO")){
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		presence.absence.summary(	DATA, 
							which.model=mod,
							main=paste(	sp,"(",
									signif(pred.prev[pred.prev$species==sp,]$Obs.Prevalence,2),
									"Prevalence ) - ",
									model.names[mod]),
							legend.cex=.9,
							opt.legend.cex=.9)}

# For all species and models:

	windows(width=9,height=6,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mar=c(3,10,4,10))
	for(i in 1:N.sp){
		sp<-species[i]
		DATA=SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:N.models){
			presence.absence.summary(	DATA, 
								which.model=mod,
								truncate.tallest=TRUE,
								main=paste(	sp,"(",
										signif(pred.prev[pred.prev$species==sp,]$Obs.Prevalence,2),
										"Prevalence ) - ",
										names(SPDATA[mod+2])),
								legend.cex=.9,
								opt.legend.cex=.9)}}



############ Effect of threshold optimization method on Predicted Prevalence ###############

# We discussed earlier how choosing a threshold so that sensitivity equals specificity
# can in theory affect the estimated prevalence of a species. This graph demonstrates
# this effect on an actual dataset.
#
# The predicted prevalence was calculated for all 13 species and three models, with
# threshold optimized by all four methods. Then the predicted prevalence was plotted 
# against observed prevalence for each method.
#
# As you can see, choosing a threshold to maximize Kappa results in unbiased estimates 
# of the true prevalence. The default method of 'threshold = 0.5' very slightly 
# underestimates the true prevalence.  
#
# 'ReqSens=0.85' overestimates prevalence, though this effect is dependent 
# on our choice of required sensitivity. If our requirements were
# less strict, this method might be unbiased, or even underestimate prevalence.
 
# The method of 'sensitivity = specificity' over estimates rare species. On the other hand,
# if there had been common species with prevalence greater than 0.50, 'sensitivity = specificity' 
# would instead underestimate their prevalence.

### Four methods - in color ###

	### calculate accuracy prevalence ###

	opt.methods<-c("Default","Sens=Spec","ReqSens","MaxKappa")

	PREVALENCE<-data.frame(matrix(0,(N.sp*N.models*length(opt.methods)),6))
	names(PREVALENCE)<-c("species","model","opt.methods","threshold","Obs.Prev","Pred.Prev")

	PREVALENCE$species<-rep(species,each=(N.models*length(opt.methods)))
	PREVALENCE$model<-rep((1:N.models),N.sp,each=length(opt.methods))
	PREVALENCE$opt.methods<-rep(opt.methods,(N.sp*N.models))
	
	for(i in 1:N.sp){
		sp<-species[i]
		DATA<-SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:N.models){
			PREVALENCE[PREVALENCE$species==sp & PREVALENCE$model==mod,]$threshold<-
				error.threshold.plot( 	DATA,
								which.model=mod,
								plot.it=FALSE,
								opt.thresholds=TRUE,
								opt.methods=opt.methods,
								req.sens=0.85)$threshold

			PREVALENCE[PREVALENCE$species==sp & PREVALENCE$model==mod,5:6]<-
				predicted.prevalence(	DATA,
								threshold=PREVALENCE[PREVALENCE$species==sp & PREVALENCE$model==mod,]$threshold,
								which.model=mod)[,2:3]
		}
	}

	### Make Graph ### 

	windows(width=9,height=6,pointsize=12,record=TRUE)

	pch.prevalence<-c(1,2,3)
	col.prevalence<-(1:N.models)+1

	op<-par(mfrow=c(2,2),pty="s",oma=c(3,2,2,0),mar=c(2,2,4,1))
	for(opt.meth in opt.methods){
		plot(	0:1,0:1,
			xlab="",ylab="",
			xlim=c(0,1),ylim=c(0,1),
			type="n")
		points(	PREVALENCE[PREVALENCE$opt.methods==opt.meth,]$Pred.Prev,
				PREVALENCE[PREVALENCE$opt.methods==opt.meth,]$Obs.Prev,
				pch=pch.prevalence[PREVALENCE[PREVALENCE$opt.methods==opt.meth,]$model],
				col=col.prevalence[PREVALENCE[PREVALENCE$opt.methods==opt.meth,]$model])
		abline(a=0,b=1)
		mtext(opt.meth,side=3,line=1,cex=1,font=2)
	}

	legend(x=.6,y=.4,legend=names(SPDATA)[2+(1:N.models)],pch=pch.prevalence,col=col.prevalence)
	mtext("Methods for Optimizing Threshold Choice",side=3,line=0,cex=1.2,outer=TRUE)
	mtext("Predicted Prevalence",side=1,line=1,cex=1,outer=TRUE)
	mtext("Observed Prevalence",side=2,line=-5,cex=1,outer=TRUE)



### All methods - in black and white ###

	opt.methods<-c("Default","Sens=Spec","MaxSens+Spec","MaxKappa","MaxPCC","PredPrev=Obs","ObsPrev","MeanProb","MinROCdist","ReqSens","ReqSpec")
	pretty.opt.methods<-c("Default","Sens=Spec","MaxSens+Spec","MaxKappa","MaxPCC","PredPrev=Obs","ObsPrev","MeanProb","MinROCdist","ReqSens=0.85","ReqSpec=0.85")

	PREVALENCE<-data.frame(matrix(0,(N.sp*N.models*length(opt.methods)),6))
	names(PREVALENCE)<-c("species","model","opt.methods","threshold","Obs.Prev","Pred.Prev")

	PREVALENCE$species<-rep(species,each=(N.models*length(opt.methods)))
	PREVALENCE$model<-rep((1:N.models),N.sp,each=length(opt.methods))
	PREVALENCE$opt.methods<-rep(opt.methods,(N.sp*N.models))
	
	for(i in 1:N.sp){
		sp<-species[i]
		DATA<-SPDATA[SPDATA$SPECIES==sp,]
		for(mod in 1:N.models){
			PREVALENCE[PREVALENCE$species==sp & PREVALENCE$model==mod,]$threshold<-
				error.threshold.plot( 	DATA,
								which.model=mod,
								plot.it=FALSE,
								opt.thresholds=TRUE,
								opt.methods=opt.methods,
								req.sens=0.85,
								req.spec=0.85)$threshold

			PREVALENCE[PREVALENCE$species==sp & PREVALENCE$model==mod,5:6]<-
				predicted.prevalence(	DATA,
								threshold=PREVALENCE[PREVALENCE$species==sp & PREVALENCE$model==mod,]$threshold,
								which.model=mod)[,2:3]
		}
	}

# Make Graph # 

	windows(width=9,height=6,pointsize=12,record=TRUE)
	op<-par(mfrow=c(3,4),cex=.8,oma=c(3,2,3,0),mar=c(3,3,3,1),pty="s",cex=.7)
	for(opt.meth in opt.methods){
		plot(	c(0,.6),c(0,.6),xlab="",ylab="",xlim=c(0,.6),ylim=c(0,.6),type="n")
		for(i in 1:N.sp){
			sp<-species[i]
			points(	PREVALENCE[PREVALENCE$opt.methods==opt.meth & PREVALENCE$species==sp,]$Pred.Prev,
					PREVALENCE[PREVALENCE$opt.methods==opt.meth & PREVALENCE$species==sp,]$Obs.Prev,
					pch=i)
			abline(a=0,b=1,col="gray")
			mtext(pretty.opt.methods[opt.methods==opt.meth],side=3,line=1,cex=1,font=2)}
	}
	par(mar=c(0,0,0,0))
	plot(0:1,0:1,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")
	legend(x="bottomright",legend=species,pch=1:N.sp,ncol=2,title="species",cex=1.1)
	mtext("Methods for Optimizing Threshold Choice",side=3,line=1,cex=1.2,outer=TRUE)
	mtext("Predicted Prevalence",side=1,line=1,cex=1,outer=TRUE)
	mtext("Observed Prevalence",side=2,line=-2,cex=1,outer=TRUE)

###########################################################################

print("The R script source code for this demo contains extensive notations on the use and interpretation of these functions. This is located in the demo folder of the package instalation.")
