
################################################################################
###############################          #######################################
###############################   DEMO   #######################################
###############################          #######################################
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

library(PresenceAbsence)



#############################################################################
####################### But first, we need a dataset...######################
#############################################################################



################### Simulating Presence/Absence Data ########################

# We will start by simulating data for two species on 1000 plot locations,
# along with model predictions from three diffent models. Just to make life
# more interesting, we will include a few NA values for one of the species.

set.seed(888)

# low prevalence species:
SPlow<-presence.absence.simulation(	n=1000,
						prevalence=.2,
						N.models=3,
						shape1.absent=c(1,1,1),
						shape2.absent=c(20,8,5), 
						shape1.present=c(6,4,2),
						shape2.present=c(1.5,4,4))

# high prevalence species:
SPhigh<-presence.absence.simulation(n=1000,
						prevalence=.8,
						N.models=3,
						shape1.absent=c(6,11,14),
						shape2.absent=c(10,8,7), 
						shape1.present=c(20,12,10),
						shape2.present=c(2,2,2))

# add NA values:

SPlow[c(3,4)  ,4]<-NA
SPlow[c(2,4,5),5]<-NA

# To make the graphs prettier, we will define a vector of pretty looking
# model names:

PRETTY.NAMES<-c("Model A", "Model B", "Model C")

# Note: the beta distribution used in this function is extremely flexible,
# and is capable of generating datasets that do not behave like typical
# presence/absence data. Some rules of thumb to get realistic distributions:
#
#	The mean of the beta distribution equals shape1/(shape1+shape2)
#
#	To get reasonable predictions (e.g. better than random), the mean for the 
#	plots where the observed value is present should be higher than that of 
#	the plots where the species is absent:
#
#		mean(present) > mean(absent)
#
#	The overall mean probability should be approximately equal to the prevalence:
#
#		prevalence*mean(present) + (1-prevalence)*mean(absent) = prevalence

################## take a look at the data format ###########################

# Let's take a look at the first few rows of the datasets:

SPlow[1:10,]
SPhigh[1:10,]

# First look at the data structure:
#
# Column 1 is the plot ID's. In this case they are sequential index numbers
# but they can be actual plot names made up of any combinations of numbers 
# or characters. They are not actually used by the functions, except for the
# fact that NA values will cause the functions to discard that row of data.
#
# Column 2 is the observed values where '0 = Absent' and '1 = Present'. If 
# the observed values are actual values (for example, basal area or tree 
# counts) any values other than zero will be treated as 'Present'.
#
# The remaining columns are for model predictions. Preferably these will
# be predicted probabilities ranging from zero to one. If all that is 
# available is predicted Presnce/Absence values, basic accuracy statistics
# can still be calculated, but you loose the ability to examine the effect
# of threshold choice, and most of the graphs are not meaningful.
#
# Note: while the demonstration datasets all have predictions from 3 models,
# the functions will work for datasets with any number of models.
#
# Note also: the functions rely on column order to identify the Observed and
# Predicted values, not on the column names. Therefor, as long as the columns 
# are in the correct order, column names can be anything you choose. The only
# time the functions use the names is for the default labels on graphs.

#############################################################################
########################## ANALYTICAL ANALYSIS ##############################
#############################################################################



###################### Presence Absence Histograms ##########################

# Let's get a visual feel for the data. The 'presence.absence.hist' function
# is a bar plot of observed values as a function predicted probability:

win.graph(width=7,height=6,pointsize=12)
op<-par(mfrow=c(2,3),oma=c(0,3,0,0))
for(i in 1:3){
presence.absence.hist(	DATA=SPlow,
				which.model=i,
				model.names=PRETTY.NAMES,
				N.bars=20,
				na.rm=TRUE)
}
for(i in 1:3){
presence.absence.hist(	DATA=SPhigh,
				which.model=i,
				model.names=PRETTY.NAMES,
				N.bars=20)
}
mtext("SPlow" ,side=2,line=1,cex=1.2,at=.75,outer=T)
mtext("SPhigh",side=2,line=1,cex=1.2,at=.25,outer=T)
par(op)

# These plots are valuable in two ways:
#
# First, these plots make the different prevalences of the two species 
# dramatically obvious. Many of the traditional accuracy measures (such as 
# percent correctly classified (PCC), sensitivity and specificity) are highly
# dependant on species prevalence. Therefor, it is important to check
# prevalence early in an investigation.
#
# Second, these plots are a good visual indicator of model quality. 
#
# A good model (such as Model A) will have a distictly double humped histogram. If
# you choose a threshold in the valley (not nessarily at 0.5!) it will divide
# the data cleanly into Observe Present and Observed Absent.
#
# A poor model (such as Model C) will have considerable overlap between
# Observed Present and Observed Absent. No threshold exists that will cleanly
# divide the data into Observe Present and Observed Absent. The challenge here
# is to find a threshold that offers the best compromise. 
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
# 'predicted.prevalence' lets us make a table of the effect of threshold
# choice on the predicted prevalence of our models:

PREV.low <-predicted.prevalence(SPlow, threshold=11, na.rm=TRUE)
PREV.high<-predicted.prevalence(SPhigh, threshold=11)

PREV.low
PREV.high

# by setting 'threshold = 11' the function calculates 11 evenly spaced
# thresholds between zero and one. If you are interested in particular 
# thresholds you can instead set 'threshold = c(.23,.6,.97)'.



########### A Digression on how the functions handle NA values ##############

# You may have noticed that for 'SPlow' we have been setting 'na.rm = TRUE' 
# and that there has been a print out in the command window stating:
# "4 rows ignored due to NA values".
#
# Let's experament a little to get a feel for what is going on:

predicted.prevalence(SPlow, threshold=11)

# In this case, the function simply returns 'NA'. The default is for 
# 'na.rm = FALSE'. To get the function to work, we must set 'na.rm = TRUE':

predicted.prevalence(SPlow, threshold=11,na.rm=TRUE)

# Next let's check how many NA's are present in each column of SPlow:

apply(is.na(SPlow),2,sum)

# Now what is going on? There are 5 NA's total, with no more than 3 in any
# particular column, but the print out stated that 4 rows were removed?
#
# The answer is that if an NA is present in any entry (including plotID!)
# the entire row is discarded. The reason only 4 rows were discarded is that
# one row has two NA entries.
#
# Now let's play with the 'which.model' argument and see how it affects the
# treatment of NA values. The Model A doesn't have any NA's, so let's see
# how it behaves:

predicted.prevalence(SPlow, threshold=11, which.model=1)

# Hmmm...it is still returning NA. Let's try setting 'na.rm=TRUE':

predicted.prevalence(SPlow, threshold=11, which.model=1,na.rm=TRUE)

# It is still discarding the 4 rows even though none of the NA's are actually
# in the predictions for Model A. If we want to calculate statistics for every
# plot from a particular model, including those plots who have NA predictions 
# from other models, we have to subset the data ouselves:

predicted.prevalence(SPlow[,c(1,2,3)],threshold=11,na.rm=TRUE)
predicted.prevalence(SPlow[,c(1,2,4)],threshold=11,na.rm=TRUE)
predicted.prevalence(SPlow[,c(1,2,5)],threshold=11,na.rm=TRUE)



###################### Accuracy Tables #######################################

# We will make some tables of the accuracy statistics.
#
# First we will look at all three models at a single threshold:

presence.absence.accuracy(SPlow,na.rm=TRUE)
presence.absence.accuracy(SPhigh)

# Note: The standard deviations in this table are only apropriate for
# examining each model individually. To do significance testing for differences
# between models it is nessassary to take into account correlation. See
# 'help(auc)' for further details.
#
# We can suppress the standard deviations by setting 'st.dev=FALSE'
#
# The default for 'presence.absence.accuracy' is to calculate accuracy measures for 
# all models at 'threshold = 0.5'. We can also substitute a different 
# threshold:

presence.absence.accuracy(SPhigh,threshold=0.4,st.dev=FALSE)

# When the table is being produced for multiple models, 'threshold' must be
# either of length 1, or it must be the same length as the number of models
# in the table. In the later case, the first threshold will be used for the 
# first model, the second threshold for the second row, etc...

presence.absence.accuracy(SPhigh, which.model=c(1,3),threshold=c(0.4,.5),st.dev=FALSE)

# 'presence.absence.accuracy' can also be used to produce a table for a single 
# model, illustrating how the accuracy measures change as a function of threshold:

presence.absence.accuracy(SPhigh, which.model=3,threshold=11,st.dev=FALSE)
presence.absence.accuracy(SPhigh, which.model=3,threshold=c(.2,.4,.5,.89),st.dev=FALSE)

# You will notice on these tables, that the area under the curve (AUC) is
# constant for all thresholds. The AUC is the only one of these basic accuracy
# statistics that is threshold independant.



################# Graphing the basic Accuracy statistics #####################

# If you need to look up specific values, tables are nessassary, but often it
# is easier to spot patterns in a graph. Lets take a graphical look at how
# these error statistics vary with threshold:

win.graph(width=7,height=6,pointsize=12)
op<-par(mfrow=c(2,3),oma=c(0,3,3,0))
for(i in 1:3){
error.threshold.plot(	DATA=SPlow,
				which.model=i,
				main=PRETTY.NAMES[i],
				na.rm=TRUE,
				color=TRUE,
				opt.method="MaxKappa")
}
for(i in 1:3){
error.threshold.plot(	DATA=SPhigh,
				which.model=i,
				main=PRETTY.NAMES[i],
				color=TRUE,
				opt.method="MaxKappa")
}
mtext("SPlow" ,side=2,line=1,cex=1.2,at=.75,outer=T)
mtext("SPhigh",side=2,line=1,cex=1.2,at=.25,outer=T)
mtext("Error Statistics vs. Threshold",side=3,line=0,cex=1.5,outer=T)
par(op)

# These plots always include sensitivity and specificity. 
# Some optimization methods from 'opt.method' have associated accuracy
# statistics, which will be plotted, even if 'opt.thresholds' = 'False'.

################### Optimizing Thresholds ##################################

# These graphs can be used in several ways to optimize threshold choice, and
# in fact, the function 'error.threshold.plot' can calculate these optimized
# thresholds as well as add them to the plots:

i<-3
win.graph(width=9,height=5,pointsize=12)
op<-par(mfrow=c(1,2),oma=c(0,0,3,0),mar=c(5,4,2,1),cex=.8)

error.threshold.plot(	DATA=SPlow,
				which.model=i,
				main="SPlow",
				na.rm=TRUE,
				color=TRUE,
				opt.thresholds=TRUE,
				opt.methods=c(1,2,10,4))

error.threshold.plot(	DATA=SPhigh,
				which.model=i,
				main="SPhigh",
				color=TRUE,
				opt.thresholds=TRUE,
				opt.methods=c(1,2,10,4))

mtext(PRETTY.NAMES[i],side=3,line=0,cex=1.5,outer=T)
par(op)

# The PresenceAbsence library provides multiple methods for optimizing 
# thresholds:
#
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
#
# Note also: 'error.threshold.plot' is a rough and ready applied function. 
# It optimizes thresholds simply by calculating a large number of evenly spaced 
# thresholds and looking for the best ones. This is good enough for graphs,
# but to find the theoretically 'best' thresholds, would require calculating
# every possible unique threshold (not necessarily evenly spaced!).

################################# ROC plots #############################################

# ROC plots, with their associated AUC (Area Under the Curve) provide a
# threshold independent method of assessing model performance. 

	windows(width=9,height=5,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,2))

	auc.roc.plot(	DATA=SPlow, 
				color=TRUE,
				legend.cex=.7,
				main="SPlow",
				model.names=PRETTY.NAMES,
				na.rm=TRUE)

	auc.roc.plot(	DATA=SPhigh, 
				color=TRUE,
				legend.cex=.7,
				main="SPhigh",
				model.names=PRETTY.NAMES)

	mtext("ROC Plots",side=3,line=1,cex=1.5,outer=TRUE)
				
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
# The argument 'mark' can be used to label particular thresholds along the ROC curve:

	windows(width=9,height=5,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,2))

	auc.roc.plot(	DATA=SPlow, 
				color=TRUE,
				legend.cex=.6,
				model.names=PRETTY.NAMES,
				main="SPlow",
				mark=.5,
				na.rm=TRUE)

	auc.roc.plot(	DATA=SPhigh, 
				color=TRUE,
				legend.cex=.6,
				model.names=PRETTY.NAMES,
				mark=.5,
				main="SPhigh")

	mtext("ROC Plots",side=3,line=1,cex=1.5,outer=TRUE)

# the argument 'optimal.thresholds' will mark the optimized thresholds from
# 'error.threshold.plot' along the ROC curves:

	windows(width=9,height=5,pointsize=12,record=TRUE)
	par(oma=c(0,0,5,0),mfrow=c(1,2))

	auc.roc.plot(	DATA=SPlow, 
				color=TRUE,
				legend.cex=.6,
				main="SPlow",
				model.names=PRETTY.NAMES,
				opt.thresholds=TRUE,
				opt.methods=c(1,2,4),
				na.rm=TRUE)

	auc.roc.plot(	DATA=SPhigh, 
				color=TRUE,
				legend.cex=.6,
				opt.thresholds=TRUE,
				opt.methods=c(1,2,4),
				main="SPhigh",
				model.names=PRETTY.NAMES)

	mtext("ROC Plots",side=3,line=1,cex=1.5,outer=TRUE)

###########################################################################

cat("The R script source code for this demo contains extensive
notations on the use and interpretation of these functions.
This is located in the demo folder of the package instalation.",
fill=TRUE,sep="")

