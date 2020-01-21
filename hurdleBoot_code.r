###################################
# Improvement Analytics Unit      #
#                                 #
# ------------------------------- #
#                                 #
# Routines for fitting Hurdle     #
# models to count data and for    #
# deriving bootstrap inferences   #
# on an intervention effect       #
#                                 #
# ------------------------------- #
#                                 #
# Author: Stefano Conti           #
# (e. stefano.conti@nhs.net)      #
###################################


##############
## Preamble ##
##############

# rm(list=ls())  # Clear work-space


library(boot)  # Load library required for bootstrap routines

library(pscl)  # Load library required for two-part modelling


bootCI_vec <- c("bca", "perc"
                )  # Set vector of estimation methods for bootstrap confidence interval


positiveHurdleModel_vec <- c("poisson", "negbin", "geometric"
                             )  # Set vector of GLM families for positive count data


##################
## User-defined ##
## parameters   ##
##################

set.seed(pi)  # Set random seed


bootSampleSize <- 5e3  # Set number of bootstrap replicates


bootCI <- bootCI_vec[1]  # Set bootstrap confidence interval estimation method


positiveHurdleModel <- positiveHurdleModel_vec[2]  # Set positive count data GLM family


## User is additionally to set the following:
## 
## data_dat: individual patient-level data-frame comprising a count outcome 
##           (namely "Outcome"), a binary intervention assignment factor 
##           (namely "Intervention") and a set of candidate outcome predictors;
##
## positiveHurdlePredictorNames_vec: character vector of predictor names, 
##                                   including any interactions among them 
##                                   and the binary intervention assignment 
##                                   factor "Intervention" (which is assumed
##                                   not to interact with any predictor) for
##                                   the positive part of the Hurdle model;
##
## 
## zeroHurdlePredictorNames_vec: character vector of predictor names, 
##                               including any interactions among them 
##                               and the binary intervention assignment 
##                               factor "Intervention" (which is assumed
##                               not to interact with any predictor) for
##                               the zero part of the Hurdle model.


####################
## Setting up the ##
## Hurdle model   ##
####################

linearPositiveHurdlePredictor <- paste0(positiveHurdlePredictorNames_vec, 
                                        collapse=" + "
                                        )  # Derive linear predictor for the Hurdle positive model part

linearZeroHurdlePredictor <- paste0(zeroHurdlePredictorNames_vec, 
                                    collapse=" + "
                                    )  # Derive linear predictor for the Hurdle zero model part


linearHurdlePredictor <- paste(linearPositiveHurdlePredictor, 
                               linearZeroHurdlePredictor, 
                               sep=" | "
                               )  # Derive linear predictor for the Hurdle model


hurdleModel <- as.formula(paste("Outcome", 
                                linearHurdlePredictor, 
                                sep=" ~ ")
                          )  # Derive formatted Hurdle model formula


###########################
## Bootstrap re-sampling ##
###########################

hurdleBoot_fn <- function(X_dat, idx_dat) {
  
  ## Function fits to continuous data a Hurdle regression model
  ## and computes inferences (point estimates and 95% confidence 
  ## intervals) on regression coefficients via bootstrap.
  ## 
  ## Arguments:
  ## 
  ## X_dat  : a data.frame with response and predictor variables
  ## idx_dat: a vector of row indices for X_dat
  ## 
  ## Output:
  ## 
  ## Predicted intervention effect, averaged acrosss the covariate 
  ## distribution, from the data sub-sample.
  
  
  stopifnot(is.data.frame(X_dat), 
            is.integer(idx_dat), 
            any(idx_dat %in% seq.int(nrow(X_dat)))
            )  # Set sanity checks against argument incongruency
  
  
  hurdleBoot_mdl <- tryCatch(hurdle(hurdleModel, data=X_dat[idx_dat, ], 
                                    weights=NULL, offset=NULL, 
                                    dist=positiveHurdleModel, zero.dist="binomial", link="logit", 
                                    control=hurdle.control(maxit=1e3, trace=FALSE)), 
                             error=function(err) 
                               cat("Supplied bootstrap sample won't fit a Logistic-", 
                                   c("Poisson", "Negative Binomial", "Geometric")[match(positiveHurdleModel_vec, 
                                                                                        table=c("poisson", "negbin", "geometric"))
                                                                                  ], 
                                   " Hurdle regression model\n")
                             )  # Fit Hurdle Binomial regression model to data-frame
  
  
  XPredict_ls <- lapply(levels(X_dat$Intervention), 
                        FUN=function(lvl) 
                          within(X_dat, 
                                 expr=assign("Intervention", value=lvl))
                        )  # Derive list of data-frames with fixed "Intervention" level
  
  
  yPredict_arr <- sapply(XPredict_ls, 
                         FUN=function(dat) 
                           predict(hurdleBoot_mdl, 
                                   newdata=dat, 
                                   type="response")
                         )  # Derive 2d-array of predicted counts by "Case", "Intervention" status
  
  dimnames(yPredict_arr) <- list(Case=seq.int(nrow(yPredict_arr)), 
                                 Intervention=levels(X_dat$Intervention)
                                 )  # Rename margins of absolute predicted costs 2d-array
  
  
  gamma_vec <- colMeans(yPredict_arr)  # Derive vector of average predicted counts by "Intervention" level
  
  
  tau <- gamma_vec[2] / gamma_vec[1]  # Derive average predicted intervention effect
  
  
  return(tau)  # Output average predicted intervention effect
  }


tauBoot_out <- boot(data_dat, statistic=hurdleBoot_fn, R=bootSampleSize, 
                    sim="ordinary", stype="i", strata=data_dat$Intervention
                    )  # Derive stratified bootstrap sample of average predicted costs by "Intervention" levels


#########################
## Post-processing     ##
## bootstrap estimates ##
#########################

tauBoot_vec <- vector("numeric", 
                      length=4)  # Initialise vector of bootstrap effect estimates by "Inference" type

names(tauBoot_vec) <- c("Point estimate", 
                        paste(c("Lower", "Upper"), 
                              "bound", 
                              sep=" "), 
                        "P-value"
                       )  # Name entries of bootstrap effect estimates vector


tauBoot_vec["Point estimate"] <- tauBoot_out$t0  # Set bootstrap effect point estimate

tauBoot_vec[paste(c("Lower", "Upper"), "bound", sep=" ")] <- boot.ci(tauBoot_out, 
                                                                     conf=.95, 
                                                                     type=bootCI, 
                                                                     index=1
                                                                    )[[4]][4:5]  # Set 95% bootstrap CI


tauBoot_lnCentred_vec <- scale(log(tauBoot_out$t), 
                               center=TRUE, scale=FALSE
                               )  # Derive vector of centred bootstrap log-effect estimates

tauBoot_vec["P-value"] <- mean(abs(tauBoot_lnCentred_vec) 
                               > abs(log(tauBoot_out$t0))
                              )  # Derive approximate bootstrap effect p-value
