# Copyright (c) 2010-2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
# 
# train_*** <- function(X, y, pars)
# e.g. train_linear or train_gp
#
# Performs a regression from X to y
# 
# INPUT:
#   X         nxp matrix of training inputs (n data points, p dimensions)
#   Y         vector of N training outputs (n data points)
#   pars      list containing parameters of the regression method
#
# OUTPUT:
#   result    list with the result of the regression
#      $model        list of learned model (e.g., weight vector)
#      $Yfit         fitted outputs for training inputs according to the learned model
#      $residuals    noise values (e.g., residuals in the additive noise case)



####
#Linear Regression
####
train_linear <- function(X,y,pars = list())
{
    mod <- lm(y ~ X)
    result <- list()
    result$Yfit = as.matrix(mod$fitted.values)
    result$residuals = as.matrix(mod$residuals)
    result$model = mod
    #for coefficients see list(mod$coef)
    return(result)
}



####
#GP Regression
####
gp_regression <- function(X,y, pars=list())
{
    options=gpOptions("ftc")
    options$kern$comp=list("rbf","white")
    #options$learnScales=TRUE
    model<-gpCreate(dim(X)[2],1,X,y,options)
    y2<-gpOut(model,X)
    model$Yfit<-y2
    model$residuals<-y-y2
    return(model)
}

train_gp <- function(X,y,pars = list())
{
    library(gptk)
    mod <- gp_regression(as.matrix(X),as.matrix(y))
    result <- list()
    result$Yfit = mod$Yfit
    result$residuals = mod$residuals
    result$model = mod
    return(result)
}



# =========
# 2. train_model
# =========

train_model <- function(f,X,y,pars = list())
{
    result <- f(X,y,pars)
}

