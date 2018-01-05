## This file is structured as follows:
## 1) required libraries are loaded,
## 2) required learners are wrapped in the R package mlr,
## 3) measures are defined,
## 4) benchmark experiments are defined for each dataset (estimated mean perfor
##    mances are computed after benchmarking),
## 5) and boxplot is generated as a .tex file, which needs to be compiled to ob-
##    tain the final figure presented in the paper.

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(survivalsvm)
library(mlr)
library(caret)
library(irace)# for iterated F-racing
library(Hmisc)# for rcorr.cens
library(Rmisc)
library(BBmisc)
library(checkmate)
library(survival)
library(survpack)
library(ggplot2)
library(mlr)
library(plyr)
require(tikzDevice)
# ------------------------------------------------------------------------------
# Wrap learners into the mlr package
# ------------------------------------------------------------------------------

#--- creation of an mlr leaner for survivalsvm
makeRLearner.surv.survivalsvm = function() {
  makeRLearnerSurv(
    cl = "surv.survivalsvm",
    package = "survivalsvm",
    par.set = makeParamSet(
      makeDiscreteLearnerParam(id = "type",  default = "regression",
                               values = c("regression", "vanbelle1", "vanbelle2",
                                          "hybrid")),
      makeDiscreteLearnerParam(id = "diff.meth",  default = "makediff3",
                               values = c("makediff1", "makediff2",
                                          "makediff3")),
      makeNumericVectorLearnerParam(id = "gamma.mu", 
                                    tunable = TRUE, lower = 2^-5, upper = 2^5),
      makeDiscreteLearnerParam(id = "opt.meth", default = "quadprog",
                               values = c("quadprog", "ipop")),
      makeDiscreteLearnerParam(id = "kernel", default = "lin_kernel",
                               values = c("lin_kernel", "add_kernel",
                                          "rbf_kernel",
                                          "rbf4_kernel", "poly_kernel")),
      makeNumericLearnerParam(id = "kernel.pars", tunable = TRUE),
      makeNumericLearnerParam(id = "sgf.sv", default = 5, tunable = FALSE),
      makeNumericLearnerParam(id = "sigf", default = 7, tunable = FALSE),
      makeNumericLearnerParam(id = "maxiter", default = 20, tunable = FALSE),
      makeNumericLearnerParam(id = "margin", default = 0.05, tunable = FALSE),
      makeNumericLearnerParam(id = "bound", default = 10, tunable = FALSE),
      makeNumericLearnerParam(id = "eig.tol", default = 1e-06,
                              tunable = FALSE),
      makeNumericLearnerParam(id = "conv.tol", default = 1e-07,
                              tunable = FALSE),
      makeNumericLearnerParam(id = "posd.tol", default = 1e-08,
                              tunable = FALSE)
    ),
    properties = c("missings", "numerics", "factors", "weights", "prob",
                   "rcens"),
    name = "survival support vector machines",
    short.name = "survivalsvm",
    note = "survivalsvm in mlr"
  )
}

#--- creation of trainer for survivalsvm
trainLearner.surv.survivalsvm = function(.learner, .task, .subset, ...) {
  f <-  getTaskFormula(.task)
  data  <-  getTaskData(.task, subset = .subset)
  mod <-  survivalsvm::survivalsvm(formula = f, data = data, ...)
  return(mod)
}

#--- creation of predictor for survivalsvm
predictLearner.surv.survivalsvm = function(.learner, .model, .newdata, ...) {
  if (.learner$predict.type == "response") {
    predict(object = .model$learner.model,
            newdata = .newdata, ...)$predicted[1,]
  }
}

#--- Wrapper for scale data when required
makePreprocWrapperScale = function(learner, center = TRUE, scale = TRUE) {
  trainfun = function(data, target, args = list(center, scale)) {
    cns = colnames(data)
    nums = setdiff(cns[sapply(data, is.numeric)], target)
    x = as.matrix(data[, nums, drop = FALSE])
    x = scale(x, center = args$center, scale = args$scale)
    control = args
    if (is.logical(control$center) && control$center)
      control$center = attr(x, "scaled:center")
    if (is.logical(control$scale) && control$scale)
      control$scale = attr(x, "scaled:scale")
    data = data[, setdiff(cns, nums), drop = FALSE]
    data = cbind(data, as.data.frame(x))
    return(list(data = data, control = control))
  }
  predictfun = function(data, target, args, control) {
    cns = colnames(data)
    nums = cns[sapply(data, is.numeric)]
    x = as.matrix(data[, nums, drop = FALSE])
    x = scale(x, center = control$center, scale = control$scale)
    data = data[, setdiff(cns, nums), drop = FALSE]
    data = cbind(data, as.data.frame(x))
    return(data)
  }
  if(!("surv.glmboost" %in% class(learner)))
    makePreprocWrapper(
      learner,
      train = trainfun,
      predict = predictfun,
      par.set = makeParamSet(
        makeLogicalLearnerParam("center"),
        makeLogicalLearnerParam("scale")
      ),
      par.vals = list(center = center, scale = scale)
    )
  else
    makePreprocWrapper(
      learner,
      train = trainfun,
      predict = predictfun,
      par.set = makeParamSet(),
      par.vals = list(center = center, scale = scale)
    )
}
# --- TuneWrapper for survivalsvm 
tuneWrapperGammaMu <- function(type, kernel, ...,
                               preproc = TRUE, tune.scale = TRUE, center = TRUE,
                               scale = TRUE, resolution = 10L, method = "CV",
                               lower = -10L, upper = 10L, iters.rep = 9L) {
  lrn <- makePreprocWrapperScale(makeLearner(cl = "surv.survivalsvm",
                                             type = type, kernel = kernel,
                                             ...),
                                 center = center, scale = scale)
  configureMlr(on.learner.error = "warn")
  if(kernel %in% c("lin_kernel", "add_kernel")) {
    # Parameters set for linear and additiv kernel
    if (type != "hybrid") {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeDiscreteParam("gamma.mu", values = 2^(lower:upper)),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale")
        )
      } else {
        makeParamSet(
          makeDiscreteParam("gamma.mu",
                            values = 2^(lower:upper))
        )
      }
    } else {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale")
        )
      } else {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x)
        )
      }
    }
  } else {
    # Parameters set for RBF kernel
    if (type != "hybrid") {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeDiscreteParam("gamma.mu", values = 2^(lower:upper)),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5)),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale")
        )
      } else {
        makeParamSet(
          makeDiscreteParam("gamma.mu",
                            values = 2^(lower:upper)),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5))
        )
      }
    } else {
      discrete_ps <- if(tune.scale) {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x),
          makeLogicalLearnerParam("center"),
          makeLogicalLearnerParam("scale"),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5))
        )
      } else {
        makeParamSet(
          makeNumericVectorParam("gamma.mu", len = 2L, lower = lower,
                                 upper = upper,
                                 trafo = function(x) 2^x),
          makeDiscreteParam("kernel.pars", values = 2^(-5:5))
        )
      }
    }
    
  }
  
  ctrl  <-  makeTuneControlGrid(resolution = resolution)
  inner  <-  makeResampleDesc(method = method, iters = iters.rep)
  survivalsvm.tuned  <-  makeTuneWrapper(lrn,  resampling = inner,
                                         par.set = discrete_ps, control = ctrl,
                                         measures = c.i)
  return(survivalsvm.tuned)
}

#--- creation of an mlr learner for survpack
makeRLearner.surv.survsvm = function() {
  makeRLearnerSurv(
    cl = "surv.survsvm",
    package = "survpack",
    par.set = makeParamSet(
      makeDiscreteLearnerParam(id = "kernel", default = "linear", 
                               values = c("linear", "gaussian")),
      makeNumericLearnerParam(id = "gamma", tunable = TRUE, lower = 2^-5, 
                              upper = 2^+5),
      makeNumericLearnerParam(id = "cost", tunable = TRUE, lower = 1e-3, 
                              upper = 1e+3),
      makeLogicalLearnerParam(id = "rescale", default = TRUE),
      makeLogicalLearnerParam(id = "cox.score", default = FALSE),
      makeLogicalLearnerParam(id = "thresh", default = TRUE),
      makeLogicalLearnerParam(id = "binary", default = FALSE), 
      makeNumericLearnerParam(id = "quant", tunable = FALSE, default = 0.9),
      makeNumericLearnerParam(id = "fudge", tunable = FALSE, default = 0.5),
      makeDiscreteLearnerParam(id = "max.iter", tunable = FALSE,
                               default = 20L, values = 0:1e+5),
      makeNumericLearnerParam(id = "alpha.tol", tunable = FALSE,
                              default = 1e-8),
      makeNumericLearnerParam(id = "gap.tol", tunable = FALSE, default = 1e-12)
    ),
    properties = c("missings", "numerics", "factors", "weights", "prob",
                   "rcens"),
    name = "survpack survival support vector machines",
    short.name = "survsvm",
    note = "survpack in mlr"
  )
}

# --- TuneWrapper for surpack 
tuneWrapperGammaCost <- function(kernel, ...,
                                 preproc = TRUE, tune.scale = TRUE,
                                 center = TRUE, scale = TRUE, resolution = 10L,
                                 method = "CV", lower = -10L, upper = 10L,
                                 iters.rep = 9L) {
  lrn <- makePreprocWrapperScale(makeLearner(cl = "surv.survsvm",
                                             kernel = kernel, ...),
                                 center = center, scale = scale)
  configureMlr(on.learner.error = "quiet")
  if(kernel == "linear") {
    survsvm_grid <- makeParamSet(makeDiscreteParam("gamma",
                                                   values = 2^(lower:upper))
    )
  } else {
    survsvm_grid <- makeParamSet(makeDiscreteParam("gamma",
                                                   values = 2^(lower:upper)),
                                 makeDiscreteParam("cost",
                                                   values = 2^(-(lower/2):
                                                                 (upper/2))))
  }
  inner <- makeResampleDesc(method = method, iters = iters.rep)
  ctrl <- makeTuneControlGrid(resolution = resolution)
  survsvm.tuned <- makeTuneWrapper(lrn, resampling = inner, 
                                   par.set = survsvm_grid,
                                   control = ctrl,
                                   measures = c.i.survsvm)
}
#--- creation of trainer for survsvm
trainLearner.surv.survsvm = function(.learner, .task, .subset, ...) {
  f <-  getTaskFormula(.task)
  data  <-  getTaskData(.task, subset = .subset)
  mod <-  survpack::survsvm(formula = f, data = data, ...)
  return(mod)
}

#--- creation of predictor for survsvm
predictLearner.surv.survsvm = function(.learner, .model, .newdata, ...) {
  if (.learner$predict.type == "response") {
    predict(object = .model$learner.model, newdata = .newdata, ...)
  }
}

# --- TuneWrapper for random forest
tuneWrapperRandomForestSRC <- function(...,
                                       center = TRUE, scale = TRUE,
                                       method = "CV", iters.rep = 9L) {
  lrn <- makePreprocWrapperScale(makeLearner(cl = "surv.randomForestSRC",
                                             ...),
                                 center = center, scale = scale)
  configureMlr(on.learner.error = "warn")
  rf_grid <- makeParamSet(makeIntegerParam("mtry", 1, 5),
                          makeDiscreteParam("nodesize", c(1, 3, 5, 8, 10)))
  
  inner  <-  makeResampleDesc(method = method, iters = iters.rep)
  ## Use 20 random parameter combinations for each learner
  ctrl <- makeTuneControlRandom(maxit = 20)
  rf_tuned <- makeTuneWrapper(lrn, resampling = inner,
                              par.set = rf_grid, control = ctrl)
  return(rf_tuned)
}

#--- TuneWrapper for cox model
wrapperPHModell <- function(center = TRUE, scale = TRUE,
                            ...) {
  #task <- makeSurvTask(data = data, target = target)
  lrn <- makePreprocWrapperScale(makeLearner(cl = "surv.coxph", ...),
                                 center = center, scale = scale)
  return(lrn)
}

#--- TuneWrapper for gradient boosting
tuneWrapperGBoost <- function(...,
                              center = TRUE, scale = TRUE,
                              method = "CV", iters.rep = 9L) {
  lrn <- makePreprocWrapperScale(makeLearner(cl = "surv.glmboost",
                                             ...),
                                 center = center, scale = scale)
  configureMlr(on.learner.error = "warn")
  gb_grid <- makeParamSet(makeIntegerParam("mstop", 1e2, 1e3),
                          makeDiscreteParam("nu", c(0.05, 0.1, 0.3, 0.5, 0.8,
                                                    1)))
  
  inner  <-  makeResampleDesc(method = method, iters = iters.rep)
  ## Use 20 random parameter combinations for each learner
  ctrl <- makeTuneControlRandom(maxit = 20)
  gb_tuned <- makeTuneWrapper(lrn, resampling = inner,
                              par.set = gb_grid, control = ctrl)
  return(gb_tuned)
}


# ------------------------------------------------------------------------------
# Wrap measures into the mlr package
# ------------------------------------------------------------------------------
#' my.ci.fun computes the C-index for "survivalsvmprediction" objects.
#'
#' @param task survival task used to fit the model.
#'        This will be later used by "makeMeasure" from the pakage "mlr".
#' @param model model of class "survivalsvm".
#'        This will be later used by "makeMeasure" from the pakage "mlr".
#' @param pred The predictions.
#' @param feats will be later used by "makeMeasure" from the pakage "mlr".
#' @param extra.args will be later used by "makeMeasure" from the pakage "mlr".
#'
#' @return The C-Index
#'
my.ci.fun <- function(task, model, pred, feats, extra.args) {
  myci = rcorr.cens(x = getPredictionResponse(pred),
                    S = getPredictionTruth(pred))
  return(myci["C Index"])
}
# constructs the C-index measure for survivalsvmprediction objects.
c.i <- makeMeasure(
  id = "ci", name = "C-Index",
  properties = c("surv"),
  minimize = FALSE, best = 1, worst = 0,
  fun = my.ci.fun
)

# # cindex for surxvm objects from survpack package
c.i.survsvm <- makeMeasure(
  id = "ci", minimize = FALSE, best = 1, worst = 0,
  properties = c("surv", "req.pred", "req.truth"),
  name = "C-Index",
  fun = function(task, model, pred, feats, extra.args) {
    truth <- getPredictionTruth(pred)
    response <- getPredictionResponse(pred)
    myci = rcorr.cens(x = response, S = truth)
    return(1 - myci["C Index"])
  }
)

# constructs the C-index measure for PH Model.
ph.cindex = makeMeasure(id = "ci", minimize = FALSE, best = 1, worst = 0,
                        properties = c("surv", "req.pred", "req.truth"),
                        name = "C-Index",
                        fun = function(task, model, pred, feats, extra.args) {
                          requirePackages("Hmisc", default.method = "load")
                          resp = pred$data$response
                          if (anyMissing(resp))
                            return(NA_real_)
                          s = Surv(pred$data$truth.time, pred$data$truth.event)
                          Hmisc::rcorr.cens(-1 * resp, s)[["C Index"]]
                        }
)

# --- Log-rank statistic ######
#
#' @param t1 Observations in group 1
#' @param d1 Predictions in group 1
#' @param t2 Observations in group 2
#' @param d2 Predictions in group 2
#'
#' @return log-rank statistik
logrank <- function(t1, d1, t2, d2){
  t1.ord <- order(t1)
  t1 <- t1[t1.ord]
  d1 <-  d1[t1.ord]
  t2.ord <- order(t2)
  t2 <- t2[t2.ord]
  d2 <-  d2[t2.ord]
  t <- c(t1, t2)
  d <- c(d1, d2)
  times <- unique(t[d == 1])
  n <- length(times)
  tel.ner <- sapply(times, function(i){
    o <- sum(t == i) # failures
    r <- sum(t >= i) # at risk
    o1 <- sum(t1 == i) # failures
    r1 <- sum(t1 >= i) # at risk
    o2 <- sum(t2 == i) # failures
    r2 <- sum(t2 >= i) # at risk
    teller <- o1 - r1*o/r
    noemer <- r2*r1*o*(r-o) / (r^2 *(r - 1))
    return(c(teller, noemer))
  })
  search.na <- colSums(tel.ner)
  na.index <- which(is.na(search.na))
  if(length(na.index) > 0){
    tel.ner <- tel.ner[, -na.index] 
  }
  chi <- sum(tel.ner[1,])^2 / sum(tel.ner[2,])
  return(chi_sq = chi)
}

#' @param task a survival task
#' @param model survival model
#' @param pred survival predictions
#' @param feats features
#' @param extra.args further arguments
#'
#' @return log-rank statistik
lr.fun <- function(task, model, pred, feats, extra.args){
  t <- getPredictionTruth(pred)[,1]
  delta <- getPredictionTruth(pred)[,2]
  u <- getPredictionResponse(pred)
  part2 <- which(u > mean(u))
  part1 <- setdiff(1:length(t), part2)
  return(logrank(t1 = t[part1], d1 = delta[part1], t2 = t[part2],
                 d2 = delta[part2]))
}

lgrk <- makeMeasure(
  id = "logrank", name = "Log-Rank",
  properties = c("surv"),
  minimize = FALSE, best = Inf, worst = 0,
  fun = lr.fun
)

# --- Hazardratio
#' hr.fun computes the hazard rate for survivalsvm objects given:
#'
#' @param task a survival task
#' @param model a survivalsvm model
#' @param pred object of class survivalsvm
#' @param feats the features
#' @param extra.args further arguments
#'
#' @return hazard rate
hr.fun <- function(task, model, pred, feats, extra.args){
  u <- getPredictionResponse(pred)
  a <- min(u)
  b <- max(u)
  u <- (u-a) / (b-a)
  mod.ph <- try(coxph(getPredictionTruth(pred) ~ u), silent = TRUE)
  if(!("list" %in% is(mod.ph))) return(1)
  return(exp(mod.ph$coef))
}

# construction of the measure
hr <- makeMeasure(
  id = "hr", name = "Hazardrate",
  properties = c("surv"),
  minimize = FALSE, best = Inf, worst = 0,
  fun = hr.fun
)

# --- Hazardrate
#' computes the hazard rate for reference models given:
#'
#' @param task a task
#' @param model a model
#' @param pred predictions
#' @param feats features
#' @param extra.args further arguments
#'
#' @return hazard rate
hrph.fun <- function(task, model, pred, feats, extra.args){
  u <- getPredictionResponse(pred)
  a <- min(u)
  b <- max(u)
  u <- (u-a) / (b-a)
  mod.ph <- try(coxph(getPredictionTruth(pred) ~ u), silent = TRUE)
  return(exp(-mod.ph$coef))
}

hrph <- makeMeasure(
  id = "hr", name = "Hazardrate",
  properties = c("surv"),
  minimize = FALSE, best = Inf, worst = 0,
  fun = hrph.fun
)

# ------------------------------------------------------------------------------
# Benchmarking for the veteran dataset
# ------------------------------------------------------------------------------

# #---------- Some pre-processing steps ----------------------------------------
data(veteran, package = "survival")
veteran.adj <- veteran
veteran.adj[, "prior"] <- as.factor(veteran$prior == 10)
veteran.adj[, "trt"] <- as.factor(veteran$trt == 2)
veteran.task <- makeSurvTask(data = veteran.adj, target = c("time", "status"))
outer <- makeResampleDesc("CV", iters = 5L)

# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.reg <- tuneWrapperGammaMu(type = "regression",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog", method = "CV",
                                    iters.rep = 5)
bench.vet.ssvm.reg <- benchmark(learners = tunwrp.gm.reg, tasks = veteran.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
lapply(bench.vet.ssvm.reg$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog",
                                    diff.meth = "makediff3")
bench.vet.ssvm.vb1 <- benchmark(learners = tunwrp.gm.vb1, tasks = veteran.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
save(bench.vet.ssvm.vb1, file = "benchVetSsvmVb1.rda")
lapply(bench.vet.ssvm.vb1$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog",
                                    diff.meth = "makediff3")
bench.vet.ssvm.vb2 <- benchmark(learners = tunwrp.gm.vb2, tasks = veteran.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
lapply(bench.vet.ssvm.vb2$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog",
                                    diff.meth = "makediff3")
bench.vet.ssvm.hyb <- benchmark(learners = tunwrp.gm.hyb, tasks = veteran.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
lapply(bench.vet.ssvm.hyb$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### additiv Kernel #########################################*
# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.ak.reg <- tuneWrapperGammaMu(type = "regression",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog", method = "CV",
                                       iters.rep = 5)
bench.vet.ssvm.ak.reg <- benchmark(learners = tunwrp.gm.ak.reg,
                                   tasks = veteran.task,
                                   resamplings = outer, measures = list(c.i,
                                                                        lgrk,
                                                                        hr))
lapply(bench.vet.ssvm.ak.reg$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.ak.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog",
                                       diff.meth = "makediff3")
bench.vet.ssvm.ak.vb1 <- benchmark(learners = tunwrp.gm.ak.vb1,
                                   tasks = veteran.task,
                                   resamplings = outer, measures = list(c.i,
                                                                        lgrk,
                                                                        hr))
lapply(bench.vet.ssvm.ak.vb1$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.ak.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog",
                                       diff.meth = "makediff3")
bench.vet.ssvm.ak.vb2 <- benchmark(learners = tunwrp.gm.ak.vb2,
                                   tasks = veteran.task,
                                   resamplings = outer, measures = list(c.i,
                                                                        lgrk,
                                                                        hr))
lapply(bench.vet.ssvm.ak.vb2$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.Ak.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog",
                                       diff.meth = "makediff3")
bench.vet.ssvm.ak.hyb <- benchmark(learners = tunwrp.gm.Ak.hyb,
                                   tasks = veteran.task,
                                   resamplings = outer, measures = list(c.i,
                                                                        lgrk,
                                                                        hr))
lapply(bench.vet.ssvm.ak.hyb$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# ################### RBF/Gaussian kernel #####################################*
# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.rbf.reg <- tuneWrapperGammaMu(type = "regression",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog", method = "CV",
                                        iters.rep = 5)
bench.vet.ssvm.rbf.reg <- benchmark(learners = tunwrp.gm.rbf.reg,
                                    tasks = veteran.task,
                                    resamplings = outer, measures = list(c.i,
                                                                         lgrk,
                                                                         hr))
lapply(bench.vet.ssvm.rbf.reg$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.rbf.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3")
bench.vet.ssvm.rbf.vb1 <- benchmark(learners = tunwrp.gm.rbf.vb1,
                                    tasks = veteran.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.vet.ssvm.rbf.vb1$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.rbf.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3")
bench.vet.ssvm.rbf.vb2 <- benchmark(learners = tunwrp.gm.rbf.vb2,
                                    tasks = veteran.task,
                                    resamplings = outer, measures = list(c.i,
                                                                         lgrk,
                                                                         hr))
lapply(bench.vet.ssvm.rbf.vb2$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.rbf.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3")
bench.vet.ssvm.rbf.hyb <- benchmark(learners = tunwrp.gm.rbf.hyb,
                                    tasks = veteran.task,
                                    resamplings = outer, measures = list(c.i,
                                                                         lgrk,
                                                                         hr))
lapply(bench.vet.ssvm.rbf.hyb$results$veteran.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# ## REFERENCE METHODS ########################################################*

# --- randomForestSRC.tuned
set.seed(123)
rfsrc.tuned <- tuneWrapperRandomForestSRC(center = TRUE, scale = TRUE,
                                          method = "CV", iters.rep = 10L)
bench.vet.rfsrc <- benchmark(rfsrc.tuned, veteran.task, outer,
                             measures = list(ph.cindex, lgrk, hrph))
lapply(bench.vet.rfsrc$results$veteran.adj$surv.randomForestSRC.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- coxph
set.seed(123)
lrn.ph <- wrapperPHModell(center = TRUE, scale = TRUE)
bench.vet.ph <- benchmark(lrn.ph, veteran.task, resamplings = outer,
                          measures = list(ph.cindex, lgrk, hrph))
lapply(bench.vet.ph$results$veteran.adj$surv.coxph.preproc$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- Gradient boosting 
set.seed(123)
lrn.gboost <- tuneWrapperGBoost(center = TRUE, scale = TRUE)
bench.vet.gboost <- benchmark(lrn.gboost, veteran.task, resamplings = outer,
                              measures = list(ph.cindex, lgrk, hrph))
lapply(bench.vet.gboost$results$veteran.adj$surv.glmboost.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# --- survsvm linear kernel
set.seed(123)
tunwrp.survsvm.linear <- tuneWrapperGammaCost(kernel = "linear",
                                              method = "CV", max.iter = 50,
                                              iters.rep = 5)
bench.vet.survsvm.lin <- benchmark(learners = tunwrp.survsvm.linear, 
                                   tasks = veteran.task,
                                   resamplings = outer,
                                   measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.vet.survsvm.lin$results$veteran.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survsvm RBF kernel
set.seed(123)
tunwrp.survsvm.gaussian <- tuneWrapperGammaCost(kernel = "gaussian",
                                                method = "CV", max.iter = 50,
                                                iters.rep = 5)
bench.vet.survsvm.gaus <- benchmark(learners = tunwrp.survsvm.gaussian, 
                                    tasks = veteran.task,
                                    resamplings = outer,
                                    measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.vet.survsvm.gaus$results$veteran.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# --- Benchmark mean runtime for survivalsvm approach models
bmr_vet_survivalsvm_reg <- list(bench.vet.ssvm.reg, bench.vet.ssvm.ak.reg,
                                bench.vet.ssvm.rbf.reg)
bmr_vet_survivalsvm_reg_times <- vector(mode = "list", length = 3L)
names(bmr_vet_survivalsvm_reg_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_vet_survivalsvm_reg[[bmr]])
  bmr_vet_survivalsvm_reg_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_vet_survivalsvm_reg_times, function(l){
  round(mean(l)/60, 2)
})

bmr_vet_survivalsvm_vb1 <- list(bench.vet.ssvm.vb1, bench.vet.ssvm.ak.vb1,
                                bench.vet.ssvm.rbf.vb1)
bmr_vet_survivalsvm_vb1_times <- vector(mode = "list", length = 3L)
names(bmr_vet_survivalsvm_vb1_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_vet_survivalsvm_vb1[[bmr]])
  bmr_vet_survivalsvm_vb1_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_vet_survivalsvm_vb1_times, function(l){
  round(mean(l)/60, 2)
})

bmr_vet_survivalsvm_vb2 <- list(bench.vet.ssvm.vb2, bench.vet.ssvm.ak.vb2, 
                                bench.vet.ssvm.rbf.vb2)
bmr_vet_survivalsvm_vb2_times <- vector(mode = "list", length = 3L)
names(bmr_vet_survivalsvm_vb2_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_survivalsvm_vb2[[bmr]])
  bmr_vet_survivalsvm_vb2_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_vet_survivalsvm_vb2_times, function(l){
  round(mean(l)/60, 2)
})

bmr_vet_survivalsvm_hyb <- list(bench.vet.ssvm.hyb, bench.vet.ssvm.ak.hyb,
                                bench.vet.ssvm.rbf.hyb)
bmr_vet_survivalsvm_hyb_times <- vector(mode = "list", length = 3L)
names(bmr_vet_survivalsvm_hyb_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_vet_survivalsvm_hyb[[bmr]])
  bmr_vet_survivalsvm_hyb_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_vet_survivalsvm_hyb_times, function(l){
  round(mean(l)/60, 2)
})

bmr_vet_survpack <- list(bench.vet.survsvm.lin, bench.vet.survsvm.gaus)
bmr_vet_survpack_times <- vector(mode = "list", length = 2L)
names(bmr_vet_survpack_times) <- c("linear", "RBF")
for(bmr in 1:2){
  models <- getBMRModels(bmr_vet_survpack[[bmr]])
  bmr_vet_survpack_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                     models[[1]][[1]][[2]][["time"]], 
                                     models[[1]][[1]][[3]][["time"]], 
                                     models[[1]][[1]][[4]][["time"]],
                                     models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_vet_survpack_times, function(l){
  round(mean(l)/60, 2)
})

# ------------------------------------------------------------------------------
# Benchmarking for the dataset leukemia with complete remission
# ------------------------------------------------------------------------------
# Read the dataset
leukemia <- read.table("https://www.stats.ox.ac.uk/pub/datasets/csb/ch14.dat",
                       na.strings = ".")

# --- Some pre-processing steps
leuk <- leukemia[, c(2, 3, 4, 5, 7, 8, 9, 10, 12, 14, 15, 16)]
names(leuk) <- c("datestudie", "treatement", "sex", "age",
                 "karnofsky", "baseline_wb", "baseline_plat",
                 "baseline_hemog", "complete_rem", "date_cr",
                 "date_last_fol_up", "status_last_fol_up")
leuk <- leuk[-42, ]
leuk$complete_rem <- as.numeric(leuk$complete_rem == "Y")
leuk$date_cr[!leuk$complete_rem] <- leuk$date_last_fol_up[!leuk$complete_rem]
# coerce dates to an appropriated format
leuk$datestudie <- sapply(leuk$datestudie, function(i){
  i <- if(nchar(i) == 5) as.character(paste("0", i, sep = ""))
  else as.character(i)
  i <- gsub("^(.{4})(.*)$", "\\119\\2",i)
  as.character(i, format = "%m%d%Y")
})
leuk$date_cr <- sapply(leuk$date_cr, function(i){
  i <- if(nchar(i) == 5) as.character(paste("0", i, sep = ""))
  else as.character(i)
  i <- gsub("^(.{4})(.*)$", "\\119\\2",i)
  as.character(i, format = "%m%d%Y")
})
leuk$date_last_fol_up <- sapply(leuk$date_last_fol_up, function(i){
  i <- if(nchar(i) == 5) as.character(paste("0", i, sep = ""))
  else as.character(i)
  i <- gsub("^(.{4})(.*)$", "\\119\\2",i)
  as.character(i, format = "%m%d%Y")
})
time.cr <- as.numeric(difftime(as.Date(leuk$date_cr, "%m%d%Y"),
                               as.Date(leuk$datestudie, "%m%d%Y"),
                               units = "days"))
time.death <- as.numeric(difftime(as.Date(leuk$date_last_fol_up, "%m%d%Y"),
                                  as.Date(leuk$datestudie, "%m%d%Y"),
                                  units = "days"))
leuk.adj <- data.frame(leuk, time.cr, time.death)
leuk.adj$datestudie <- NULL
leuk.adj$date_cr <- NULL
leuk.adj$date_last_fol_up <- NULL
leukemia.task <- makeSurvTask(data = leuk.adj,
                              target = c("time.cr", "complete_rem"))
outer <- makeResampleDesc("CV", iters = 5L)

# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.leuk.reg <- tuneWrapperGammaMu(type = "regression",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog", method = "CV",
                                         iters.rep = 5, upper = 2)
bench.leuk.cr.ssvm.reg <- benchmark(learners = tunwrp.gm.leuk.reg,
                                 tasks = leukemia.task,
                                 resamplings = outer, measures = list(c.i, lgrk,
                                                                      hr))
lapply(bench.leuk.cr.ssvm.reg$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.leuk.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog",
                                         diff.meth = "makediff3", upper = 2)
bench.leuk.cr.ssvm.vb1 <- benchmark(learners = tunwrp.gm.leuk.vb1,
                                 tasks = leukemia.task,
                                 resamplings = outer, measures = list(c.i, lgrk,
                                                                      hr))
lapply(bench.leuk.cr.ssvm.vb1$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.leuk.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog",
                                         diff.meth = "makediff3", upper = 2)
bench.leuk.cr.ssvm.vb2 <- benchmark(learners = tunwrp.gm.leuk.vb2,
                                 tasks = leukemia.task,
                                 resamplings = outer, measures = list(c.i, lgrk,
                                                                      hr))
lapply(bench.leuk.cr.ssvm.vb2$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.leuk.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog",
                                         diff.meth = "makediff3", upper = 2)
bench.leuk.cr.ssvm.hyb <- benchmark(learners = tunwrp.gm.leuk.hyb,
                                 tasks = leukemia.task,
                                 resamplings = outer, measures = list(c.i, lgrk,
                                                                      hr))
lapply(bench.leuk.cr.ssvm.hyb$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### additiv Kernel #######################################*
set.seed(123)
tunwrp.gm.leuk.ak.reg <- tuneWrapperGammaMu(type = "regression",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            method = "CV",
                                            iters.rep = 5, upper = 2)
bench.leuk.cr.ak.ssvm.reg <- benchmark(learners = tunwrp.gm.leuk.ak.reg,
                                    tasks = leukemia.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.ak.ssvm.reg$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.leuk.ak.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.leuk.cr.ak.ssvm.vb1 <- benchmark(learners = tunwrp.gm.leuk.ak.vb1,
                                    tasks = leukemia.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.ak.ssvm.vb1$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.leuk.ak.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.leuk.cr.ak.ssvm.vb2 <- benchmark(learners = tunwrp.gm.leuk.ak.vb2,
                                    tasks = leukemia.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.ak.ssvm.vb2$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.leuk.ak.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.leuk.cr.ak.ssvm.hyb <- benchmark(learners = tunwrp.gm.leuk.ak.hyb,
                                    tasks = leukemia.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.ak.ssvm.hyb$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### rbf Kernel #############################################*
set.seed(123)
tunwrp.gm.leuk.rbf.reg <- tuneWrapperGammaMu(type = "regression",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             method = "CV",
                                             iters.rep = 5, upper = 2)
bench.leuk.cr.rbf.ssvm.reg <- benchmark(learners = tunwrp.gm.leuk.rbf.reg,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.rbf.ssvm.reg$results$leuk.adj$
         surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.leuk.rbf.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             diff.meth = "makediff3", upper = 2)
bench.leuk.cr.rbf.ssvm.vb1 <- benchmark(learners = tunwrp.gm.leuk.rbf.vb1,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.rbf.ssvm.vb1$results$leuk.adj$
         surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.leuk.rbf.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             diff.meth = "makediff3", upper = 2)
bench.leuk.cr.rbf.ssvm.vb2 <- benchmark(learners = tunwrp.gm.leuk.rbf.vb2,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.rbf.ssvm.vb2$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.leuk.rbf.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             diff.meth = "makediff3", upper = 2)
bench.leuk.cr.rbf.ssvm.hyb <- benchmark(learners = tunwrp.gm.leuk.rbf.hyb,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.cr.rbf.ssvm.hyb$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# ## REFERENCE METHODS ########################################################*

# --- randomForestSRC.tuned
set.seed(123)
rfsrc.tuned <- tuneWrapperRandomForestSRC(center = TRUE, scale = TRUE,
                                          method = "CV", iters.rep = 10L)
bench.leuk.cr.rfsrc <- benchmark(rfsrc.tuned, leukemia.task, outer,
                              measures = list(ph.cindex, lgrk, hrph))
lapply(bench.leuk.cr.rfsrc$results$leuk.adj$surv.randomForestSRC.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- coxph
set.seed(123)
lrn.ph <- wrapperPHModell(center = TRUE, scale = TRUE)
bench.leuk.cr.ph <- benchmark(lrn.ph, leukemia.task,
                           resamplings = outer,
                           measures = list(ph.cindex, lgrk, hrph))
lapply(bench.leuk.cr.ph$results$leuk.adj$surv.coxph.preproc$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- Gradient Boosting
set.seed(123)
lrn.gboost <- tuneWrapperGBoost(center = TRUE, scale = TRUE)
bench.leuk.cr.gboost <- benchmark(lrn.gboost, leukemia.task,
                                  resamplings = outer,
                                  measures = list(ph.cindex, lgrk, hrph))

# --- survsvm linear
set.seed(123)
tunwrp.survsvm.linear <- tuneWrapperGammaCost(kernel = "linear",
                                              method = "CV", max.iter = 50,
                                              iters.rep = 5)
bench.lkr.survsvm.lin <- benchmark(learners = tunwrp.survsvm.linear, 
                                   tasks = leukemia.task,
                                   resamplings = outer,
                                   measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.lkr.survsvm.lin$results$leuk.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survsvm RBF
set.seed(123)
tunwrp.survsvm.gaussian <- tuneWrapperGammaCost(kernel = "gaussian",
                                                method = "CV", max.iter = 50,
                                                iters.rep = 5)
bench.lkr.survsvm.gaus <- benchmark(learners = tunwrp.survsvm.gaussian, 
                                    tasks = leukemia.task,
                                    resamplings = outer,
                                    measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.lkr.survsvm.gaus$results$leuk.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# --- Benchmark mean runtime for survivalsvm approach models
bmr_lkr_survivalsvm_reg <- list(bench.leuk.cr.ssvm.reg,
                                bench.leuk.cr.ak.ssvm.reg,
                                bench.leuk.cr.rbf.ssvm.reg)
bmr_lkr_survivalsvm_reg_times <- vector(mode = "list", length = 3L)
names(bmr_lkr_survivalsvm_reg_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkr_survivalsvm_reg[[bmr]])
  bmr_lkr_survivalsvm_reg_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkr_survivalsvm_reg_times, function(l){
  round(mean(l)/60, 2)
})

bmr_lkr_survivalsvm_vb1 <- list(bench.leuk.cr.ssvm.vb1,
                                bench.leuk.cr.ak.ssvm.vb1,
                                bench.leuk.cr.rbf.ssvm.vb1)
bmr_lkr_survivalsvm_vb1_times <- vector(mode = "list", length = 3L)
names(bmr_lkr_survivalsvm_vb1_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkr_survivalsvm_vb1[[bmr]])
  bmr_lkr_survivalsvm_vb1_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkr_survivalsvm_vb1_times, function(l){
  round(mean(l)/60, 2)
})


bmr_lkr_survivalsvm_vb2 <- list(bench.leuk.cr.ssvm.vb2,
                                bench.leuk.cr.ak.ssvm.vb2, 
                                bench.leuk.cr.rbf.ssvm.vb2)
bmr_lkr_survivalsvm_vb2_times <- vector(mode = "list", length = 3L)
names(bmr_lkr_survivalsvm_vb2_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkr_survivalsvm_vb2[[bmr]])
  bmr_lkr_survivalsvm_vb2_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkr_survivalsvm_vb2_times, function(l){
  round(mean(l)/60, 2)
})

bmr_lkr_survivalsvm_hyb <- list(bench.leuk.cr.ssvm.hyb,
                                bench.leuk.cr.ak.ssvm.hyb,
                                bench.leuk.cr.rbf.ssvm.hyb)
bmr_lkr_survivalsvm_hyb_times <- vector(mode = "list", length = 3L)
names(bmr_lkr_survivalsvm_hyb_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkr_survivalsvm_hyb[[bmr]])
  bmr_lkr_survivalsvm_hyb_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkr_survivalsvm_hyb_times, function(l){
  round(mean(l)/60, 2)
})

bmr_lkr_survpack <- list(bench.lkr.survsvm.lin, bench.lkr.survsvm.gaus)
bmr_lkr_survpack_times <- vector(mode = "list", length = 2L)
names(bmr_lkr_survpack_times) <- c("linear", "RBF")
for(bmr in 1:2){
  models <- getBMRModels(bmr_lkr_survpack[[bmr]])
  bmr_lkr_survpack_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                     models[[1]][[1]][[2]][["time"]], 
                                     models[[1]][[1]][[3]][["time"]], 
                                     models[[1]][[1]][[4]][["time"]],
                                     models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkr_survpack_times, function(l){
  round(mean(l)/60, 2)
})

# ------------------------------------------------------------------------------
# Benchmarking for the dataset leukemia with death
# ------------------------------------------------------------------------------
leukemia.task <- makeSurvTask(data = leuk.adj, target = c("time.death",
                                                          "status_last_fol_up"))
outer <- makeResampleDesc("CV", iters = 5L)

# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.leuk.reg <- tuneWrapperGammaMu(type = "regression",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog", method = "CV",
                                         iters.rep = 5, upper = 2)
bench.leuk.death.ssvm.reg <- benchmark(learners = tunwrp.gm.leuk.reg,
                                 tasks = leukemia.task,
                                 resamplings = outer, measures = list(c.i, lgrk,
                                                                      hr))
lapply(bench.leuk.death.ssvm.reg$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.leuk.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog",
                                         diff.meth = "makediff3", upper = 2)
bench.leuk.death.ssvm.vb1 <- benchmark(learners = tunwrp.gm.leuk.vb1,
                                 tasks = leukemia.task,
                                 resamplings = outer,
                                 measures = list(c.i, lgrk, hr))
lapply(bench.leuk.death.ssvm.vb1$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.leuk.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog",
                                         diff.meth = "makediff3", upper = 2)
bench.leuk.death.ssvm.vb2 <- benchmark(learners = tunwrp.gm.leuk.vb2,
                                 tasks = leukemia.task,
                                 resamplings = outer, measures = list(c.i, lgrk,
                                                                      hr))
lapply(bench.leuk.death.ssvm.vb2$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.leuk.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                         kernel = "lin_kernel",
                                         opt.meth = "quadprog",
                                         diff.meth = "makediff3", upper = 2)
bench.leuk.death.ssvm.hyb <- benchmark(learners = tunwrp.gm.leuk.hyb,
                                 tasks = leukemia.task,
                                 resamplings = outer, measures = list(c.i, lgrk,
                                                                      hr))
lapply(bench.leuk.death.ssvm.hyb$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### additiv Kernel #########################################*
set.seed(123)
tunwrp.gm.leuk.ak.reg <- tuneWrapperGammaMu(type = "regression",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            method = "CV",
                                            iters.rep = 5, upper = 2)
bench.leuk.death.ak.ssvm.reg <- benchmark(learners = tunwrp.gm.leuk.ak.reg,
                                    tasks = leukemia.task,
                                    resamplings = outer, measures = list(c.i,
                                                                         lgrk,
                                                                         hr))
lapply(bench.leuk.death.ak.ssvm.reg$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.leuk.ak.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.leuk.death.ak.ssvm.vb1 <- benchmark(learners = tunwrp.gm.leuk.ak.vb1,
                                    tasks = leukemia.task,
                                    resamplings = outer, measures = list(c.i,
                                                                         lgrk,
                                                                         hr))
lapply(bench.leuk.death.ak.ssvm.vb1$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.leuk.ak.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.leuk.death.ak.ssvm.vb2 <- benchmark(learners = tunwrp.gm.leuk.ak.vb2,
                                    tasks = leukemia.task,
                                    resamplings = outer, measures = list(c.i,
                                                                         lgrk,
                                                                         hr))
lapply(bench.leuk.death.ak.ssvm.vb2$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.leuk.ak.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                            kernel = "add_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.leuk.death.ak.ssvm.hyb <- benchmark(learners = tunwrp.gm.leuk.ak.hyb,
                                    tasks = leukemia.task,
                                    resamplings = outer, measures = list(c.i,
                                                                         lgrk,
                                                                         hr))
lapply(bench.leuk.death.ak.ssvm.hyb$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### rbf Kernel #############################################*
set.seed(123)
tunwrp.gm.leuk.rbf.reg <- tuneWrapperGammaMu(type = "regression",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             method = "CV",
                                             iters.rep = 5, upper = 2)
bench.leuk.death.rbf.ssvm.reg <- benchmark(learners = tunwrp.gm.leuk.rbf.reg,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.death.rbf.ssvm.reg$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.leuk.rbf.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             diff.meth = "makediff3", upper = 2)
bench.leuk.death.rbf.ssvm.vb1 <- benchmark(learners = tunwrp.gm.leuk.rbf.vb1,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.death.rbf.ssvm.vb1$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.leuk.rbf.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             diff.meth = "makediff3", upper = 2)
bench.leuk.death.rbf.ssvm.vb2 <- benchmark(learners = tunwrp.gm.leuk.rbf.vb2,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.death.rbf.ssvm.vb2$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.leuk.rbf.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                             kernel = "rbf_kernel",
                                             opt.meth = "quadprog",
                                             diff.meth = "makediff3", upper = 2)
bench.leuk.death.rbf.ssvm.hyb <- benchmark(learners = tunwrp.gm.leuk.rbf.hyb,
                                     tasks = leukemia.task,
                                     resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.leuk.death.rbf.ssvm.hyb$results$leuk.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# ## REFERENCE METHODS ########################################################*

# --- randomForestSRC.tuned
set.seed(123)
rfsrc.tuned <- tuneWrapperRandomForestSRC(center = TRUE, scale = TRUE,
                                          method = "CV", iters.rep = 10L)
bench.leuk.death.rfsrc <- benchmark(rfsrc.tuned, leukemia.task, outer,
                              measures = list(ph.cindex, lgrk, hrph))
lapply(bench.leuk.death.rfsrc$results$leuk.adj$surv.randomForestSRC.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- coxph
set.seed(123)
lrn.ph <- wrapperPHModell(center = TRUE, scale = TRUE)
bench.leuk.death.ph <- benchmark(lrn.ph, leukemia.task, resamplings = outer,
                           measures = list(ph.cindex, lgrk, hrph))
lapply(bench.leuk.death.ph$results$leuk.adj$surv.coxph.preproc$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- Gradient boosting
set.seed(123)
lrn.gboost <- tuneWrapperGBoost(center = TRUE, scale = TRUE)
bench.leuk.death.gboost <- benchmark(lrn.gboost, leukemia.task,
                                     resamplings = outer,
                                     measures = list(ph.cindex, lgrk, hrph))

# --- survsvm linear
set.seed(123)
tunwrp.survsvm.linear <- tuneWrapperGammaCost(kernel = "linear",
                                              method = "CV", max.iter = 50,
                                              iters.rep = 5)
bench.lkd.survsvm.lin <- benchmark(learners = tunwrp.survsvm.linear, 
                                   tasks = leukemia.task,
                                   resamplings = outer,
                                   measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.lkd.survsvm.lin$results$leuk.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survsvm RBF
set.seed(123)
tunwrp.survsvm.gaussian <- tuneWrapperGammaCost(kernel = "gaussian",
                                                method = "CV", max.iter = 50,
                                                iters.rep = 5)
bench.lkd.survsvm.gaus <- benchmark(learners = tunwrp.survsvm.gaussian, 
                                    tasks = leukemia.task,
                                    resamplings = outer,
                                    measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.lkd.survsvm.gaus$results$leuk.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

##--- Benchmark mean runtime for survivalsvm approach models
bmr_lkd_survivalsvm_reg <- list(bench.leuk.death.ssvm.reg,
                                bench.leuk.death.ak.ssvm.reg,
                                bench.leuk.death.rbf.ssvm.reg)
bmr_lkd_survivalsvm_reg_times <- vector(mode = "list", length = 3L)
names(bmr_lkd_survivalsvm_reg_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkd_survivalsvm_reg[[bmr]])
  bmr_lkd_survivalsvm_reg_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkd_survivalsvm_reg_times, function(l){
  round(mean(l)/60, 2)
})

bmr_lkd_survivalsvm_vb1 <- list(bench.leuk.death.ssvm.vb1,
                                bench.leuk.death.ak.ssvm.vb1,
                                bench.leuk.death.rbf.ssvm.vb1)
bmr_lkd_survivalsvm_vb1_times <- vector(mode = "list", length = 3L)
names(bmr_lkd_survivalsvm_vb1_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkd_survivalsvm_vb1[[bmr]])
  bmr_lkd_survivalsvm_vb1_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkd_survivalsvm_vb1_times, function(l){
  round(mean(l)/60, 2)
})

bmr_lkd_survivalsvm_vb2 <- list(bench.leuk.death.ssvm.vb2,
                                bench.leuk.death.ak.ssvm.vb2, 
                                bench.leuk.death.rbf.ssvm.vb2)
bmr_lkd_survivalsvm_vb2_times <- vector(mode = "list", length = 3L)
names(bmr_lkd_survivalsvm_vb2_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkd_survivalsvm_vb2[[bmr]])
  bmr_lkd_survivalsvm_vb2_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkd_survivalsvm_vb2_times, function(l){
  round(mean(l)/60, 2)
})

bmr_lkd_survivalsvm_hyb <- list(bench.leuk.death.ssvm.hyb,
                                bench.leuk.death.ak.ssvm.hyb,
                                bench.leuk.death.rbf.ssvm.hyb)
bmr_lkd_survivalsvm_hyb_times <- vector(mode = "list", length = 3L)
names(bmr_lkd_survivalsvm_hyb_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_lkd_survivalsvm_hyb[[bmr]])
  bmr_lkd_survivalsvm_hyb_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkd_survivalsvm_hyb_times, function(l){
  round(mean(l)/60, 2)
})

bmr_lkd_survpack <- list(bench.lkd.survsvm.lin, bench.lkd.survsvm.gaus)
bmr_lkd_survpack_times <- vector(mode = "list", length = 2L)
names(bmr_lkd_survpack_times) <- c("linear", "RBF")
for(bmr in 1:2){
  models <- getBMRModels(bmr_lkd_survpack[[bmr]])
  bmr_lkd_survpack_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                     models[[1]][[1]][[2]][["time"]], 
                                     models[[1]][[1]][[3]][["time"]], 
                                     models[[1]][[1]][[4]][["time"]],
                                     models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_lkd_survpack_times, function(l){
  round(mean(l)/60, 2)
})


# ------------------------------------------------------------------------------
# Benchmarking for the dataset GBSG2
# ------------------------------------------------------------------------------
require(pec)
data(GBSG2, package = "pec")
gbcsg.adj <- GBSG2
gbcsg.task <- makeSurvTask(data = gbcsg.adj, target = c("time", "cens"))
outer <- makeResampleDesc("CV", iters = 5L)

# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.reg <- tuneWrapperGammaMu(type = "regression",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog", method = "CV",
                                    iters.rep = 5, tune.scale = FALSE,
                                    lower = -5, upper = 5)
bench.gbcsg.ssvm.reg <- benchmark(learners = tunwrp.gm.reg, tasks = gbcsg.task,
                                  resamplings = outer, measures = list(c.i,
                                                                       lgrk,
                                                                       hr))
lapply(bench.gbcsg.ssvm.reg$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog",
                                    diff.meth = "makediff3", method = "CV",
                                    iters.rep = 5, tune.scale = FALSE,
                                    lower = -5, upper = 5)
bench.gbcsg.ssvm.vb1 <- benchmark(learners = tunwrp.gm.vb1, tasks = gbcsg.task,
                                  resamplings = outer, measures = list(c.i,
                                                                       lgrk,
                                                                       hr))
lapply(bench.gbcsg.ssvm.vb1$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog",
                                    diff.meth = "makediff3", method = "CV",
                                    iters.rep = 5, tune.scale = FALSE,
                                    lower = -10, upper = 5)
bench.gbcsg.ssvm.vb2 <- benchmark(learners = tunwrp.gm.vb2, tasks = gbcsg.task,
                                  resamplings = outer, measures = list(c.i,
                                                                       lgrk,
                                                                       hr))
lapply(bench.gbcsg.ssvm.vb2$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                    kernel = "lin_kernel",
                                    opt.meth = "quadprog",
                                    diff.meth = "makediff3", method = "CV",
                                    iters.rep = 5, tune.scale = FALSE,
                                    lower = -5, upper = 5)
bench.gbcsg.ssvm.hyb <- benchmark(learners = tunwrp.gm.hyb, tasks = gbcsg.task,
                                  resamplings = outer, measures = list(c.i,
                                                                       lgrk))
lapply(bench.gbcsg.ssvm.hyb$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(na.omit(i))
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### additiv Kernel #########################################*
# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.ak.reg <- tuneWrapperGammaMu(type = "regression",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog", method = "CV",
                                       iters.rep = 5, tune.scale = FALSE,
                                       lower = -5, upper = 5)
bench.bgcsg.ssvm.ak.reg <- benchmark(learners = tunwrp.gm.ak.reg,
                                     tasks = gbcsg.task, resamplings = outer, 
                                     measures = list(c.i, lgrk, hr))
lapply(bench.bgcsg.ssvm.ak.reg$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.ak.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog",
                                       diff.meth = "makediff3", method = "CV",
                                       iters.rep = 5, tune.scale = FALSE,
                                       lower = -5, upper = 5)
bench.gbcsg.ssvm.ak.vb1 <- benchmark(learners = tunwrp.gm.ak.vb1,
                                     tasks = gbcsg.task, resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.gbcsg.ssvm.ak.vb1$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.ak.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog",
                                       diff.meth = "makediff3", method = "CV",
                                       iters.rep = 5, tune.scale = FALSE,
                                       lower = -5, upper = 5)
bench.gbcsg.ssvm.ak.vb2 <- benchmark(learners = tunwrp.gm.ak.vb2,
                                     tasks = gbcsg.task, resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.gbcsg.ssvm.ak.vb2$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.Ak.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                       kernel = "add_kernel",
                                       opt.meth = "quadprog",
                                       diff.meth = "makediff3", method = "CV",
                                       iters.rep = 5, tune.scale = FALSE,
                                       lower = -5, upper = 5, resolution = 5L)
bench.gbcsg.ssvm.ak.hyb <- benchmark(learners = tunwrp.gm.Ak.hyb,
                                     tasks = gbcsg.task, resamplings = outer,
                                     measures = list(c.i, lgrk, hr))
lapply(bench.gbcsg.ssvm.Ak.hyb$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### rbf Kernel #############################################*
# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.rbf.reg <- tuneWrapperGammaMu(type = "regression",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog", method = "CV",
                                        iters.rep = 5, tune.scale = FALSE,
                                        lower = -5, upper = 5)
bench.bgcsg.ssvm.rbf.reg <- benchmark(learners = tunwrp.gm.rbf.reg,
                                      tasks = gbcsg.task, resamplings = outer, 
                                      measures = list(c.i, lgrk, hr))
lapply(bench.bgcsg.ssvm.rbf.reg$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.rbf.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3", method = "CV",
                                        iters.rep = 5, tune.scale = FALSE,
                                        lower = -5, upper = 5)
bench.gbcsg.ssvm.rbf.vb1 <- benchmark(learners = tunwrp.gm.rbf.vb1,
                                      tasks = gbcsg.task, resamplings = outer,
                                      measures = list(c.i, lgrk, hr))
lapply(bench.gbcsg.ssvm.rbf.vb1$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.rbf.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3", method = "CV",
                                        iters.rep = 5, tune.scale = FALSE,
                                        lower = -5, upper = 5)
bench.gbcsg.ssvm.rbf.vb2 <- benchmark(learners = tunwrp.gm.rbf.vb2,
                                      tasks = gbcsg.task, resamplings = outer,
                                      measures = list(c.i, lgrk, hr))
lapply(bench.gbcsg.ssvm.rbf.vb2$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.rbf.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                        kernel = "rbf_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3", method = "CV",
                                        iters.rep = 5, tune.scale = FALSE,
                                        lower = -5, upper = 5, resolution = 5L)
bench.gbcsg.ssvm.rbf.hyb <- benchmark(learners = tunwrp.gm.rbf.hyb,
                                      tasks = gbcsg.task, resamplings = outer,
                                      measures = list(c.i, lgrk, hr))
lapply(bench.gbcsg.ssvm.rbf.hyb$results$gbcsg.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# ## REFERENCE METHODS ########################################################*

# --- randomForestSRC.tuned
set.seed(123)
rfsrc.tuned <- tuneWrapperRandomForestSRC(center = TRUE, scale = TRUE,
                                          method = "CV", iters.rep = 10L)
bench.gbcsg.rfsrc <- benchmark(rfsrc.tuned, gbcsg.task, outer,
                               measures = list(ph.cindex, lgrk, hrph))
lapply(bench.gbcsg.rfsrc$results$gbcsg.adj$surv.randomForestSRC.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- coxph
set.seed(123)
lrn.ph <- wrapperPHModell(center = TRUE, scale = TRUE)
bench.gbcsg.ph <- benchmark(lrn.ph, gbcsg.task, resamplings = outer,
                            measures = list(ph.cindex, lgrk, hrph))
lapply(bench.gbcsg.ph$results$gbcsg.adj$surv.coxph.preproc$measures.test,
       function(i){
         a <- CI(i)
         err <- a[1] - a[2]
         mittelwert <- a[2]
         r <- round(c(mittelwert, err), 2)
         names(r) <- c("mean", "error")
         return(r)
       })
# --- Gradient Boosting
set.seed(123)
lrn.gboost <- tuneWrapperGBoost(center = TRUE, scale = TRUE)
bench.gbcsg.gboost <- benchmark(lrn.gboost, gbcsg.task, resamplings = outer,
                                measures = list(ph.cindex, lgrk, hrph))

# --- survsvm linear
set.seed(123)
tunwrp.survsvm.linear <- tuneWrapperGammaCost(kernel = "linear",
                                              method = "CV", max.iter = 50,
                                              iters.rep = 5)
bench.gbcsg.survsvm.lin <- benchmark(learners = tunwrp.survsvm.linear, 
                                     tasks = gbcsg.task,
                                     resamplings = outer,
                                     measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.gbcsg.survsvm.lin$results$gbcsg.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survsvm RBF
set.seed(123)
tunwrp.survsvm.gaussian <- tuneWrapperGammaCost(kernel = "gaussian",
                                                method = "CV", max.iter = 50,
                                                iters.rep = 5)
bench.gbcsg.survsvm.gaus <- benchmark(learners = tunwrp.survsvm.gaussian, 
                                      tasks = gbcsg.task,
                                      resamplings = outer,
                                      measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.gbcsg.survsvm.gaus$results$gbcsg.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# --- Benchmark mean runtime for survivalsvm approach models
bmr_gbcsg_survivalsvm_reg <- list(bench.gbcsg.ssvm.reg, bench.bgcsg.ssvm.ak.reg,
                                  bench.bgcsg.ssvm.rbf.reg)
bmr_gbcsg_survivalsvm_reg_times <- vector(mode = "list", length = 3L)
names(bmr_gbcsg_survivalsvm_reg_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_gbcsg_survivalsvm_reg[[bmr]])
  bmr_gbcsg_survivalsvm_reg_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                              models[[1]][[1]][[2]][["time"]], 
                                              models[[1]][[1]][[3]][["time"]], 
                                              models[[1]][[1]][[4]][["time"]],
                                              models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_gbcsg_survivalsvm_reg_times, function(l){
  round(mean(l)/60, 2)
})

bmr_gbcsg_survivalsvm_vb1 <- list(bench.gbcsg.ssvm.vb1, bench.gbcsg.ssvm.ak.vb1,
                                  bench.gbcsg.ssvm.rbf.vb1)
bmr_gbcsg_survivalsvm_vb1_times <- vector(mode = "list", length = 3L)
names(bmr_gbcsg_survivalsvm_vb1_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_gbcsg_survivalsvm_vb1[[bmr]])
  bmr_gbcsg_survivalsvm_vb1_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                              models[[1]][[1]][[2]][["time"]], 
                                              models[[1]][[1]][[3]][["time"]], 
                                              models[[1]][[1]][[4]][["time"]],
                                              models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_gbcsg_survivalsvm_vb1_times, function(l){
  round(mean(l)/60, 2)
})
bmr_gbcsg_survivalsvm_vb2 <- list(bench.gbcsg.ssvm.vb2, bench.gbcsg.ssvm.ak.vb2, 
                                  bench.gbcsg.ssvm.rbf.vb2)
bmr_gbcsg_survivalsvm_vb2_times <- vector(mode = "list", length = 3L)
names(bmr_gbcsg_survivalsvm_vb2_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_gbcsg_survivalsvm_vb2[[bmr]])
  bmr_gbcsg_survivalsvm_vb2_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                              models[[1]][[1]][[2]][["time"]], 
                                              models[[1]][[1]][[3]][["time"]], 
                                              models[[1]][[1]][[4]][["time"]],
                                              models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_gbcsg_survivalsvm_vb2_times, function(l){
  round(mean(l)/60, 2)
})
# --- bench.gbcsg.ssvm.rbf.hyb
bmr_gbcsg_survivalsvm_hyb <- list(bench.gbcsg.ssvm.hyb,
                                  bench.gbcsg.ssvm.ak.hyb,
                                  bench.gbcsg.ssvm.rbf.hyb)
bmr_gbcsg_survivalsvm_hyb_times <- vector(mode = "list", length = 3L)
names(bmr_gbcsg_survivalsvm_hyb_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_gbcsg_survivalsvm_hyb[[bmr]])
  bmr_gbcsg_survivalsvm_hyb_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                              models[[1]][[1]][[2]][["time"]], 
                                              models[[1]][[1]][[3]][["time"]], 
                                              models[[1]][[1]][[4]][["time"]],
                                              models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_gbcsg_survivalsvm_hyb_times, function(l){
  round(mean(l)/60, 2)
})

# --- bench.gbcsg.survsvm.gaus
bmr_gbcsg_survpack <- list(bench.gbcsg.survsvm.lin, NULL)
bmr_gbcsg_survpack_times <- vector(mode = "list", length = 2L)
names(bmr_gbcsg_survpack_times) <- c("linear", "RBF")
for(bmr in 1:2){
  models <- getBMRModels(bmr_gbcsg_survpack[[bmr]])
  bmr_gbcsg_survpack_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                       models[[1]][[1]][[2]][["time"]], 
                                       models[[1]][[1]][[3]][["time"]], 
                                       models[[1]][[1]][[4]][["time"]],
                                       models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_gbcsg_survpack_times, function(l){
  round(mean(l)/60, 2)
})

# ------------------------------------------------------------------------------
# Benchmarking for the dataset MCLC
# ------------------------------------------------------------------------------
# --- Some pre-processing steps
data("lung", package = "survival")
mlc <- lung
mlc$inst <- NULL
mlc.adj <- na.omit(mlc)
mlc.adj$status <- as.numeric(mlc.adj$status == 2)
mlc.adj$sex <- as.factor(mlc.adj$sex)

mlc.task <- makeSurvTask(data = mlc.adj, target = c("time", "status"))
outer <- makeResampleDesc("CV", iters = 5L)
# --- survivalsvm regression.tuned
set.seed(123)
tunwrp.gm.mlc.reg <- tuneWrapperGammaMu(type = "regression",
                                        kernel = "lin_kernel",
                                        opt.meth = "quadprog", method = "CV",
                                        iters.rep = 5, upper = 2)
bench.mlc.ssvm.reg <- benchmark(learners = tunwrp.gm.mlc.reg, tasks = mlc.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
lapply(bench.mlc.ssvm.reg$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.mlc.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                        kernel = "lin_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3", upper = 2)
bench.mlc.ssvm.vb1 <- benchmark(learners = tunwrp.gm.mlc.vb1, tasks = mlc.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
lapply(bench.mlc.ssvm.vb1$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.mlc.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                        kernel = "lin_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3", upper = 2)
bench.mlc.ssvm.vb2 <- benchmark(learners = tunwrp.gm.mlc.vb2, tasks = mlc.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
lapply(bench.mlc.ssvm.vb2$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.mlc.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                        kernel = "lin_kernel",
                                        opt.meth = "quadprog",
                                        diff.meth = "makediff3", upper = 2)
bench.mlc.ssvm.hyb <- benchmark(learners = tunwrp.gm.mlc.hyb, tasks = mlc.task,
                                resamplings = outer, measures = list(c.i, lgrk,
                                                                     hr))
lapply(bench.mlc.ssvm.hyb$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### additiv  Kernel ########################################*
set.seed(123)
tunwrp.gm.mlc.ak.reg <- tuneWrapperGammaMu(type = "regression",
                                           kernel = "add_kernel",
                                           opt.meth = "quadprog", method = "CV",
                                           iters.rep = 5, upper = 2)
bench.mlc.ak.ssvm.reg <- benchmark(learners = tunwrp.gm.mlc.ak.reg,
                                   tasks = mlc.task, resamplings = outer,
                                   measures = list(c.i, lgrk, hr))
lapply(bench.mlc.ak.ssvm.reg$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.mlc.ak.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                           kernel = "add_kernel",
                                           opt.meth = "quadprog",
                                           diff.meth = "makediff3", upper = 2)
bench.mlc.ak.ssvm.vb1 <- benchmark(learners = tunwrp.gm.mlc.ak.vb1,
                                   tasks = mlc.task, resamplings = outer,
                                   measures = list(c.i, lgrk, hr))
lapply(bench.mlc.ak.ssvm.vb1$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.mlc.ak.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                           kernel = "add_kernel",
                                           opt.meth = "quadprog",
                                           diff.meth = "makediff3", upper = 2)
bench.mlc.ak.ssvm.vb2 <- benchmark(learners = tunwrp.gm.mlc.ak.vb2,
                                   tasks = mlc.task, resamplings = outer,
                                   measures = list(c.i, lgrk, hr))
lapply(bench.mlc.ak.ssvm.vb2$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.mlc.ak.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                           kernel = "add_kernel",
                                           opt.meth = "quadprog",
                                           diff.meth = "makediff3",
                                           upper = 2, tune.scale = FALSE)
bench.mlc.ak.ssvm.hyb <- benchmark(learners = tunwrp.gm.mlc.ak.hyb,
                                   tasks = mlc.task, resamplings = outer,
                                   measures = list(c.i, lgrk, hr))
lapply(bench.mlc.ak.ssvm.hyb$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# #################### rbf Kernel #############################################*
set.seed(123)
tunwrp.gm.mlc.rbf.reg <- tuneWrapperGammaMu(type = "regression",
                                            kernel = "rbf_kernel",
                                            opt.meth = "quadprog",
                                            method = "CV",
                                            iters.rep = 5, upper = 2)
bench.mlc.rbf.ssvm.reg <- benchmark(learners = tunwrp.gm.mlc.rbf.reg,
                                    tasks = mlc.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.mlc.rbf.ssvm.reg$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle1.tuned
set.seed(123)
tunwrp.gm.mlc.rbf.vb1 <- tuneWrapperGammaMu(type = "vanbelle1",
                                            kernel = "rbf_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.mlc.rbf.ssvm.vb1 <- benchmark(learners = tunwrp.gm.mlc.rbf.vb1,
                                    tasks = mlc.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.mlc.rbf.ssvm.vb1$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm vanbelle2.tuned
set.seed(123)
tunwrp.gm.mlc.rbf.vb2 <- tuneWrapperGammaMu(type = "vanbelle2",
                                            kernel = "rbf_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.mlc.rbf.ssvm.vb2 <- benchmark(learners = tunwrp.gm.mlc.rbf.vb2,
                                    tasks = mlc.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.mlc.rbf.ssvm.vb2$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survivalsvm hybrid.tuned
set.seed(123)
tunwrp.gm.mlc.rbf.hyb <- tuneWrapperGammaMu(type = "hybrid",
                                            kernel = "rbf_kernel",
                                            opt.meth = "quadprog",
                                            diff.meth = "makediff3", upper = 2)
bench.mlc.rbf.ssvm.hyb <- benchmark(learners = tunwrp.gm.mlc.rbf.hyb,
                                    tasks = mlc.task,
                                    resamplings = outer,
                                    measures = list(c.i, lgrk, hr))
lapply(bench.mlc.rbf.ssvm.hyb$results$mlc.adj$surv.survivalsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# ## REFERENCE METHODS ########################################################*

# --- randomForestSRC.tuned
set.seed(123)
rfsrc.tuned <- tuneWrapperRandomForestSRC(center = TRUE, scale = TRUE,
                                          method = "CV", iters.rep = 10L)
bench.mlc.rfsrc <- benchmark(rfsrc.tuned, mlc.task, outer,
                             measures = list(ph.cindex, lgrk, hrph))
lapply(bench.mlc.rfsrc$results$mlc.adj$surv.randomForestSRC.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- coxph
set.seed(123)
lrn.ph <- wrapperPHModell(center = TRUE, scale = TRUE)
bench.mlc.ph <- benchmark(lrn.ph, mlc.task, resamplings = outer,
                          measures = list(ph.cindex, lgrk, hrph))
lapply(bench.mlc.ph$results$mlc.adj$surv.coxph.preproc$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- Gradient boosting
set.seed(123)
lrn.gboost <- tuneWrapperGBoost(center = TRUE, scale = TRUE)
bench.mlc.gboost <- benchmark(lrn.gboost, mlc.task, resamplings = outer,
                              measures = list(ph.cindex, lgrk, hrph))


# --- survsvm linear
set.seed(123)
tunwrp.survsvm.linear <- tuneWrapperGammaCost(kernel = "linear",
                                              method = "CV", max.iter = 50,
                                              iters.rep = 5)
bench.mlc.survsvm.lin <- benchmark(learners = tunwrp.survsvm.linear, 
                                   tasks = mlc.task,
                                   resamplings = outer,
                                   measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.mlc.survsvm.lin$results$mlc.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })
# --- survsvm RBF
set.seed(123)
tunwrp.survsvm.gaussian <- tuneWrapperGammaCost(kernel = "gaussian",
                                                method = "CV", max.iter = 50,
                                                iters.rep = 5)
bench.mlc.survsvm.gaus <- benchmark(learners = tunwrp.survsvm.gaussian, 
                                    tasks = mlc.task,
                                    resamplings = outer,
                                    measures = list(c.i.survsvm, lgrk, hr))
lapply(bench.mlc.survsvm.gaus$results$mlc.adj$surv.survsvm.preproc.tuned$
         measures.test, function(i){
           a <- CI(i)
           err <- a[1] - a[2]
           mittelwert <- a[2]
           r <- round(c(mittelwert, err), 2)
           names(r) <- c("mean", "error")
           return(r)
         })

# --- Benchmark mean runtime for survivalsvm approach models
bmr_mlc_survivalsvm_reg <- list(bench.mlc.ssvm.reg, bench.mlc.ak.ssvm.reg,
                                bench.mlc.rbf.ssvm.reg)
bmr_mlc_survivalsvm_reg_times <- vector(mode = "list", length = 3L)
names(bmr_mlc_survivalsvm_reg_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_mlc_survivalsvm_reg[[bmr]])
  bmr_mlc_survivalsvm_reg_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_mlc_survivalsvm_reg_times, function(l){
  round(mean(l)/60, 2)
})

bmr_mlc_survivalsvm_vb1 <- list(bench.mlc.ssvm.vb1, bench.mlc.ak.ssvm.vb1,
                                bench.mlc.rbf.ssvm.vb1)
bmr_mlc_survivalsvm_vb1_times <- vector(mode = "list", length = 3L)
names(bmr_mlc_survivalsvm_vb1_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_mlc_survivalsvm_vb1[[bmr]])
  bmr_mlc_survivalsvm_vb1_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_mlc_survivalsvm_vb1_times, function(l){
  round(mean(l)/60, 2)
})


bmr_mlc_survivalsvm_vb2 <- list(bench.mlc.ssvm.vb2, bench.mlc.ak.ssvm.vb2, 
                                bench.mlc.rbf.ssvm.vb2)
bmr_mlc_survivalsvm_vb2_times <- vector(mode = "list", length = 3L)
names(bmr_mlc_survivalsvm_vb2_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_mlc_survivalsvm_vb2[[bmr]])
  bmr_mlc_survivalsvm_vb2_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_mlc_survivalsvm_vb2_times, function(l){
  round(mean(l)/60, 2)
})

bmr_mlc_survivalsvm_hyb <- list(bench.mlc.ssvm.hyb, bench.mlc.ak.ssvm.hyb,
                                bench.mlc.rbf.ssvm.hyb)
bmr_mlc_survivalsvm_hyb_times <- vector(mode = "list", length = 3L)
names(bmr_mlc_survivalsvm_hyb_times) <- c("linear", "additive", "RBF")
for(bmr in 1:3){
  models <- getBMRModels(bmr_mlc_survivalsvm_hyb[[bmr]])
  bmr_mlc_survivalsvm_hyb_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                            models[[1]][[1]][[2]][["time"]], 
                                            models[[1]][[1]][[3]][["time"]], 
                                            models[[1]][[1]][[4]][["time"]],
                                            models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_mlc_survivalsvm_hyb_times, function(l){
  round(mean(l)/60, 2)
})

bmr_mlc_survpack <- list(bench.mlc.survsvm.lin, bench.mlc.survsvm.gaus)
bmr_mlc_survpack_times <- vector(mode = "list", length = 2L)
names(bmr_mlc_survpack_times) <- c("linear", "RBF")
for(bmr in 1:2){
  models <- getBMRModels(bmr_mlc_survpack[[bmr]])
  bmr_mlc_survpack_times[[bmr]] <- c(models[[1]][[1]][[1]][["time"]], 
                                     models[[1]][[1]][[2]][["time"]], 
                                     models[[1]][[1]][[3]][["time"]], 
                                     models[[1]][[1]][[4]][["time"]],
                                     models[[1]][[1]][[5]][["time"]])
}
lapply(bmr_mlc_survpack_times, function(l){
  round(mean(l)/60, 2)
})

# ------------------------------------------------------------------------------
# Prepare data.frame of measures to generate the  box plots
# ------------------------------------------------------------------------------
bench.vet.ci <- data.frame(c(getBMRPerformances(bench.vet.ssvm.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.ak.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.rbf.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.ak.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.rbf.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.ak.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.rbf.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.ak.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ssvm.rbf.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.survsvm.lin,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.survsvm.gaus,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.ph,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.vet.rfsrc,
                                                as.df = TRUE)[, "ci"], 
                             getBMRPerformances(bench.vet.gboost,
                                                as.df = TRUE)[, "ci"]),
                           c(rep("vanbelle1.lin", 5),
                             rep("vanbelle1.add", 5),
                             rep("vanbelle1.RBF", 5),
                             rep("vanbelle2.lin", 5),
                             rep("vanbelle2.add", 5),
                             rep("vanbelle2.RBF", 5),
                             rep("SSVR.lin", 5),
                             rep("SSVR.add", 5),
                             rep("SSVR.RBF", 5),
                             rep("hybrid.lin", 5),
                             rep("hybrid.add", 5),
                             rep("hybrid.RBF", 5),
                             rep("evers.lin", 5),
                             rep("evers.RBF", 5),
                             rep("PH", 5), rep("RSF", 5), rep("GBoost", 5))
)
names(bench.vet.ci) <- c("CIndex", "Model")
bench.vet.ci$Model <- factor(bench.vet.ci$Model, levels = c("vanbelle1.lin",
                                                            "vanbelle1.add",
                                                            "vanbelle1.RBF",
                                                            "vanbelle2.lin",
                                                            "vanbelle2.add",
                                                            "vanbelle2.RBF",
                                                            "SSVR.lin",
                                                            "SSVR.add",
                                                            "SSVR.RBF",
                                                            "hybrid.lin",
                                                            "hybrid.add",
                                                            "hybrid.RBF",
                                                            "evers.lin",
                                                            "evers.RBF",
                                                            "PH", "RSF",
                                                            "GBoost"),
                             ordered = TRUE)
ddply(bench.vet.ci, ~Model, summarise, mean=round(mean(CIndex), 2),
      sd= round(sd(CIndex), 2))

# ## Leukemia complete remission Benchmark objects ############################*
bench.lkr.ci <- data.frame(c(getBMRPerformances(bench.leuk.cr.ssvm.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ak.ssvm.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.rbf.ssvm.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ssvm.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ak.ssvm.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.rbf.ssvm.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ssvm.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ak.ssvm.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.rbf.ssvm.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ssvm.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ak.ssvm.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.rbf.ssvm.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.lkr.survsvm.lin,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.lkr.survsvm.gaus,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.ph,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.leuk.cr.rfsrc,
                                                as.df = TRUE)[, "ci"], 
                             getBMRPerformances(bench.leuk.cr.gboost,
                                                as.df = TRUE)[, "ci"]),
                           c(rep("vanbelle1.lin", 5),
                             rep("vanbelle1.add", 5),
                             rep("vanbelle1.RBF", 5),
                             rep("vanbelle2.lin", 5),
                             rep("vanbelle2.add", 5),
                             rep("vanbelle2.RBF", 5),
                             rep("SSVR.lin", 5),
                             rep("SSVR.add", 5),
                             rep("SSVR.RBF", 5),
                             rep("hybrid.lin", 5),
                             rep("hybrid.add", 5),
                             rep("hybrid.RBF", 5),
                             rep("evers.lin", 5),
                             rep("evers.RBF", 5),
                             rep("PH", 5), rep("RSF", 5), rep("GBoost", 5))
)
names(bench.lkr.ci) <- c("CIndex", "Model")
bench.lkr.ci$Model <- factor(bench.lkr.ci$Model, levels = c("vanbelle1.lin",
                                                            "vanbelle1.add",
                                                            "vanbelle1.RBF",
                                                            "vanbelle2.lin",
                                                            "vanbelle2.add",
                                                            "vanbelle2.RBF",
                                                            "SSVR.lin",
                                                            "SSVR.add",
                                                            "SSVR.RBF",
                                                            "hybrid.lin",
                                                            "hybrid.add",
                                                            "hybrid.RBF",
                                                            "evers.lin",
                                                            "evers.RBF",
                                                            "PH", "RSF",
                                                            "GBoost"),
                             ordered = TRUE)
ddply(bench.lkr.ci, ~Model, summarise, mean=mean(CIndex),
      sd= round(sd(CIndex), 2))

# ## Leukemia death Benchmark objects #########################################*
bench.lt.ci <- data.frame(c(getBMRPerformances(bench.leuk.death.ssvm.vb1,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ak.ssvm.vb1,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.rbf.ssvm.vb1,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ssvm.vb2,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ak.ssvm.vb2,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.rbf.ssvm.vb2,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ssvm.reg,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ak.ssvm.reg,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.rbf.ssvm.reg,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ssvm.hyb,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ak.ssvm.hyb,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.rbf.ssvm.hyb,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.lkd.survsvm.lin,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.lkd.survsvm.gaus,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.ph,
                                               as.df = TRUE)[, "ci"],
                            getBMRPerformances(bench.leuk.death.rfsrc,
                                               as.df = TRUE)[, "ci"], 
                            getBMRPerformances(bench.leuk.death.gboost,
                                               as.df = TRUE)[, "ci"]),
                          c(rep("vanbelle1.lin", 5),
                            rep("vanbelle1.add", 5),
                            rep("vanbelle1.RBF", 5),
                            rep("vanbelle2.lin", 5),
                            rep("vanbelle2.add", 5),
                            rep("vanbelle2.RBF", 5),
                            rep("SSVR.lin", 5),
                            rep("SSVR.add", 5),
                            rep("SSVR.RBF", 5),
                            rep("hybrid.lin", 5),
                            rep("hybrid.add", 5),
                            rep("hybrid.RBF", 5),
                            rep("evers.lin", 5),
                            rep("evers.RBF", 5),
                            rep("PH", 5), rep("RSF", 5), rep("GBoost", 5))
)
names(bench.lt.ci) <- c("CIndex", "Model")
bench.lt.ci$Model <- factor(bench.lt.ci$Model, levels = c("vanbelle1.lin",
                                                          "vanbelle1.add",
                                                          "vanbelle1.RBF",
                                                          "vanbelle2.lin",
                                                          "vanbelle2.add",
                                                          "vanbelle2.RBF",
                                                          "SSVR.lin",
                                                          "SSVR.add",
                                                          "SSVR.RBF",
                                                          "hybrid.lin",
                                                          "hybrid.add",
                                                          "hybrid.RBF",
                                                          "evers.lin",
                                                          "evers.RBF",
                                                          "PH", "RSF",
                                                          "GBoost"),
                            ordered = TRUE)
ddply(bench.lt.ci, ~Model, summarise, mean=mean(CIndex),
      sd= round(sd(CIndex), 2))


# ## GCSGB2 Benchmark' results ################################################*
bench.gbcsg.ci <- data.frame(c(getBMRPerformances(bench.gbcsg.ssvm.vb1,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.ak.vb1,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.rbf.vb1,
                                                  as.df = TRUE)[,"ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.vb2,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.ak.vb2,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.rbf.vb2,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.reg,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.bgcsg.ssvm.ak.reg,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.bgcsg.ssvm.rbf.reg,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.hyb,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.ssvm.ak.hyb,
                                                  as.df = TRUE)[, "ci"],
                               #getBMRPerformances(bench.gbcsg.ssvm.rbf.hyb,
                               # as.df = TRUE)[,"ci"],
                               # we interrupted execution after 10 days
                               rep(NA, 5),
                               getBMRPerformances(bench.gbcsg.survsvm.lin,
                                                  as.df = TRUE)[, "ci"],
                               #getBMRPerformances(bench.gbcsg.survsvm.gaus,
                               # as.df = TRUE)[,"ci"],
                               #we interrupted execution after 10 days
                               rep(NA, 5),
                               getBMRPerformances(bench.gbcsg.ph,
                                                  as.df = TRUE)[, "ci"],
                               getBMRPerformances(bench.gbcsg.rfsrc,
                                                  as.df = TRUE)[, "ci"], 
                               getBMRPerformances(bench.gbcsg.gboost,
                                                  as.df = TRUE)[, "ci"]),
                             c(rep("vanbelle1.lin", 5),
                               rep("vanbelle1.add", 5),
                               rep("vanbelle1.RBF", 5),
                               rep("vanbelle2.lin", 5),
                               rep("vanbelle2.add", 5),
                               rep("vanbelle2.RBF", 5),
                               rep("SSVR.lin", 5),
                               rep("SSVR.add", 5),
                               rep("SSVR.RBF", 5),
                               rep("hybrid.lin", 5),
                               rep("hybrid.add", 5),
                               rep("hybrid.RBF", 5),
                               rep("evers.lin", 5),
                               rep("evers.RBF", 5),
                               rep("PH", 5), rep("RSF", 5), rep("GBoost", 5))
)
names(bench.gbcsg.ci) <- c("CIndex", "Model")
bench.gbcsg.ci$Model <- factor(bench.gbcsg.ci$Model, levels = c("vanbelle1.lin",
                                                                "vanbelle1.add",
                                                                "vanbelle1.RBF",
                                                                "vanbelle2.lin",
                                                                "vanbelle2.add",
                                                                "vanbelle2.RBF",
                                                                "SSVR.lin",
                                                                "SSVR.add",
                                                                "SSVR.RBF",
                                                                "hybrid.lin",
                                                                "hybrid.add",
                                                                "hybrid.RBF",
                                                                "evers.lin",
                                                                "evers.RBF",
                                                                "PH", "RSF",
                                                                "GBoost"),
                               ordered = TRUE)
ddply(bench.gbcsg.ci, ~Model, summarise, mean=round(mean(CIndex), 2),
      sd= round(sd(CIndex), 2))

bench.mlc.ci <- data.frame(c(getBMRPerformances(bench.mlc.ssvm.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ak.ssvm.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.rbf.ssvm.vb1,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ssvm.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ak.ssvm.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.rbf.ssvm.vb2,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ssvm.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ak.ssvm.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.rbf.ssvm.reg,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ssvm.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ak.ssvm.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.rbf.ssvm.hyb,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.survsvm.lin,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.survsvm.gaus,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.ph,
                                                as.df = TRUE)[, "ci"],
                             getBMRPerformances(bench.mlc.rfsrc,
                                                as.df = TRUE)[, "ci"], 
                             getBMRPerformances(bench.mlc.gboost,
                                                as.df = TRUE)[, "ci"]),
                           c(rep("vanbelle1.lin", 5),
                             rep("vanbelle1.add", 5),
                             rep("vanbelle1.RBF", 5),
                             rep("vanbelle2.lin", 5),
                             rep("vanbelle2.add", 5),
                             rep("vanbelle2.RBF", 5),
                             rep("SSVR.lin", 5),
                             rep("SSVR.add", 5),
                             rep("SSVR.RBF", 5),
                             rep("hybrid.lin", 5),
                             rep("hybrid.add", 5),
                             rep("hybrid.RBF", 5),
                             rep("evers.lin", 5),
                             rep("evers.RBF", 5),
                             rep("PH", 5), rep("RSF", 5), rep("GBoost", 5))
)
names(bench.mlc.ci) <- c("CIndex", "Model")
bench.mlc.ci$Model <- factor(bench.mlc.ci$Model, levels = c("vanbelle1.lin",
                                                            "vanbelle1.add",
                                                            "vanbelle1.RBF",
                                                            "vanbelle2.lin",
                                                            "vanbelle2.add",
                                                            "vanbelle2.RBF",
                                                            "SSVR.lin",
                                                            "SSVR.add",
                                                            "SSVR.RBF",
                                                            "hybrid.lin",
                                                            "hybrid.add",
                                                            "hybrid.RBF",
                                                            "evers.lin",
                                                            "evers.RBF",
                                                            "PH", "RSF",
                                                            "GBoost"),
                             ordered = TRUE)
ddply(bench.mlc.ci, ~Model, summarise, mean=round(mean(CIndex), 2),
      sd= round(sd(CIndex), 2))

all.ci.result <- rbind(bench.vet.ci, bench.lkr.ci, bench.lt.ci,
                       bench.gbcsg.ci, bench.mlc.ci)
datasets <- factor(c(rep("veteran", 85), rep("leuk\\_cr", 85),
                     rep("leuk\\_death", 85), rep("GBSG2", 85), rep("MCLC", 85)), 
                   ordered = TRUE, levels = c("veteran", "leuk\\_cr",
                                              "leuk\\_death", "GBSG2", "MCLC"))
Approach <- NULL
for(i in 1:nrow(all.ci.result)){
  Approach[i] <- if(all.ci.result[i, 2] %in% c("vanbelle1.lin",
                                               "vanbelle1.add",
                                               "vanbelle2.lin",
                                               "vanbelle2.add",
                                               "vanbelle1.RBF",
                                               "vanbelle2.RBF")){
    "ranking"} else {
      if(all.ci.result[i, 2] %in% c("SSVR.lin", "SSVR.add",
                                    "SSVR.RBF")){
        "regression"
      } else {
        if(all.ci.result[i, 2] %in% c("hybrid.lin", "hybrid.add",
                                      "hybrid.RBF")){
          "hybrid"
        } else {
          "reference"
        }
      }
    }
}
Approach <- as.factor(Approach)
Approach <- factor(Approach, levels(Approach)[c(2, 4, 1, 3)])
all.ci.result <- data.frame(all.ci.result, datasets, Approach)
all.ci.result$Model <- factor(all.ci.result$Model, levels = c("vanbelle1.lin",
                                                              "vanbelle1.add",
                                                              "vanbelle1.RBF",
                                                              "vanbelle2.lin",
                                                              "vanbelle2.add",
                                                              "vanbelle2.RBF",
                                                              "SSVR.lin",
                                                              "SSVR.add",
                                                              "SSVR.RBF",
                                                              "hybrid.lin",
                                                              "hybrid.add",
                                                              "hybrid.RBF",
                                                              "evers.lin",
                                                              "evers.RBF",
                                                              "PH", "RSF", 
                                                              "GBoost"),
                              ordered = TRUE)

# ------------------------------------------------------------------------------
# Generate the box plots in a pdf file with a latex form.
# ------------------------------------------------------------------------------
tikz('boxplotci.tex',
     standAlone = TRUE,
     width = 4.5, height = 7,
     packages=c(options()$tikzLatexPackages, "\\usepackage{amsfonts}",
                "\\usepackage[T1]{fontenc}", "\\usepackage{times}"))
p <- ggplot(data = all.ci.result, mapping = aes(y = CIndex, x = Model,
                                                fill = Approach)) + 
  geom_boxplot(outlier.size = .25) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.title=element_text("Approach"),
        panel.background = element_rect(colour = "red")) + 
  labs(x = "Model", y = "C-index")
p + facet_wrap(~datasets, nrow = 5)
dev.off()