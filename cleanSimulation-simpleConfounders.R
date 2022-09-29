
# install all of these libraries if you haven't yet
# first: install.packages("BiocManager")
# then: 
# BiocManager::install(c("tmle", "SuperLearner", "kableExtra", "tidyverse", "ggplot2", "earth", "ranger", "dagitty"), dependencies = TRUE, force = TRUE)
#
library(tmle) 
library(SuperLearner)
library(kableExtra)
library(tidyverse)
library(ggplot2)
library(earth) 
library(ranger)
library(dagitty)

set.seed(7) # so results are reproducible

generate_data <- function(n, adjSetType) {
  #n = 1000BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
  W1 <- rbinom(n, size = 1, prob=0.2) # binary confounder
  W2 <- rbinom(n, size = 1, prob=0.5) #binary confounder
  W3 <- rbinom(n, size = 1, prob=0.3)
  W4 <- rbinom(n, size = 1, prob=0.2)
  A <- rbinom(n, size=1, prob = plogis(-2+0.2*W1+0.1*W2+0.1*W3+0.1*W4)) # binary treatment
  P1 <- rbinom(n, size=1, prob=0.4)
  P2 <- rbinom(n, size=1, prob=0.3)
  P3 <- rbinom(n, size=1, prob=0.1)
#  Y <- rbinom(n, size =1, prob=plogis(-1+A-0.1*W1+0.2*W2+0.1*W3+0.1*W4+0.1*P1+0.3*P2+0.1*P3))
  # counterfactual
  Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  # Observed outcome
  Y <- Y.1*A + Y.0*(1 - A)
   M1 <- rbinom(n, size=1, prob=plogis(1+0.4*A))
  C1 <- rbinom(n, size=1, prob=plogis(-1+0.3*A+0.4*Y))
  M2 <- rbinom(n, size = 1, prob = plogis(-1+0.8*A))
  WM <- rbinom(n, size=1, prob=plogis(1+0.4*W3+0.7*M2))
  C2 <- rbinom(n, size=1, prob=plogis(-1+0.3*A+0.4*Y))
  M3 <- rbinom(n, size = 1, prob = plogis(-1+0.3*A))
  Chimera <- rbinom(n, size=1, prob=plogis(2+0.3*W4+0.1*M3+.4*C2))
  if (adjSetType == "all") {
    dat <- data.frame(W1, W2, M1, C1, W3, M2, WM, C2, W4, M3, P1, P2, P3, Chimera, A, Y, Y.1, Y.0)
  } else if (adjSetType == "allPure") {
    dat <- data.frame(W1, W2, M1, C1, WM, P1, P2, P3, Chimera, A, Y, Y.1, Y.0)
  } else if (adjSetType == "onlyConfounders") { dat <- data.frame(W1, W2, W3, W4, A, Y, Y.1, Y.0) 
  } else if (adjSetType == "onlyConfoundersAndPrecisionVars") {
    dat <- data.frame(W1, W2, P1, P2, P3, A, Y, Y.1, Y.0)
  } else if (adjSetType == "Confounders_C1_WM_Chimera") {
    dat <- data.frame(W1, W2, W3, W4, C1, WM, Chimera, A, Y, Y.1, Y.0)
  } else if (adjSetType == "onlyBadCovariates") {
    dat <- data.frame(C1, M1, M2, C2, WM, Chimera, A, Y, Y.1, Y.0)
  }
  return(dat)
}

generateData<- function(n){
  w1 <- rbinom(n, size=1, prob=0.5)
  w2 <- rbinom(n, size=1, prob=0.65)
  w3 <- round(runif(n, min=0, max=4), digits=0)
  w4 <- round(runif(n, min=0, max=5), digits=0)
  A <- rbinom(n, size=1, prob= plogis(-5 + 0.05*w2 + 0.25*w3 + 0.6*w4 + 0.4*w2*w4))
  # the counterfactuals aka potential outcomes
  Y.1 <- rbinom(n, size=1, prob= plogis(-1 + 1 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  Y.0 <- rbinom(n, size=1, prob= plogis(-1 + 0 -0.1*w1 + 0.35*w2 + 0.25*w3 + 0.20*w4 + 0.15*w2*w4))
  # Observed outcome
  Y <- Y.1*A + Y.0*(1 - A)
  # return data.frame
  data.frame(w1, w2, w3, w4, A, Y, Y.1, Y.0)
}
#quantifyEffectWithSubset()

n = 20000
#dat_obs <- generate_data(n, "onlyBadCovariates")
dat_obs <- generate_data(n, "onlyConfounders")
#dat_obs <- generate_data(n, "onlyConfoundersAndPrecisionVars")
#dat_obs <- generate_data(n, "all")
#dat_obs <- generate_data(n, "allPure")

#dat <- data.frame(dat_obs)
#
#kable(head(dat_obs), digits=2, caption = "Simulated dataset.")

# So we will start out with a simulation containing only
# four "well-behaved" / "pure" confounders (W1, W2, W3, and W4)
# plus the exposure (A) and outcome variables (Y)
#############

# Now starting with where the original tutorial left off
#  with a more solid data generating process

# Source: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fsim.7628&file=sim7628-sup-0003-Appendix.pdf

# the gold standard: 

# True ATE in the population
set.seed(7777)


n = 20000
#dat_obs <- generate_data(n, "onlyBadCovariates")
#dat_obs <- generate_data(n, "onlyConfounders")

ObsDataTrueATE <- generate_data(n = 500000, "onlyConfounders")
True_EY.1 <- mean(ObsDataTrueATE$Y.1)
True_EY.0 <- mean(ObsDataTrueATE$Y.0)
True_ATE <- True_EY.1-True_EY.0 ;True_ATE
True_MOR <- (True_EY.1*(1-True_EY.0))/((1-True_EY.1)*True_EY.0);True_MOR
cat("\n True_ATE:", abs(True_ATE))
#

