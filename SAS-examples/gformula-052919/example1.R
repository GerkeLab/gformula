# Set up and load ---------------------------------------------------------

library(tidyverse)

setwd(here::here("SAS-examples", "gformula-052919"))

example1 <- read_csv("example1.csv")

# SAS macro call ----------------------------------------------------------
# output for this call is showwn on page 8-9 of gformula3_documentation.pdf

# %gformula(
#    data= sample,
#    id=id,
#    time=time,
#    timepoints = 6,
#    outc=dia,
#    outctype=binsurv,
#    comprisk =  dead  ,
#
#    fixedcov = baseage,
#    timeptype= concat,
#    timeknots = 1 2 3 4 5,
#
#    ncov=2,
#    cov1  = hbp,    cov1otype  = 2, cov1ptype = tsswitch1,
#    cov2  = act,    cov2otype  = 4, cov2ptype = lag2cub,
#
#    hazardratio = 1 ,
#    intcomp = 0 1 ,
#    seed= 9458, nsamples = 0, numint=1
# );


# R call ------------------------------------------------------------------

# time 1: we will use variables measured at time 0 to predict the levels of
# variables measured at time 1
t1 <- example1 %>% filter(time %in% c(0,1))

# verify that baseage remains constant across t0 and t1
t1 %>% distinct(id, baseage) %>% count(id, baseage, sort = TRUE) %>% count(n)

# make a flat file to permit use of glm() and similar models
t1_flat <- t1 %>% tidyr::pivot_wider(
   id_cols = c(id, baseage),
   names_from = time,
   values_from = c(-id, -baseage))

# note that the user's ordering of cov1, cov2, etc impacts the covariate
# adjustment/Markov decomposition ordering!
# details in Step 1a of algorithm outline in documentation.pdf

# fit model for hbp_1
# note that this is covXotype = 2:
# Fits model only to records where the first lagged value of covX=0.
# Should be used for binary covariates that, once they switch from 0 to 1,
# they stay 1 (e.g. indicator of diabetes diagnosis)
hbp1_fit <- glm(hbp_1 ~ baseage,
                family = "binomial",
                data = t1_flat %>% filter(hbp_0 == 0))

# fit model for act_1
# note that this is covXotype = 4: see p 10 of documentation.pdf
# In short, it's a two-stage model, where the first estimates
# whether act_1==0. If not, it's a log-linear model.
act1_fit_z <- glm(act_1 > 0 ~ baseage + hbp_0,
                  family = "binomial",
                  data = t1_flat)
act1_fit_I <- lm(log(act_1) ~ baseage + hbp_0,
                  data = t1_flat %>% filter(act_1 > 0))




