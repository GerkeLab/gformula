# Set up and load ---------------------------------------------------------

library(tidyverse)

setwd(here::here("SAS-examples", "gformula-052919"))

example1 <- read_csv("example1.csv")

# output for this call is showwn on page 8-9 of gformula3_documentation.pdf

# SAS macro call ----------------------------------------------------------
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
# Note: somewhere along the way, we need to incorporate the lag option which
# will be variable-specific (e.g. some will have lag==2, others lag==1)
t1 <- example1 %>% filter(time %in% c(0, 1))

# verify that baseage remains constant across t0 and t1
t1 %>%
  distinct(id, baseage) %>%
  count(id, baseage, sort = TRUE) %>%
  count(n)

# make a flat file to permit use of glm() and similar models
t1_flat <- t1 %>% tidyr::pivot_wider(
  id_cols = c(id, baseage),
  names_from = time,
  values_from = c(-id, -baseage)
)

# note that the user's ordering of cov1, cov2, etc impacts the covariate
# adjustment/Markov decomposition ordering!
# details in Step 1a of algorithm outline in documentation.pdf

# fit model for hbp_1
# note that this is covXotype = 2:
# Fits model only to records where the first lagged value of covX=0.
# Should be used for binary covariates that, once they switch from 0 to 1,
# they stay 1 (e.g. indicator of diabetes diagnosis)
hbp1_fit <- glm(
  hbp_1 ~ baseage,
  family = "binomial",
  data = t1_flat %>% filter(hbp_0 == 0)
)

# fit model for act_1
# note that this is covXotype = 4: see p 10 of documentation.pdf
# In short, it's a two-stage model, where the first estimates
# whether act_1==0. If not, it's a log-linear model.
act1_fit_z <- glm(
  act_1 > 0 ~ baseage + hbp_0,
  family = "binomial",
  data = t1_flat
)
act1_fit_I <- lm(
  log(act_1) ~ baseage + hbp_0,
  data = t1_flat %>% filter(act_1 > 0)
)

# fit model for competing risk (death)
dead1_fit <- glm(
  dead_1 ~ baseage + hbp_0 + act_0,
  family = "binomial",
  data = t1_flat
)

# fit model for event of interest (diabetes)
# Note: subset to those without competing event at time 1,
# need to circle back to methods papers to ensure this is correct
# but it seems to be from Keil, 2014 Figure 2
dia1_fit <- glm(
  dia_1 ~ baseage + hbp_0 + act_0,
  family = "binomial",
  data = t1_flat %>% filter(dead_1 == 0)
)


# Garrick -----------------------------------------------------------------


# ---- Functions ----

model_formula <- function(outcome, predictors) {
  rlang::new_formula(
    rlang::parse_expr(outcome),
    rlang::parse_expr(predictors)
  )
}

append_time_suffix <- function(x, time = 0, ...) UseMethod("append_time_suffix")

append_time_suffix.default <- function(x, time = 0) {
  paste0(x, sprintf("_T%02d", time))
}

append_time_suffix.data.frame <- function(x, time = 0, exclude = NULL) {
  idx_names <- which(!names(x) %in% exclude)
  names(x)[idx_names] <- append_time_suffix(names(x)[idx_names], time)
  x
}

fit_step <- function(
  data,
  data_next,
  time,
  ...,
  outcome = "hbp",
  predictors = "baseage",
  filter_at_step = NULL,
  cov_type = 2,
  id = c("id", "baseage")
) {
  # apply filter to data (current step)
  # append time suffix to outcome (+1) and predictors
  # append time suffix to data, data_next then left_join using id
  # construct model formula
  # run model

  if (is.null(data_next) || is.na(data_next) || !nrow(data_next)) return(NULL)

  for (col in names(data)) {
    if (col %in% id) next
    outcome <- sub(col, append_time_suffix(col, time + 1), outcome, fixed = TRUE)
    predictors <- gsub(col, append_time_suffix(col, time), predictors, fixed = TRUE)
  }

  filter_at_step <- rlang::enexpr(filter_at_step)
  if (!is.null(filter_at_step)) {
    data <- data %>% filter(!!filter_at_step)
  }

  data <- data %>% append_time_suffix(time, exclude = id)
  data_next <- data_next %>% append_time_suffix(time + 1, exclude = id)

  data <- left_join(data, data_next, by = id)

  if (cov_type == 2) {
    # note that this is covXotype = 2:
    # Fits model only to records where the first lagged value of covX=0.
    # Should be used for binary covariates that, once they switch from 0 to 1,
    # they stay 1 (e.g. indicator of diabetes diagnosis)

    formula <- model_formula(outcome, predictors)
    model_call <- rlang::expr(glm(!!formula, data = data, family = "binomial"))
    eval(model_call)

  } else if (cov_type == 4) {
    # note that this is covXotype = 4: see p 10 of documentation.pdf
    # In short, it's a two-stage model, where the first estimates
    # whether act_1==0. If not, it's a log-linear model.
    outcome_gt_0 <- paste(outcome, "> 0")

    formulas <- list(
      zero = model_formula(outcome_gt_0, predictors),
      nzero = model_formula(paste0("log(", outcome, ")"), predictors)
    )

    data_nzero <- data %>% filter(!!rlang::parse_expr(outcome_gt_0))

    list(
      zero = rlang::expr(glm(!!formulas$zero, data = data, family = "binomial")),
      nzero = rlang::expr(
        lm(!!formulas$nzero, data = data_nzero)
      )
    ) %>%
      map(~ eval(.))
  }
}

# ---- Example ----

#' This is an example of the expected data structure using nested data frames.
#' The data at each step are nested in to the `data` column, indexed with the
#' `time` column. In this example we only need `data` (the current step) and
#' `data_next` (data at the next time step) but we can similarly use `lag()` to
#' access data at previous time steps as needed.
example_nested <-
  example1 %>%
  select(id, baseage, time, everything()) %>%
  nest(data = -time) %>%
  mutate(
    data = map(data, as_tibble),
    data_next = lead(data, 1)
  )

#' These variables do not vary (and are used to match subjects between time steps)
example_ids <- c("id", "baseage")

#' This example repeats the code above, but for all time steps in the example.
#' Models are returned into columns named after the variables they predict. They
#' are stored in the row corresponding to the model's input data (rather than
#' the output time step).
example_nested$hbp <-
  example_nested %>%
  pmap(fit_step, outcome = "hbp", predictors = "baseage", filter_at_step = hbp == 0, id = example_ids)

#' The `act` covariate is of type 4 (borrowing terminology from `gformula.sas`),
#' and `fit_step()` returns a list of models. The model in `.$zero` estimates
#' whether `act` is non-zero and the model `.$nzero` estimates `act` given that
#' it is non-zero. (We can change this structure/naming as needed.)
example_nested$act <-
  example_nested %>%
  pmap(fit_step, outcome = "act", predictors = "baseage + hbp", cov_type = 4, id = example_ids)

example_nested$dead <-
  example_nested %>%
  pmap(fit_step, outcome = "dead", predictors = "baseage + hbp + act", id = example_ids)

example_nested$dia <-
  example_nested %>%
  pmap(fit_step, outcome = "dia", predictors = "baseage + hbp + act", id = example_ids, filter_at_step = dead == 0)

example_nested
