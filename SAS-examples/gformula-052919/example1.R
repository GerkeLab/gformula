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

# ---- Functions ----

model_formula <- function(outcome,predictors) {
  rlang::new_formula(
    rlang::parse_expr(outcome),
    rlang::parse_expr(predictors)
  )
}

append_time_suffix <- function(x,time=0,...) UseMethod("append_time_suffix")

append_time_suffix.default <- function(x,time=0) {
  paste0(x, sprintf("_T%02d", time))
}

append_time_suffix.data.frame <- function(x,time=0,exclude=NULL) {
  idx_names <- which(!names(x) %in% exclude)
  names(x)[idx_names] <- append_time_suffix(names(x)[idx_names], time)
  x
}

select_cols_in_string <- function(data, string) {
  cols_keep <- c()
  for (col in colnames(data)) {
    col_rgx <- paste0("\\b", regex_escape_dot(col), "\\b")
    if (grepl(col_rgx, string)) {
      cols_keep <- c(cols_keep, col)
    }
  }
  data[, intersect(colnames(data), cols_keep)]
}

regex_escape_dot <- function(str) {
  gsub("([.])", "\\\\.", str)
}

fit_step <- function(
  data,
  data_next,
  time,
  ...,
  outcome="hbp",
  predictors="baseage",
  filter_at_step=NULL,
  cov_type=2,
  id=c("id", "baseage")
) {
  # apply filter to data (current step)
  # append time suffix to outcome (+1) and predictors
  # append time suffix to data, data_next then left_join using id
  # construct model formula
  # run model

  if (is.null(data_next) || is.na(data_next) || !nrow(data_next)) {
    return(NULL)
  }

  for (col in names(data)) {
    if (col %in% id) next
    col_rgx <- paste0("\\b", regex_escape_dot(col), "\\b")
    outcome <- sub(col_rgx, append_time_suffix(col, time + 1), outcome)
    predictors <- gsub(col_rgx, append_time_suffix(col, time), predictors)
  }

  filter_at_step <- rlang::enexpr(filter_at_step)
  if (!is.null(filter_at_step)) {
    data <- data %>% filter(!!filter_at_step)
  }

  data <- data %>% append_time_suffix(time, exclude = id)
  data_next <- data_next %>% append_time_suffix(time + 1, exclude = id)

  str_formula <- paste(outcome, predictors)
  data <- left_join(data, data_next, by = id) %>%
    select_cols_in_string(str_formula)

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
example_nested$m_hbp <-
  example_nested %>%
  pmap(fit_step, outcome = "hbp", predictors = "baseage", filter_at_step = hbp == 0, id = example_ids)

#' The `act` covariate is of type 4 (borrowing terminology from `gformula.sas`),
#' and `fit_step()` returns a list of models. The model in `.$zero` estimates
#' whether `act` is non-zero and the model `.$nzero` estimates `act` given that
#' it is non-zero. (We can change this structure/naming as needed.)
example_nested$m_act <-
  example_nested %>%
  pmap(fit_step, outcome = "act", predictors = "baseage + hbp + act", cov_type = 4, id = example_ids)

example_nested$m_dead <-
  example_nested %>%
  pmap(fit_step, outcome = "dead", predictors = "baseage + hbp + act", id = example_ids)

example_nested$m_dia <-
  example_nested %>%
  pmap(fit_step, outcome = "dia", predictors = "baseage + hbp + act", id = example_ids, filter_at_step = dead == 0)

example_nested


# Monte Carlo simulation --------------------------------------------------

# set up a function to get a vector of 0/1 from a vector of probabilities
sample_bin <- function(p) {
  map_int(p, ~ sample(0:1, 1, prob = c(1 - .x, .x)))
}

sample_model <- function(model, data) {
  predict(model, newdata = data, type = "response") %>%
    sample_bin()
}

sample_model_bilevel <- function(model, data) {
  pred <- predict(model$zero, newdata = data, type = "response") %>%
    sample_bin()

  idx_non_zero <- which(pred > 0)

  pred[idx_non_zero] <- predict(model$nzero, newdata = data[idx_non_zero, ]) +
    rnorm(length(idx_non_zero), 0, sd(model$nzero$residuals))

  pred[idx_non_zero] <- exp(pred[idx_non_zero])

  pred
}

sim_next_step <- function(data, model, data_next, var, ..., type = "single") {
  data_next[[var]] <- switch(
    type,
    bilevel = sample_model_bilevel(model, data),
    sample_model(model, data)
  )
  data_next
}

initial_data <-
  example1 %>%
  filter(time == 0) %>%
  { .[complete.cases(.), ] } %>%
  select(-contains("_l"), -censlost) %>%
  sample_n(size = 1e4, replace = TRUE) %>%
  mutate(id = row_number()) %>%
  nest(data = -time)

initial_data_new <-
  initial_data %>%
  mutate(data_next = map(data, ~ select(., id, baseage))) %>%
  select(-data) %>%
  expand(time = 0:5, data_next)

pop <-
  crossing(time = 0:5) %>%
  left_join(example_nested %>% select(-contains('data')), by = "time") %>%
  left_join(initial_data, by = "time") %>%
  left_join(initial_data_new, by = "time")

for (idx_time in seq_along(pop$time)) {
  if (idx_time == nrow(pop)) next

  tstep <- pop$time[idx_time]

  data_this <- pop$data[[idx_time]] %>%
    append_time_suffix(tstep, exclude = c("id", "baseage"))

  # hbp
  pop$data_next[[idx_time]] <- sim_next_step(
    data = data_this,
    model = pop$m_hbp[[idx_time]],
    data_next = pop$data_next[[idx_time]],
    var = "hbp"
  )

  # act
  pop$data_next[[idx_time]] <- sim_next_step(
    data = data_this,
    model = pop$m_act[[idx_time]],
    data_next = pop$data_next[[idx_time]],
    var = "act",
    type = "bilevel"
  )

  # dead
  pop$data_next[[idx_time]] <- sim_next_step(
    data = data_this,
    model = pop$m_dead[[idx_time]],
    data_next = pop$data_next[[idx_time]],
    var = "dead"
  )

  # dia
  pop$data_next[[idx_time]] <- sim_next_step(
    data = data_this,
    model = pop$m_dia[[idx_time]],
    data_next = pop$data_next[[idx_time]],
    var = "dia"
  )

  pop$data_next[[idx_time]][pop$data_next[[idx_time]]$dead == 1L, 'dia'] <- 0L

  # done move data_next to next row of data
  pop$data[[idx_time + 1]] <- pop$data_next[[idx_time]]
}
