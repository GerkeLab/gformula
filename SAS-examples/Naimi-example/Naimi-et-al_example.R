# Set up and generate data ------------------------------------------------

library(tidyverse)

setwd(here::here("SAS-examples", "Naimi-example"))

# these data are generated according to
# https://academic.oup.com/ije/article/46/2/756/2760169

raw_data <-
"0 0 0 87.288 209271
0 0 1 112.107 93779
0 1 0 119.654 60657
0 1 1 144.842 136293
1 0 0 105.282 134781
1 0 1 130.184 60789
1 1 0 137.720 93903
1 1 1 162.832 210527"

naimi <- read_delim(raw_data, col_names = c("a0", "z1", "a1", "y", "n"), delim = " ")

naimi <- naimi %>%
  rename(y2 = y) %>%
  mutate(row = row_number()) %>%
  nest(data = c(-n, -row)) %>%
  mutate(data = map2(data, n, function(d, n) map_dfr(1:n, ~d))) %>%
  select(-n, -row) %>%
  unnest(data) %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = matches("^[ay]"), names_to = "time", values_to = "a") %>%
  mutate(time = str_remove(time, "^\\w") %>% as.integer(),
         z = if_else(time == 1, z1, NA_real_),
         y = if_else(time == 2, a, NA_real_),
         a = if_else(time == 2, NA_real_, a)) %>%
  select(-z1) %>%
  mutate_at(vars(a, z), as.integer)

# Estimate parameters -----------------------------------------------------

# z1 is caused by a0
z1_fit <- glm(z_1 ~ a_0,
              family = "binomial",
              data = naimi %>%
                filter(time %in% 0:1) %>%
                select(-y) %>%
                pivot_wider(id_cols = id, names_from = time, values_from = c(a, z)) %>%
                select(-a_1, -z_0)
              )
# from the SAS output in the supplement, this should return
# beta_0 = -0.4309 and beta_1 = 0.8735
coef(z1_fit)

## pull functions from example1.R
# naimi %>%
#   nest(data = -time) %>%
#   filter(time == 0) %>%
#   pmap(fit_step, outcome = "z", predictors = "a", id = "id")

# a1 is caused by z1 (and *not* a0; we read this from the DAG. This
# is an important option to consider: we don't always use all previously
# measured variables to predict future variables)
a1_fit <- glm(a_1 ~ z_1,
              family = "binomial",
              data = naimi %>%
                filter(time == 1) %>%
                select(-y) %>%
                pivot_wider(id_cols = id, names_from = time, values_from = c(a, z))
              )
# from the SAS output in the supplement, this should return
# beta_0 = -0.8002 and beta_1 = 1.6084
coef(a1_fit)

# y is caused by a1, z1, a1, *and* an a0*a1 interaction
# another important option to include: permit variable interactions
y_fit <- lm(y_2 ~ a_0 + z_1 + a_1 + a_0*a_1,
              data = naimi %>%
                pivot_wider(id_cols = id, names_from = time, values_from = c(a, z, y)) %>%
                select(-a_2, -z_0, -z_2, -y_0, -y_1)
          )
# from supplement:
# 87.2466 + 24.9699*a1 + 32.5502*z1 + 17.9893*a0 + 0.0541*a1*a0
coef(y_fit)

# MC simulation -----------------------------------------------------------

# set up a function to get a vector of 0/1 from a vector of probabilities
sample_bin <- function(p) {
  map_int(p, ~ sample(0:1, 1, prob = c(1 - .x, .x)))
}

# sample 1000000 a0 entries to initialize the population
intial_data <- sample...

# simulate z1 values

# simulate a1 values

# simulate y values

# match solutions for specific regimes from supplement and paper:
# %gform(run=1,a0=a0,a1=a1,z1=z1);*natural course; y not shown
# %gform(run=2,a0=1,a1=1,z1=z1);*fully exposed; y=150.0
# %gform(run=3,a0=1,a1=0,z1=z1);*time zero only; y=125.0
# %gform(run=4,a0=0,a1=1,z1=z1);*time one only; y=125.0
# %gform(run=5,a0=0,a1=0,z1=z1);*fully unexposed; y=100.0


