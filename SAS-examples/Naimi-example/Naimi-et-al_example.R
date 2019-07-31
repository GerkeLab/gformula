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

naimi %>%
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

