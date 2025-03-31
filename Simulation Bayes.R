library(future)
options(future.globals.maxSize = +Inf)
library(DeclareDesign)
library(pwr)
library(rstanarm)
library(broom.mixed)
library(tidyverse)
plan(multisession, workers = 7)

N = 200
g = 2

#ifelse(g==1, 0, ifelse(g==2, 0.25*di, ifelse(g==3,di,2.5*di)))

bay_mod <-
  declare_model(N = N,
                g = g,
                di = pwr::pwr.t.test(N/2,NULL,.01,.8)$d*2,
                dz = ifelse(g==1, 0, ifelse(g==2, 0.25*di, ifelse(g==3,di,2.5*di))),
                X = rbinom(N, 1, 0.5) - 0.5,
                fabricatr::potential_outcomes(Y ~ rnorm(N, dz * Z + di * Z * X))) +
  declare_inquiry(
    ATEZ = mean(Y_Z_1 - Y_Z_0),
    diff_in_CATEs = mean(Y_Z_1[X == 0.5] - Y_Z_0[X == 0.5]) - mean(Y_Z_1[X == -0.5] - Y_Z_0[X == -0.5])
  ) +
  declare_assignment(Z = randomizr::complete_ra(N, prob = 0.5)) +
  declare_measurement(Y = fabricatr::reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z + X + Z:X, term = c("Z","Z:X"), inquiry = c("ATEZ","diff_in_CATEs"), label = "OLS") +
  declare_estimator(Y ~ Z + X + Z:X, term = c("Z","Z:X"), .method = rstanarm::stan_glm, family = gaussian(), refresh = 0, .summary = ~broom.mixed::tidy(.), inquiry = c("ATEZ","diff_in_CATEs"), label = "Bayes")
bay_mod

bay_test = bay_mod |> redesign(N = c(100,200,400,800,1600,3200), g = c(1, 2, 3, 4)) |>
  diagnose_designs(sims=1000, bootstrap_sims = FALSE, diagnosands = declare_diagnosands(power = mean(p.value <= 0.01), mse = mean((estimate - estimand) ^ 2),
                                                                                         rmse = sqrt(mean((estimate - estimand) ^ 2)), coverage = mean(estimand <= conf.high & estimand >= conf.low)))

bay_test$diagnosands_df %>% filter(inquiry=="diff_in_CATEs", g == 3)

## Archiving
save(bay_test,bay_mod,file = "./Bayes.RData",compress = "xz")
rm(bay_test,bay_mod)
load(file = "./Bayes.RData")