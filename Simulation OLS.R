library(future)
options(future.globals.maxSize = +Inf)
library(DeclareDesign)
library(pwr)
library(tidyverse)
plan(multisession, workers = 7)

N = 200
g = 2

#ifelse(g==1, 0, ifelse(g==2, 0.25*di, ifelse(g==3,di,2.5*di)))

ols_mod <-
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
  declare_estimator(Y ~ Z + X + Z:X, term = c("Z","Z:X"), inquiry = c("ATEZ","diff_in_CATEs"), label = "OLS")
ols_mod

ols_test = ols_mod |> redesign(N = c(100,200,400,800,1600,3200), g = c(1, 2, 3, 4)) |>
  diagnose_designs(sims=10000, bootstrap_sims = FALSE, diagnosands = declare_diagnosands(power = mean(p.value <= 0.01), mse = mean((estimate - estimand) ^ 2),
                                                                                      rmse = sqrt(mean((estimate - estimand) ^ 2)), coverage = mean(estimand <= conf.high & estimand >= conf.low)))

ols_test$diagnosands_df %>% filter(inquiry=="diff_in_CATEs", g == 3)

# N = 200
# di = pwr::pwr.t.test(N/2,NULL,.01,.8)$d*2
# dz = di
# 
# #ifelse(g==1, 0, ifelse(g==2, 0.25*di, ifelse(g==3,di,2.5*di)))
# 
# ols_mod <-
#   declare_model(N = N,
#                 di = di,
#                 dz = dz,
#                 X = rbinom(N, 1, 0.5) - 0.5,
#                 fabricatr::potential_outcomes(Y ~ rnorm(N, dz * Z + di * Z * X))) +
#   declare_inquiry(
#     ATEZ = mean(Y_Z_1 - Y_Z_0),
#     diff_in_CATEs = mean(Y_Z_1[X == 0.5] - Y_Z_0[X == 0.5]) - mean(Y_Z_1[X == -0.5] - Y_Z_0[X == -0.5])
#   ) +
#   declare_assignment(Z = randomizr::complete_ra(N, prob = 0.5)) +
#   declare_measurement(Y = fabricatr::reveal_outcomes(Y ~ Z)) +
#   declare_estimator(Y ~ Z + X + Z:X, term = c("Z","Z:X"), inquiry = c("ATEZ","diff_in_CATEs"), label = "OLS")
# ols_mod
# 
# # test = list(N = N2[1], di = di2[1], d = d2[1])
# # N2 = c(100,200,400,800,1600,3200)
# # di2 = map(N2,~pwr::pwr.t.test(.x/2,NULL,.01,.8)$d*2)
# # d2 = di2
# 
# ols_test = ols_mod |> redesign(N = rep(c(100,200,400,800,1600,3200),4), 
#                                di = rep(c(1.3905896, 0.9747619, 0.6863514, 0.4842827, 0.3420996, 0.2417603),4),
#                                dz = c(0, 0, 0, 0, 0, 0, 0.3476474, 0.243690475, 0.17158785, 0.121070675, 0.0855249, 0.060440075, 1.3905896, 0.9747619, 0.6863514, 0.4842827, 0.3420996, 0.2417603, 3.476474, 2.43690475, 1.7158785, 1.21070675, 0.855249, 0.60440075), 
#                                expand = FALSE) |> 
#   diagnose_designs(sims=1000, bootstrap_sims = FALSE, diagnosands = declare_diagnosands(power = mean(p.value <= 0.01), mse = mean((estimate - estimand) ^ 2), 
#                                                                                         rmse = sqrt(mean((estimate - estimand) ^ 2)), coverage = mean(estimand <= conf.high & estimand >= conf.low)))

# test |> get_simulations() |> group_by(estimator,base,inquiry,N) |> summarise(power = mean(p.value<0.01)) |> print(n = 1000)

## Archiving
save(ols_test,ols_mod,file = "./OLS.RData",compress = "xz")
rm(ols_test,ols_mod)
load(file = "./OLS.RData")
