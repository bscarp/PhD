library(future)
options(future.globals.maxSize = +Inf)
library(DeclareDesign)
library(tidyverse)
plan(multisession)

ES.h.inv = function(np1,h) {
  (sin((2 * asin(sqrt(np1)) - h)/2))^2
  } 

N = 400
base = 0.5
# dire = -1

lpm_mod <-
  declare_model(N = N,
                base = base,
                d = ES.h.inv(base,h = -(pwr::pwr.2p.test(NULL,N/2,0.01,0.8)$h)) - base,
                di = ifelse((base+2*d)>1,1-d-base,2*d),
                X = rbinom(N, 1, 0.5) - 0.5,
                fabricatr::potential_outcomes(Y ~ rbinom(N, 1, prob = base + d * Z + di * Z * X))) +
  declare_inquiry(
    ATEZ = mean(Y_Z_1 - Y_Z_0),
    diff_in_CATEs = mean(Y_Z_1[X == 0.5] - Y_Z_0[X == 0.5]) - mean(Y_Z_1[X == -0.5] - Y_Z_0[X == -0.5]),
    ATEZOR = log(mean(Y_Z_1)/mean(Y_Z_0)),
    diff_in_CATEOR = log((mean(Y_Z_1[X == 0.5])/mean(Y_Z_0[X == 0.5]))/(mean(Y_Z_1[X == -0.5])/mean(Y_Z_0[X == -0.5])))
  ) +
  declare_assignment(Z = randomizr::complete_ra(N, prob = 0.5)) +
  declare_measurement(Y = fabricatr::reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z + X + Z:X, term = c("Z","Z:X"), inquiry = c("ATEZ","diff_in_CATEs"), label = "OLS") +
  declare_estimator(Y ~ Z + X + Z:X, .method = glm, family = binomial("logit"), term = c("Z","Z:X"), inquiry = c("ATEZOR","diff_in_CATEOR"), label = "Logit")
lpm_mod

test = lpm_mod |> redesign(N = c(100,200,400,800,1600,3200), base = c(0.1,0.3,0.5)) |> diagnose_designs(sims=10000, diagnosands = declare_diagnosands(
  power = mean(p.value <= 0.01), mse = mean((estimate - estimand) ^ 2), rmse = sqrt(mean((estimate - estimand) ^ 2)), coverage = mean(estimand <= conf.high & estimand >= conf.low)))

test |> get_simulations() |> group_by(estimator,base,inquiry,N) |> summarise(power = mean(p.value<0.01)) |> print(n = 1000)

lpm_tab1 = test |> get_simulations() |> filter(estimator=="OLS",inquiry=="diff_in_CATEs"&base==0.3) |> group_by(N) |> summarise(d = mean(estimate), cov = mean(coverage), rmse = mean(rmse))

# Power graph
i_list = c(75,100,150,200,300,400,600,800,1200,1600,2000)
d_list = 0.5
reps=10000
T=as.data.frame(cbind(n=rep(NA,length(i_list)),d=rep(NA,length(i_list))))
for(i in 1:length(i_list)) {
  T$n[i] = i_list[i]
  T$d[i] = ES.h.inv(0.3,h = -pwr::pwr.2p.test(NULL,i_list[i],0.01,0.8)$h) - 0.3
}
T$n = T$n*2
T$di = T$d*2
ggplot() +
  geom_line(data=T,aes(n,d,color="Main"),size=1.25) + 
  geom_line(data=T,aes(n,di,color="Int."),size=1.25) +
  geom_point(data=lpm_tab1[lpm_tab1$N>150,],aes(N,d, shape="Simulation estimates"),size=2) +
  scale_x_continuous(breaks = c(200,500,1000,2000,3000), trans = scales::log10_trans()) +
  scale_shape_discrete("") +
  annotation_logticks(sides = "b") +
  labs(colour = "Type of effect") +
  xlab("Sample size (log scale)") + ylab("Effect size (diffence in proportions from 0.3)")
ggsave(filename = "Binary.png",width = 6, height = 4, device='png', dpi=1000)

N = 2000
base = 0.5
r = 0.4
r2 = 0.3

#Imbalance check
lpm_mod_2 <-
  declare_model(N = N,
                base = base,
                d0 = ES.h.inv(base,h = -(pwr::pwr.2p.test(NULL,N/2,0.01,0.8)$h)) - base,
                d1 = ES.h.inv(base,h = -(pwr::pwr.2p2n.test(NULL,N*r,N*(1-r),0.01,0.8)$h)) - base,
                d2 = ES.h.inv(base,h = -(pwr::pwr.2p2n.test(NULL,N*r2,N*(1-r2),0.01,0.8)$h)) - base,
                d = d1*d2/d0,
                di = ifelse((base+2*d)>1,1-d-base,2*d),
                X = rbinom(N, 1, r) - 0.5,
                fabricatr::potential_outcomes(Y ~ rbinom(N, 1, prob = base + d * Z + di * Z * X))) +
  declare_inquiry(
    ATEZ = mean(Y_Z_1 - Y_Z_0),
    diff_in_CATEs = mean(Y_Z_1[X == 0.5] - Y_Z_0[X == 0.5]) - mean(Y_Z_1[X == -0.5] - Y_Z_0[X == -0.5]),
    ATEZOR = log(mean(Y_Z_1)/mean(Y_Z_0)),
    diff_in_CATEOR = log((mean(Y_Z_1[X == 0.5])/mean(Y_Z_0[X == 0.5]))/(mean(Y_Z_1[X == -0.5])/mean(Y_Z_0[X == -0.5])))
  ) +
  declare_assignment(Z = randomizr::complete_ra(N, prob = r2)) +
  declare_measurement(Y = fabricatr::reveal_outcomes(Y ~ Z)) +
  declare_estimator(Y ~ Z + X + Z:X, term = c("Z","Z:X"), inquiry = c("ATEZ","diff_in_CATEs"), label = "OLS") +
  declare_estimator(Y ~ Z + X + Z:X, .method = glm, family = binomial("logit"), term = c("Z","Z:X"), inquiry = c("ATEZOR","diff_in_CATEOR"), label = "Logit")
lpm_mod_2

lpm_sims = lpm_mod_2 |> redesign(N = c(2000), base = c(0.5), r = c(0.1,0.2,0.3,0.4,0.5), r2 = c(0.1,0.2,0.3,0.4,0.5)) |> diagnose_designs(sims=10000, diagnosands = declare_diagnosands(
  power = mean(p.value <= 0.01), mse = mean((estimate - estimand) ^ 2), rmse = sqrt(mean((estimate - estimand) ^ 2)), coverage = mean(estimand <= conf.high & estimand >= conf.low)))

lpm_sims |> get_simulations() |> group_by(estimator,r,r2,inquiry) |> summarise(power = mean(p.value<0.01), b = sprintf("%0.8f",mean(estimate))) |> print(n = 1000)

lpm_tab2 = lpm_sims |> get_simulations() |> filter(estimator=="OLS",inquiry=="diff_in_CATEs"&r==0.5&r2==0.5) |> summarise(d = mean(estimate))

## Archiving
save(test,lpm_mod,file = "./LPM.RData",compress = "xz")
rm(test,lpm_mod)
load(file = "./LPM.RData")

save(lpm_mod_2,lpm_sims,file = "./LPM2.RData",compress = "xz")
rm(lpm_mod_2,lpm_sims)
load(file = "./LPM2.RData")
