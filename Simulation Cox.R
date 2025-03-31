library(future)
options(future.globals.maxSize = +Inf)
library(DeclareDesign)
library(simsurv)
library(rms)
library(survival)
library(tidyverse)
plan(multisession)

#ES based on powerSurvEpi package
cox_sim = function(N) {
  participant_data <- fabricatr::fabricate(
    N = N,
    G = fabricatr::draw_binary(prob = 0.5, N = N) - 0.5,
    X = fabricatr::draw_binary(prob = 0.5, N = N),
    I = ifelse(X==1&G==0.5,1,0)
  )
  p = 0.5
  es = log(uniroot(function(RR) (RR - 1)/(RR + 1) - (qnorm(0.8) + qnorm(0.995))/sqrt(N/2 * p + N/2 * p),c(0,10^5))$root)
  survival_data <- simsurv::simsurv(lambdas = 0.31, gammas = 0.5, x = participant_data, betas = c(I = 2*es), maxt = 5)
  fabricatr::fabricate(survival_data,participant_data) |> dplyr::select(id,eventtime,status,X,G)
}

N = 1000

declaration_18.1 <-
  declare_model(N = N, handler = cox_sim) +
  declare_inquiry(
    ATEXHR = log(weighted.mean(status[X==1],eventtime[X==1])/weighted.mean(status[X==0],eventtime[X==0])),
    diff_in_CATEHR = log((weighted.mean(status[X==1&G==0.5],eventtime[X==1&G==0.5])/weighted.mean(status[X==1&G==-0.5],eventtime[X==1&G==-0.5]))/
                         (weighted.mean(status[X==0&G==0.5],eventtime[X==0&G==0.5])/weighted.mean(status[X==0&G==-0.5],eventtime[X==0&G==-0.5])))
  ) +
  # declare_measurement(Y = fabricatr::reveal_outcomes(Y ~ Z)) +
  declare_estimator(survival::Surv(eventtime,status) ~ X * G, .method = survival::coxph, .summary=broom::tidy, term = c("X","X:G"), inquiry = c("ATEXHR","diff_in_CATEHR"))
declaration_18.1

diagnosis_18.1 <- declaration_18.1 |> redesign(N = c(100,200,400,800,1600,3200)) |> diagnose_design(sims = 10000, diagnosands = declare_diagnosands(
  power = mean(p.value <= 0.01), mse = mean((estimate - estimand) ^ 2), rmse = sqrt(mean((estimate - estimand) ^ 2)), coverage = mean(estimand <= conf.high & estimand >= conf.low)))

diagnosis_18.1 |> get_simulations() |> group_by(inquiry,N) |> summarise(power = mean(p.value<0.01)) |> print(n = 1000)

cox_tab1 = diagnosis_18.1 |> get_simulations() |> filter(inquiry=="diff_in_CATEHR") |> group_by(N) |> summarise(d = exp(mean(estimate)))

# Power graph
i_list = c(75,100,150,200,300,400,600,800,1200,1600,2000)
d_list = 0.5
#for(l in 1:20) {i_list[l]= 1.43^l+8.75}
reps=10000
T=as.data.frame(cbind(n=rep(NA,length(i_list)),d=rep(NA,length(i_list))))
for(i in 1:length(i_list)) {
  T$n[i] = i_list[i]
  T$d[i] = log(uniroot(function(RR) (RR - 1)/(RR + 1) - (qnorm(0.8) + qnorm(0.995))/sqrt(i_list[i] * 0.5 + i_list[i] * 0.5),c(0,10^5))$root)
}
T$n = T$n*2
T$di = exp(T$d*2)
T$d = exp(T$d)
ggplot() +
  geom_line(data=T,aes(n,d,color="Main"),size=1.25) + 
  geom_line(data=T,aes(n,di,color="Int."),size=1.25) +
  geom_point(data=cox_tab1[cox_tab1$N>150,],aes(N,d, shape="Simulation estimates"),size=2) +
  scale_x_continuous(breaks = c(200,500,1000,2000,3000), trans = scales::log10_trans()) +
  scale_shape_discrete("") +
  annotation_logticks(sides = "b") +
  labs(colour = "Type of effect") +
  xlab("Sample size (log scale)") + ylab("Effect size (Hazard ratio)")
ggsave(filename = "Cox.png",width = 6, height = 4, device='png', dpi=1000)

#ES based on powerSurvEpi package (pE = 1 - (1-pC)^HR)
cox_sim_2 = function(N,r,r2) {
  participant_data <- fabricatr::fabricate(
    N = N,
    G = fabricatr::draw_binary(prob = r, N = N) - 0.5,
    X = fabricatr::draw_binary(prob = r2, N = N),
    I = ifelse(X==1&G==0.5,1,0)
  )
  n1 = N * r
  k1 = n1/(N - n1)
  n2 = N * r2
  k2 = n2/(N - n2)
  p0 = 0.5
  es0 = log(uniroot(function(RR) sum(powerSurvEpi::ssizeCT.default(0.8,1,0.5,p0,RR,0.01)) - N,c(0.1,1))$root)
  es1 = log(uniroot(function(RR) sum(powerSurvEpi::ssizeCT.default(0.8,k1,0.5,p0,RR,0.01)) - N,c(0.1,1))$root)
  es2 = log(uniroot(function(RR) sum(powerSurvEpi::ssizeCT.default(0.8,k2,0.5,p0,RR,0.01)) - N,c(0.1,1))$root)
  es0_2 = log(uniroot(function(RR) -(RR - 1)/(RR + 1) - (qnorm(0.8) + qnorm(0.995))/sqrt((N/2 * 0.5 + N/2 * p0)),c(0,1))$root)
  es1_2 = log(uniroot(function(RR) -(RR - 1)/(k1 * RR + 1) - (qnorm(0.8) + qnorm(0.995))/sqrt(k1 * ((n1 * 0.5) + (N-n1) * p0)),c(0,1))$root)
  es2_2 = log(uniroot(function(RR) -(RR - 1)/(k2 * RR + 1) - (qnorm(0.8) + qnorm(0.995))/sqrt(k2 * ((n2 * 0.5) + (N-n2) * p0)),c(0,1))$root)
  esi = ifelse(r==0.5&r2==0.5,es0,ifelse(r==0.5,es2,ifelse(r2==0.5,es1,es1*es2/es0)))
  esi_2 = ifelse(r==0.5&r2==0.5,es0_2,ifelse(r==0.5,es2_2,ifelse(r2==0.5,es1_2,es1_2*es2_2/es0_2)))
  c(exp(esi),exp(esi*2),exp(esi_2),exp(esi_2*2))
  survival_data <- simsurv::simsurv(lambdas = 0.31, gammas = 0.5, x = participant_data, betas = c(I = 2*esi), maxt = 5)
  fabricatr::fabricate(survival_data,participant_data) |> dplyr::select(id,eventtime,status,X,G)
}

N = 4000
r = 0.4
r2 = 0.3

cox_mod_2 <-
  declare_model(N = N, r = r, r2 = r2, handler = cox_sim_2) +
  declare_inquiry(
    ATEXHR = log(1/log(mean(status[X==1]),mean(status[X==0]))),
    diff_in_CATEHR = log(1/log(mean(status[X==1&G==0.5]),mean(status[X==1&G==-0.5])))-log(1/log(mean(status[X==0&G==0.5]),mean(status[X==0&G==-0.5])))
  ) +
  # declare_measurement(Y = fabricatr::reveal_outcomes(Y ~ Z)) +
  declare_estimator(survival::Surv(eventtime,status) ~ X * G, .method = survival::coxph, .summary=broom::tidy, term = c("X","X:G"), inquiry = c("ATEXHR","diff_in_CATEHR"))
cox_mod_2

cox_sims <- cox_mod_2 |> redesign(N = c(4000), r = c(0.1,0.2,0.3,0.4,0.5), r2 = c(0.1,0.2,0.3,0.4,0.5)) |> diagnose_design(sims = 10000, diagnosands = declare_diagnosands(
  power = mean(p.value <= 0.01), mse = mean((estimate - estimand) ^ 2), rmse = sqrt(mean((estimate - estimand) ^ 2)), coverage = mean(estimand <= conf.high & estimand >= conf.low)))

cox_sims |> get_simulations() |> group_by(inquiry,r,r2) |> summarise(power = mean(p.value<0.01), b = sprintf("%0.8f",exp(mean(estimate)))) |> print(n = 1000)

cox_tab2 = cox_sims |> get_simulations() |> filter(inquiry=="diff_in_CATEHR") |> group_by(N) |> summarise(d = exp(mean(estimate)))

## Archiving
save(diagnosis_18.1,declaration_18.1,file = "./Cox.RData",compress = "xz")
rm(diagnosis_18.1,declaration_18.1)
load(file = "./Cox.RData")

save(cox_mod_2,cox_sims,file = "./Cox2.RData",compress = "xz")
rm(cox_mod_2,cox_sims)
load(file = "./Cox2.RData")
