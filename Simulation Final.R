ifelse(Sys.info()[1]=="Windows",library(doParallel),library(doMC))
if(Sys.info()[1]=="Windows") registerDoParallel(cores=4) else registerDoMC(cores=4)
library(foreach)
library(tidyverse)
# library(flexiblas)
# library(BayesFactor)
# library(BAS)
# library(brms)
library(rms)
library(pwr)
library(ggpubr)
library(car)

n_list = c(60,100,200,400,800,1600,3200)
g_list = c(0:5)
reps=10000
B3 = brm(Y~X*G, df, prior = c(prior(normal(0, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b"), prior(cauchy(0, 1), class = "sigma")),
         sample_prior = TRUE,save_pars = save_pars(all=T))
start_time=Sys.time()
d3 = as.data.frame(foreach(n = n_list, .combine="rbind") %:% foreach(g = g_list, .combine="rbind") %:% foreach(i = 1:reps, .combine="rbind",.packages = c("pwr","BAS","MASS"))  %dopar% {
  # t0 = Sys.time()
  d = pwr.t.test(n/2,NULL,.01,.8)$d*2
  a = ifelse(g==1,.5*d,ifelse(g==2,.5*d,ifelse(g==3,.75*d,ifelse(g==4,d,ifelse(g==5,2.5*d,0)))))
  b = ifelse(g==1,.5*d,ifelse(g==2,-.5*d,ifelse(g==3,-.25*d,ifelse(g==4,0,ifelse(g==5,1.5*d,0)))))
  X = rep(c(0,1,0,1),each=n/4)
  G = rep(c(-0.5,0.5),each=n/2)
  df = as.data.frame(mvrnorm(n,c(0),matrix(c(1),nrow=1)))
  names(df) = c("Y")
  df$Y[(n/4+1):(n/2)] = df$Y[(n/4+1):(n/2)] + b
  df$Y[(3*n/4+1):n] = df$Y[(3*n/4+1):n] + a
  df = cbind(df,X,G)
  # t1 = Sys.time() - t0
  Z1 = lm(Y~X+G+X:G,df)
  Z2 = lm(Y~X,df[df$G==-0.5,])
  Z3 = lm(Y~X,df[df$G==0.5,])
  # Z4 = lm(Y~G,df[df$X==1,])
  # t2 = Sys.time() - t0 - t1
  # B1 = generalTestBF(Y~X+G+X:G,df)
  # B1r = extractBF(B1[4])$bf/extractBF(B1[3])$bf
  # P1 = posterior(B1[4],iterations = 10000)
  # t3 = Sys.time() - t0 - t1 - t2
  B2 = bas.lm(Y~X*G,df,prior="JZS",include.always = ~1+X+G)
  B2r = ifelse(summary(B2)[4,2]==1,summary(B2)[5,2]/summary(B2)[5,3],summary(B2)[5,3])
  B2b = bas.lm(Y~X*G,df,prior="JZS",include.always = ~1+X*G)
  # t4 = Sys.time() - t0 - t1 - t2 - t3
  # B3f = update(B3,newdata=df)
  # B3r = hypothesis(B3f, "X:G = 0")$hypothesis$Evid.Ratio
  # B3p = fixef(B3f)[4,]
  # t5 = Sys.time() - t0 - t1 - t2 - t3 - t4
  # bt = bartlett.test(Y~X,df)$p.value
  # ks = ks.test(df$Y[df$X==1]-mean(df$Y[df$X==1]),df$Y[df$X==0])$p.value
  # c(i,g,n,a,b,t1,t2,t3,t4,t5,B1r,B2r,B3r,anova(Z1)[c(1:3),5],anova(Z2)[1,5],anova(Z3)[1,5],Z1$coefficients[4],confint(Z1)[4,],summary(P1)$statistics[4,1],summary(P1)$quantiles[4,c(1,5)],confint(coef(B2b))[4,c(3,1,2)],B3p[c(1,3,4)])
  c(i,g,n,a,b,B2r,anova(Z1)[c(1:3),5],anova(Z2)[1,5],anova(Z3)[1,5],Z2$coefficients[2],Z3$coefficients[2],Z1$coefficients[2:4],confint(Z1)[4,],confint(coef(B2b))[4,c(3,1,2)])
})
end_time=Sys.time()
(time=end_time-start_time)
# names(d3) = c("i","g","n","a","b","t1","t2","t3","t4","t5","B1r","B2r","B3r","px","pg","pi","p1","p2","bzxg","bzlci","bzuci","bb1xg","bb1lci","bb1uci","bb2xg","bb2lci","bb2uci","bb3xg","bb3lci","bb3uci")
names(d3) = c("i","g","n","a","b","bfi","px","pg","pi","p1","p2","bzx1","bzx2","bzx","bzg","bzi","bzlci","bzuci","bbi","bblci","bbuci")

#Paper 4, Table 1
t1 = data.frame(g=aggregate(bfi>3~a+b+g+n,d3,mean)[,3],
                n=aggregate(bfi>3~a+b+g+n,d3,mean)[,4],
                a=aggregate(bfi>3~a+b+g+n,d3,mean)[,1],
                b=aggregate(bfi>3~a+b+g+n,d3,mean)[,2],
                p=aggregate(pi<.01~a+b+g+n,d3,mean)[,5],
                bf=aggregate(bfi>3~a+b+g+n,d3,mean)[,5],
                p1=aggregate(p1<.01~a+b+g+n,d3,mean)[,5],
                p2=aggregate(p2<.01~a+b+g+n,d3,mean)[,5])
round(cbind(t1[7:12,c(1,5:8)],t1[13:18,5:8],t1[19:24,5:8],t1[25:30,5:8],t1[31:36,5:8],t1[37:42,5:8]),3)[,c(2,6,10,14,18,22,3,7,11,15,19,23,4,8,12,16,20,5,9,13,17,21)]

#Paper 4, Table 2
d3b = d3
d3b$bfi = d3b$bfi>3
d3b$eq = case_when(
  d3b$n == "60" ~ mean(d3b$bfi[d3b$n==60&d3b$g<2]),
  d3b$n == "100" ~ mean(d3b$bfi[d3b$n==100&d3b$g<2]),
  d3b$n == "200" ~ mean(d3b$bfi[d3b$n==200&d3b$g<2]),
  d3b$n == "400" ~ mean(d3b$bfi[d3b$n==400&d3b$g<2]),
  d3b$n == "800" ~ mean(d3b$bfi[d3b$n==800&d3b$g<2]),
  d3b$n == "1600" ~ mean(d3b$bfi[d3b$n==1600&d3b$g<2]),
  d3b$n == "3200" ~ mean(d3b$bfi[d3b$n==3200&d3b$g<2])
)
d3b$pi = d3b$pi<d3b$eq

(t1b = data.frame(g=aggregate(bfi~a+b+g+n,d3b,mean)[,3],
                  n=aggregate(bfi~a+b+g+n,d3b,mean)[,4],
                  a=aggregate(bfi~a+b+g+n,d3b,mean)[,1],
                  b=aggregate(bfi~a+b+g+n,d3b,mean)[,2],
                  p=aggregate(pi~a+b+g+n,d3b,mean)[,5],
                  bf=aggregate(bfi~a+b+g+n,d3b,mean)[,5]))
# round(cbind(t1b[7:12,c(1,2,5,6)],t1b[13:18,c(1,2,5,6)],t1b[19:24,c(1,2,5,6)],t1b[25:30,c(1,2,5,6)],t1b[31:36,c(1,2,5,6)],t1b[37:42,c(1,2,5,6)]),3)
round(cbind(t1b[1:6,c(1,5,6)],t1b[13:18,c(5,6)],t1b[25:30,c(5,6)],t1b[37:42,c(5,6)]),3)

t2 = data.frame(g=aggregate(bfi>3~a+b+g+n,d3,mean)[,3],
                n=aggregate(bfi>3~a+b+g+n,d3,mean)[,4],
                a=aggregate(bfi>3~a+b+g+n,d3,mean)[,1],
                b=aggregate(bfi>3~a+b+g+n,d3,mean)[,2],
                bz=aggregate(bzi~a+b+g+n,d3,mean)[,5],
                bb=aggregate(bbi~a+b+g+n,d3,mean)[,5])

test = data.frame(foreach(i=1:10000,.combine="rbind") %dopar% {
  X = rep(c(0,1,0,1),each=n/4)
  G = rep(c(-.5,.5),each=n/2)
  Y = rep(0,n)
  Y[1:(3*n/4)] = rnorm(3*n/4)
  Y[(3*n/4+1):n] = rnorm(n/4,a,1)
  df = as.data.frame(cbind(Y,X,G))
  Z1 = lm(Y~X*G,df)
  Z2 = lm(Y~X+G,df)
  Z3 = lm(Y~X,df)
  bt = bartlett.test(Y~X,df)
  c(coefficients(Z1)[-1],coefficients(Z2)[-1],coefficients(Z3)[-1],summary(Z1)$coefficients[-1,4],summary(Z2)$coefficients[-1,4],
    summary(Z3)$coefficients[-1,4],bartlett.test(Y~X,df)$p.value)
})
names(test) = c("Z1X","Z1G","Z1XG","Z2X","Z2G","Z3X","Z1Xp","Z1Gp","Z1XGp","Z2Xp","Z2Gp","Z3Xp","btp")
test2 = test %>% mutate(across(ends_with("p"),~.<0.05))
apply(test2,2,mean)

n = 200
g = 4
d = pwr::pwr.t.test(n/2,NULL,.05,.8)$d*2
a = ifelse(g==1,.5*d,ifelse(g==2,.5*d,ifelse(g==3,.75*d,ifelse(g==4,d,ifelse(g==5,2.5*d,0)))))
b = ifelse(g==1,.5*d,ifelse(g==2,-.5*d,ifelse(g==3,-.25*d,ifelse(g==4,0,ifelse(g==5,1.5*d,0)))))
df = data.frame(Y = c(rep(0,n)))
df$Y[1:(n/4)] = rnorm(n/4)
df$Y[(n/4+1):(n/2)] = rnorm(n/4,b,1)
df$Y[(n/2+1):(3*n/4)] = rnorm(n/4)
df$Y[(3*n/4+1):n] = rnorm(n/4,a,1)
df$X = rep(c(0,1,0,1),each=n/4)
df$G = rep(c(-.5,.5),each=n/2)

# library(rms)
rm(dd);options(datadist=NULL)
dd=datadist(df2);options(datadist="dd")
m=ols(Y~X*G,df)
ma=anova(m)
ggplot(Predict(m,X,G),anova = ma, adj.subtitle=FALSE) + geom_smooth(method = "lm")

new = data.frame(foreach(i=1:10000,.combine="rbind") %dopar% {
  X = rep(c(0,1),each=n/2)
  Y = rep(0,n)
  Y[1:(n/2)] = rnorm(n/2)
  Y[(n/2+1):n] = rnorm(n/2,a/2,1)
  df = as.data.frame(cbind(Y,X))
  df2 = as.data.frame(matrix(rnorm(n*10),ncol = 10))
  df = bind_cols(df,df2)
  Z1 = lm(Y~X,df)
  df2 = df[,-2]
  df2$Y = as.numeric(df2$Y>mean(df2$Y))
  p1 = cor.test(df2$Y,df2$V1)$p.value
  p2 = cor.test(df2$Y,df2$V2)$p.value
  p3 = cor.test(df2$Y,df2$V3)$p.value
  p4 = cor.test(df2$Y,df2$V4)$p.value
  p5 = cor.test(df2$Y,df2$V5)$p.value
  p6 = cor.test(df2$Y,df2$V6)$p.value
  p7 = cor.test(df2$Y,df2$V7)$p.value
  p8 = cor.test(df2$Y,df2$V8)$p.value
  p9 = cor.test(df2$Y,df2$V9)$p.value
  p10 = cor.test(df2$Y,df2$V10)$p.value
  c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
})

table(as.numeric(foreach(i=1:10000) %dopar% {ks.test(c(rnorm(100),rnorm(100,mean=.8)),"pnorm",.4,1)$p.value})<.05)
ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + stat_function(fun = ~ dnorm(.x)+dnorm(.x,mean=2))

s = diag(20) + matrix(c(0,rep(0.223606,19),rep(c(0.223606,rep(0,19)),19)),nrow=20)
cov = data.frame(foreach(g=c(2,3,4,5),.combine="rbind") %:% foreach(be=c(0.2,0.5,0.8),.combine="rbind") %:% foreach(al=c(0.01),.combine="rbind") %:% foreach(n=200,.combine="rbind") %:% foreach(i=1:10000,.combine="rbind",.packages = "pwr") %dopar% {
  d = pwr.t.test(n/2,NULL,al,be)$d*2
  a = ifelse(g==1,.5*d,ifelse(g==2,.5*d,ifelse(g==3,.75*d,ifelse(g==4,d,ifelse(g==5,2.5*d,0)))))
  b = ifelse(g==1,.5*d,ifelse(g==2,-.5*d,ifelse(g==3,-.25*d,ifelse(g==4,0,ifelse(g==5,1.5*d,0)))))
  X = rep(c(0,1,0,1),each=n/4)
  G = rep(c(-.5,.5),each=n/2)
  df = as.data.frame(MASS::mvrnorm(n,rep(0,20),s))
  names(df) = c("Y","Z1","Z2","Z3","Z4","Z5","Z6","Z7","Z8","Z9","Z10","Z11","Z12","Z13","Z14","Z15","Z16","Z17","Z18","Z19")
  df$Y[(n/4+1):(n/2)] = df$Y[(n/4+1):(n/2)] + b
  df$Y[(3*n/4+1):n] = df$Y[(3*n/4+1):n] + a
  df = cbind(df,X,G)
  Z0 = lm(Y~X*G,df)
  Z1 = lm(Y~X*G+Z1,df)
  Z2 = lm(Y~X*G+Z1+Z2,df)
  Z3 = lm(Y~X*G+Z1+Z2+Z3,df)
  Z4 = lm(Y~X*G+Z1+Z2+Z3+Z4,df)
  Z5 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5,df)
  Z6 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6,df)
  Z7 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7,df)
  Z8 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8,df)
  Z9 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9,df)
  Z10 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10,df)
  Z11 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11,df)
  Z12 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12,df)
  Z13 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13,df)
  Z14 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13+Z14,df)
  Z15 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13+Z14+Z15,df)
  Z16 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13+Z14+Z15+Z16,df)
  Z17 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13+Z14+Z15+Z16+Z17,df)
  Z18 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13+Z14+Z15+Z16+Z17+Z18,df)
  Z19 = lm(Y~X*G+Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13+Z14+Z15+Z16+Z17+Z18+Z19,df)
  c(g,a,b,al,n,be,anova(Z0)[3,5],anova(Z1)[4,5],anova(Z2)[5,5],anova(Z3)[6,5],anova(Z4)[7,5],anova(Z5)[8,5],anova(Z6)[9,5],anova(Z7)[10,5],
    anova(Z8)[11,5],anova(Z9)[12,5],anova(Z10)[13,5],anova(Z11)[14,5],anova(Z12)[15,5],anova(Z13)[16,5],anova(Z14)[17,5],
    anova(Z15)[18,5],anova(Z16)[19,5],anova(Z17)[20,5],anova(Z18)[21,5],anova(Z19)[22,5])
})
names(cov) = c("g","a","b","al","n","be",paste0(seq(0,95,5),"%"))
cov2 = pivot_longer(cov,contains("%"),names_to = "adj",values_to = "p")
names(cov2) = c("g","a","b","al","n","be","adj","p")
cov3 = cov2
cov3$p[cov3$al==0.01] = cov3$p[cov3$al==0.01]<0.01
cov3 = cov3 %>% filter(al==0.01)
cov3$adj = as.numeric(gsub("%$","",cov3$adj))
cov3$al = as.factor(cov3$al)
cov3$be = as.factor(cov3$be)
cov4 = aggregate(p~al+be+adj,cov3,mean)
ggplot(cov4,aes(adj,p,colour=be,linetype=be)) + geom_line(size=1.25) + xlab("Percent of variance in outcome adjusted for") + 
  ylab("Adjusted power") + labs(linetype = "Power",colour = "Power") + xlim(0,100)
ggsave(filename = "Rplot01.png",width = 6, height = 4, device='png', dpi=1000)

# ggplot(cov3,aes(adj,p,colour=al,linetype=be)) + geom_smooth(se = F) + xlab("Percent of variance in outcome adjusted for") + 
#   ylab("Adjusted power") + scale_linetype_discrete(name="Power") + 
#   ggtitle("Power gained due to adjusting for covariates")

# Seems to be a test of covariate adjustment
# test = data.frame(foreach(i=1:10000,.combine="rbind") %dopar% {
#   X = rep(c(0,1,0,1),each=n/4)
#   G = rep(c(-.5,.5),each=n/2)
#   df = as.data.frame(MASS::mvrnorm(n,rep(0,20),s))
#   names(df) = c("Y","Z1","Z2","Z3","Z4")
#   df$Y[(n/4+1):(n/2)] = df$Y[(n/4+1):(n/2)] + b
#   df$Y[(3*n/4+1):n] = df$Y[(3*n/4+1):n] + a
#   df = cbind(df,X,G)
#   Z1 = lm(Y~X*G+Z1+Z2+Z3+Z4,df)
#   Z1b = lm(Y~X*G+Z1+Z2+Z3,df)
#   Z1c = lm(Y~X*G+Z1+Z2,df)
#   Z1d = lm(Y~X*G+Z1,df)
#   Z2 = lm(Y~X*G,df)
#   Z3 = lm(Y~X+G,df)
#   Z4 = lm(Y~X,df)
#   a12 = anova(Z1,Z2)$`Pr(>F)`[-1]
#   a23 = anova(Z2,Z3)$`Pr(>F)`[-1]
#   a34 = anova(Z3,Z4)$`Pr(>F)`[-1]
#   a1 = anova(Z1)
#   a1b = anova(Z1b)
#   a1c = anova(Z1c)
#   a1d = anova(Z1d)
#   a2 = anova(Z2)
#   a3 = anova(Z3)
#   a4 = anova(Z4)
#   c(coefficients(Z1)[-1],coefficients(Z2)[-1],coefficients(Z3)[-1],coefficients(Z4)[-1],summary(Z1)$coefficients[-1,4],
#     summary(Z2)$coefficients[-1,4],summary(Z3)$coefficients[-1,4],summary(Z4)$coefficients[-1,4],a12,a23,a34,
#     a1[7,5],a1b[6,5],a1c[5,5],a1d[4,5],a2[3,5])
# })
# names(test) = c("bz1.x","bz1.g","bz1.z1","bz1.z2","bz1.z3","bz1.z4","bz1.xg","bz2.x","bz2.g","bz2.xg","bz3.x","bz3.g","bz4.x",
#                 "pz1.x","pz1.g","pz1.z1","pz1.z2","pz1.z3","pz1.z4","pz1.xg","pz2.x","pz2.g","pz2.xg","pz3.x","pz3.g","pz4.x",
#                 "pz1v2","pz2v3","pz3v4","pz1.i","pz1b.i","pz1c.i","pz1d.i","pz2.i")

pa = data.frame(d=seq(-.5,2.5,.01),B=pwr.t.test(50,seq(-.5,2.5,.01),.01)$power)
pa2 = data.frame(g=seq(1,201,1),B2=pa$B[pa$d>=0.5],B1=pa$B[pa$d>=-0.5&pa$d<=1.5],Bx=pwr.t.test(100,seq(0,2,.01),.01)$power,Bi=rep(0.8,length(pa$B[pa$d>=0.5])))
lab = expression(atop(atop(d[S1]~"="~-0.5,d[S2]~"="~0.5),atop(d[MT]~"="~0,d[IT]~"="~1)),
                 atop(atop(d[S1]~"="~0,d[S2]~"="~1),atop(d[MT]~"="~0.5,d[IT]~"="~1)),
                 atop(atop(d[S1]~"="~0.5,d[S2]~"="~1.5),atop(d[MT]~"="~1,d[IT]~"="~1)),
                 atop(atop(d[S1]~"="~1,d[S2]~"="~2),atop(d[MT]~"="~1.5,d[IT]~"="~1)),
                 atop(atop(d[S1]~"="~1.5,d[S2]~"="~2.5),atop(d[MT]~"="~2,d[IT]~"="~1)))

ggplot(pa2) + annotate("rect",xmin=0,xmax=(pwr.t.test(50,NULL,.01,.8)$d*100-49),ymin=0,ymax=1,alpha=.2,fill="cyan") +
  annotate("rect",xmin=(pwr.t.test(50,NULL,.01,.8)$d*100-49),xmax=(pwr.t.test(50,NULL,.01,.8)$d*100+51),ymin=0,ymax=1,alpha=.2,fill="deepskyblue4") +
  annotate("rect",xmin=(pwr.t.test(50,NULL,.01,.8)$d*100+51),xmax=201,ymin=0,ymax=1,alpha=.2,fill="blue") +
  geom_line(aes(g,Bx, linetype = "Main effect test (MT)"),linewidth=1) + geom_line(aes(g,Bi, linetype = "Interaction test (IT)"),linewidth=1) +
  geom_line(aes(g,B2, linetype = "Group 2 only test (S2)"),linewidth=1) + geom_line(aes(g,B1, linetype = "Group 1 only test (S1)"),linewidth=1) +
  scale_x_continuous("Effect sizes",breaks = c(1,51,101,151,201),labels = lab) +
  scale_y_continuous("Power",breaks = seq(0,1,.1)) +
  scale_linetype_manual("Type of\nsignificance test",breaks = c("Main effect test (MT)","Interaction test (IT)","Group 1 only test (S1)","Group 2 only test (S2)"),
                        values = c("solid","dotted","dashed","1342"),) +
  theme(axis.text.x=element_text(size=14))
ggsave(filename = "Rplot03.png",width = 6, height = 4, device='png', dpi=1000)


# test3 = as.data.frame(foreach(n = 10000, .combine="rbind") %:% foreach(g = g_list, .combine="rbind") %:% foreach(i = 1:reps, .combine="rbind",.packages = c("pwr","brms","BAS","BayesFactor","MASS"))  %dopar% {
#   d = pwr.t.test(n/2,NULL,.01,.8)$d*2
#   a = ifelse(g==1,.5*d,ifelse(g==2,.5*d,ifelse(g==3,.75*d,ifelse(g==4,d,ifelse(g==5,2.5*d,0)))))
#   b = ifelse(g==1,.5*d,ifelse(g==2,-.5*d,ifelse(g==3,-.25*d,ifelse(g==4,0,ifelse(g==5,1.5*d,0)))))
#   X = rep(c(0,1,0,1),each=n/4)
#   G = rep(c(-0.5,0.5),each=n/2)
#   df = as.data.frame(mvrnorm(n,c(0),matrix(c(1),nrow=1)))
#   names(df) = c("Y")
#   df$Y[(n/4+1):(n/2)] = df$Y[(n/4+1):(n/2)] + b
#   df$Y[(3*n/4+1):n] = df$Y[(3*n/4+1):n] + a
#   df = cbind(df,X,G)
#   Z1 = lm(Y~X+G+X:G,df)
#   B = bas.lm(Y~X*G,df,prior="JZS",include.always = ~1+X+G)
#   Br = ifelse(summary(B)[4,2]==1,summary(B)[5,2]/summary(B)[5,3],summary(B)[5,3])
#   Bb = bas.lm(Y~X*G,df,prior="JZS",include.always = ~1+X*G)
#   c(i,g,n,a,b,Br,anova(Z1)[3,5],confint(coef(Bb))[4,c(3,1,2)])
# })

#Distribution of effect sizes
g1 = list()
g1$a = ggplot() + geom_function(fun=~(dnorm(.x,-.1)+dnorm(.x,.1))/2) + geom_function(fun = ~.25*dnorm(.x,-.1),aes(linetype="1")) + geom_function(fun = ~.25*dnorm(.x,.1),aes(linetype="2")) +
  scale_linetype_manual("Subgroups",values=c('longdash','dotted')) + xlim(-4,4) + ylim(0,.45) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g1$b = ggplot() + geom_function(fun=~(dnorm(.x,-.25)+dnorm(.x,.25))/2) + geom_function(fun = ~.25*dnorm(.x,-.25),aes(linetype="1")) + geom_function(fun = ~.25*dnorm(.x,.25),aes(linetype="2")) +
  scale_linetype_manual("Subgroups",values=c('longdash','dotted')) + xlim(-4,4) + ylim(0,.45) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g1$c = ggplot() + geom_function(fun=~(dnorm(.x,-.4)+dnorm(.x,.4))/2) + geom_function(fun = ~.25*dnorm(.x,-.4),aes(linetype="1")) + geom_function(fun = ~.25*dnorm(.x,.4),aes(linetype="2")) +
  scale_linetype_manual("Subgroups",values=c('longdash','dotted')) + xlim(-4,4) + ylim(0,.45) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g1$d = ggplot() + geom_function(fun=~(dnorm(.x,-.7)+dnorm(.x,.7))/2) + geom_function(fun = ~.25*dnorm(.x,-.7),aes(linetype="1")) + geom_function(fun = ~.25*dnorm(.x,.7),aes(linetype="2")) +
  scale_linetype_manual("Subgroups",values=c('longdash','dotted')) + xlim(-4,4) + ylim(0,.45) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g1$e = ggplot() + geom_function(fun=~(dnorm(.x,-.98)+dnorm(.x,.98))/2) + geom_function(fun = ~.25*dnorm(.x,-.98),aes(linetype="1")) + geom_function(fun = ~.25*dnorm(.x,.98),aes(linetype="2")) +
  scale_linetype_manual("Subgroups",values=c('longdash','dotted')) + xlim(-4,4) + ylim(0,.45) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
g1$f = ggplot() + geom_function(fun=~(dnorm(.x,-1.25)+dnorm(.x,1.25))/2) + geom_function(fun = ~.25*dnorm(.x,-1.25),aes(linetype="1")) + geom_function(fun = ~.25*dnorm(.x,1.25),aes(linetype="2")) +
  scale_linetype_manual("Subgroups",values=c('longdash','dotted')) + xlim(-4,4) + ylim(0,.45) +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
ggarrange(g1$a, g1$b, g1$c, g1$d, g1$e, g1$f, common.legend = TRUE, legend = "right", ncol = 3, nrow = 2,
          labels = c("A: d = 0.2", "B: d = 0.5", "C: d = 0.8", "D: d = 1.4", "E: d = 1.96", "F: d = 2.5"))
ggsave(filename = "Rplot05.png",width = 6, height = 4, device='png', dpi=1000)

#Simple slopes
g2 = list()
g2$a = ggplot(data.frame(Time=factor(rep(c("Pre","Post"),each=2),levels = c("Pre","Post")),y=c(3,5,6,8)),aes(Time,y)) + geom_point() +
  geom_segment(x=1,y=5,xend=2,yend=8,aes(colour="Subgroup 1\nEffect",linetype="Subgroup 1\nEffect"),linewidth=1) + geom_segment(x=1,y=3,xend=2,yend=6,aes(colour="Subgroup 2\nEffect",linetype="Subgroup 2\nEffect"),linewidth=1) +
  geom_segment(x=1,y=4,xend=2,yend=7,aes(colour="Average\nEffect",linetype="Average\nEffect"),linewidth=1) + ylim(2,8) +
  scale_colour_manual(name="Effect",values=c("black","blue","red")) + scale_linetype_manual(name="Effect",values=c("dotted","solid","solid")) +
  theme(axis.text.x=element_text(size=10),axis.title.x=element_blank()) + guides(colour = guide_legend(byrow = TRUE))
g2$b = ggplot(data.frame(Time=factor(rep(c("Pre","Post"),each=2),levels = c("Pre","Post")),y=c(5,5,6.5,8)),aes(Time,y)) + geom_point() +
  geom_segment(x=1,y=5,xend=2,yend=8,aes(colour="Subgroup 1\nEffect",linetype="Subgroup 1\nEffect"),linewidth=1) + geom_segment(x=1,y=5,xend=2,yend=6.5,aes(colour="Subgroup 2\nEffect",linetype="Subgroup 2\nEffect"),linewidth=1) +
  geom_segment(x=1,y=5,xend=2,yend=7.25,aes(colour="Average\nEffect",linetype="Average\nEffect"),linewidth=1) + ylim(2,8) +
  scale_colour_manual(name="Effect",values=c("black","blue","red")) + scale_linetype_manual(name="Effect",values=c("dotted","solid","solid")) +
  theme(axis.text.x=element_text(size=10),axis.title.x=element_blank()) + guides(colour = guide_legend(byrow = TRUE))
g2$c = ggplot(data.frame(Time=factor(rep(c("Pre","Post"),each=2),levels = c("Pre","Post")),y=c(5,5,5,8)),aes(Time,y)) + geom_point() +
  geom_segment(x=1,y=5,xend=2,yend=8,aes(colour="Subgroup 1\nEffect",linetype="Subgroup 1\nEffect"),linewidth=1) + geom_segment(x=1,y=5,xend=2,yend=5,aes(colour="Subgroup 2\nEffect",linetype="Subgroup 2\nEffect"),linewidth=1) +
  geom_segment(x=1,y=5,xend=2,yend=6.5,aes(colour="Average\nEffect",linetype="Average\nEffect"),linewidth=1) + ylim(2,8) +
  scale_colour_manual(name="Effect",values=c("black","blue","red")) + scale_linetype_manual(name="Effect",values=c("dotted","solid","solid")) +
  theme(axis.text.x=element_text(size=10)) + guides(colour = guide_legend(byrow = TRUE))
g2$d = ggplot(data.frame(Time=factor(rep(c("Pre","Post"),each=2),levels = c("Pre","Post")),y=c(5,5,2,8)),aes(Time,y)) + geom_point() +
  geom_segment(x=1,y=5,xend=2,yend=8,aes(colour="Subgroup 1\nEffect",linetype="Subgroup 1\nEffect"),linewidth=1) + geom_segment(x=1,y=5,xend=2,yend=2,aes(colour="Subgroup 2\nEffect",linetype="Subgroup 2\nEffect"),linewidth=1) +
  geom_segment(x=1,y=5,xend=2,yend=5,aes(colour="Average\nEffect",linetype="Average\nEffect"),linewidth=1) + ylim(2,8) +
  scale_colour_manual(name="Effect",values=c("black","blue","red")) + scale_linetype_manual(name="Effect",values=c("dotted","solid","solid")) +
  theme(axis.text.x=element_text(size=10)) + guides(colour = guide_legend(byrow = TRUE))
ggarrange(g2$a, g2$b, g2$c, g2$d, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend = "right", ncol = 2, nrow = 2)
ggsave(filename = "Rplot04.png",width = 6, height = 4, device='png', dpi=1000)

#Overview table
t2 = round(data.frame(me=c(0,1,0,.5,1,4),
                      ie=c(0,0,2,2,2,2),
                      px=aggregate(px<.01~a+b+g,d3[d3$n==200,],mean)[,4]*100,
                      pg=aggregate(pg<.01~a+b+g,d3[d3$n==200,],mean)[,4]*100,
                      pi=aggregate(pi<.01~a+b+g,d3[d3$n==200,],mean)[,4]*100,
                      dx=(aggregate(bfi>3~a+b+g,d3[d3$n==200,],mean)[,1]+aggregate(bfi>3~a+b+g,d3[d3$n==200,],mean)[,2])/2,
                      di=aggregate(bfi>3~a+b+g,d3[d3$n==200,],mean)[,1]-aggregate(bfi>3~a+b+g,d3[d3$n==200,],mean)[,2],
                      d1=aggregate(bfi>3~a+b+g,d3[d3$n==200,],mean)[,2],
                      d2=aggregate(bfi>3~a+b+g,d3[d3$n==200,],mean)[,1]),2)
#round(cbind(c(0,1,0,.5,1,4),c(0,0,2,2,2,2),t1$p[t1$n==200]*100,(t1$a[t1$n==200]+t1$b[t1$n==200])/2,t1$a[t1$n==200]+t1$b[t1$n==200],t1$a[t1$n==200],t1$b[t1$n==200]),2)

# Power graph
i_list = c(75,100,150,200,300,400,600,800,1200,1600,2000)
d_list = 0.5
#for(l in 1:20) {i_list[l]= 1.43^l+8.75}
reps=10000
T=as.data.frame(cbind(n=rep(NA,length(i_list)),d=rep(NA,length(i_list))))
for(i in 1:length(i_list)) {
  T$n[i] = i_list[i]
  T$d[i] = pwr.t.test(i_list[i],NULL,.01,.8)$d
}
T$n = T$n*2
T$di = T$d*2
p_list = c(pwr.t.test(NULL,d_list/4,.01,.8)$n,pwr.t.test(NULL,d_list/2,.01,.8)$n,pwr.t.test(NULL,d_list,.01,.8)$n,pwr.t.test(NULL,d_list*2,.01,.8)$n)
ggplot() +
  geom_line(data=T,aes(n,d,color="Main"),size=1.25) + 
  geom_line(data=T,aes(n,di,color="Int."),size=1.25) +
  geom_point(data=t1[t1$g==4&t1$n>150,],aes(n,a, shape="Simulation estimates"),size=2) +
  scale_x_continuous(breaks = c(200,500,1000,2000,3000), trans = scales::log10_trans()) +
  scale_shape_discrete("") +
  annotation_logticks(sides = "b") +
  labs(colour = "Type of effect") +
  xlab("Sample size (log scale)") + ylab("Effect size (Cohen's d)")
ggsave(filename = "Linear.png",width = 6, height = 4, device='png', dpi=1000)

#Power subgroup ratios
r_list = c(2,seq(10,100,10))
test4 = as.data.frame(foreach(r = r_list, .combine="rbind") %:% foreach(g = 2:5, .combine="rbind") %:% foreach(i = 1:100000, .combine="rbind",.packages = c("MASS"))  %dopar% {
  n = 200
  d = pwr.t.test(n/2,NULL,.01,.8)$d*2
  a = ifelse(g==1,.5*d,ifelse(g==2,.5*d,ifelse(g==3,.75*d,ifelse(g==4,d,ifelse(g==5,2.5*d,0)))))
  b = ifelse(g==1,.5*d,ifelse(g==2,-.5*d,ifelse(g==3,-.25*d,ifelse(g==4,0,ifelse(g==5,1.5*d,0)))))
  X = rep(c(0,1),n/2)
  G = c(rep(-0.5,r),rep(0.5,n-r))
  df = as.data.frame(mvrnorm(n,c(0),matrix(c(1),nrow=1)))
  names(df) = c("Y")
  df = cbind(df,X,G)
  df$Y[df$X==1&df$G==-.5] = df$Y[df$X==1&df$G==-.5] + b
  df$Y[df$X==1&df$G==.5] = df$Y[df$X==1&df$G==.5] + a
  Z1 = lm(Y~X+G+X:G,df)
  c(i,g,r,a,b,anova(Z1)[1:3,5])
})
names(test4) = c("i","g","r","a","b","px","pg","pi")
test4$si = test4$pi<0.01

ggplot() + geom_line(data=data.frame(x=seq(198,2,-2)/2,y=pwr.t2n.test(seq(2,198,2),seq(198,2,-2),pwr.t.test(100,NULL,.01,.8)$d,.01,NULL)$power),aes(x,y),size=1.25,colour="#F8766D") + ylim(0,1) + xlab("% of total sample in one subgroup") + ylab("Power") + labs(title="Power curve for different subgroup ratios") +
  geom_point(data=aggregate(si~r,test4,mean),aes(r/2,si,colour = "Simulation results"),size=2) + scale_colour_manual("",values=c("black"))
ggsave(filename = "Rplot07.png",width = 6, height = 4, device='png', dpi=1000)

#Power treatment ratios
r_list = c(2,seq(10,100,10))
test5 = as.data.frame(foreach(r = r_list, .combine="rbind") %:% foreach(g = 2:5, .combine="rbind") %:% foreach(i = 1:100000, .combine="rbind",.packages = c("MASS"))  %dopar% {
  n = 200
  d = pwr.t.test(n/2,NULL,.01,.8)$d*2
  a = ifelse(g==1,.5*d,ifelse(g==2,.5*d,ifelse(g==3,.75*d,ifelse(g==4,d,ifelse(g==5,2.5*d,0)))))
  b = ifelse(g==1,.5*d,ifelse(g==2,-.5*d,ifelse(g==3,-.25*d,ifelse(g==4,0,ifelse(g==5,1.5*d,0)))))
  X = c(rep(0,r),rep(1,n-r))
  G = rep(c(-0.5,0.5),n/2)
  df = as.data.frame(mvrnorm(n,c(0),matrix(c(1),nrow=1)))
  names(df) = c("Y")
  df = cbind(df,X,G)
  df$Y[df$X==1&df$G==-.5] = df$Y[df$X==1&df$G==-.5] + b
  df$Y[df$X==1&df$G==.5] = df$Y[df$X==1&df$G==.5] + a
  Z1 = lm(Y~X+G+X:G,df)
  c(i,g,r,a,b,anova(Z1)[1:3,5])
})
names(test5) = c("i","g","r","a","b","px","pg","pi")
test5$si = test5$pi<0.01

ggplot() + geom_line(data=data.frame(x=seq(198,2,-2)/2,y=pwr.t2n.test(seq(2,198,2),seq(198,2,-2),pwr.t.test(100,NULL,.01,.8)$d,.01,NULL)$power),aes(x,y),size=1.25,colour="#F8766D") + ylim(0,1) + xlab("% of total sample in treatment group") + ylab("Power") + labs(title="Power curve for different treatment ratios") +
  geom_point(data=aggregate(si~r,test5,mean),aes(r/2,si,colour = "Simulation results"),size=2) + scale_colour_manual("",values=c("black"))
ggsave(filename = "Rplot08.png",width = 6, height = 4, device='png', dpi=1000)

#Double imbalance
r_list = c(seq(.1,0.5,.1))
r2_list = c(seq(.1,0.5,.1))
df_dbl = as.data.frame(foreach(r = r_list, .combine="rbind") %:% foreach(r2 = r2_list, .combine="rbind") %:% foreach(i = 1:10000, .combine="rbind", .options.future = list(packages = c("MASS"),seed = TRUE))  %dofuture% {
  n = 4000
  a = pwr.t.test(n/2,NULL,.01,.8)$d*2
  b = pwr.t2n.test(n*r,n*(1-r),NULL,.01,.8)$d*2
  c = pwr.t2n.test(n*r2,n*(1-r2),NULL,.01,.8)$d*2
  d = b*c/a
  X = rbinom(n,1,r)
  G = rbinom(n,1,r2) - 0.5
  df = as.data.frame(mvrnorm(n,c(0),matrix(c(1),nrow=1)))
  names(df) = c("Y")
  df = cbind(df,X,G)
  df$Y[df$X==1&df$G==.5] = df$Y[df$X==1&df$G==.5] + d
  Z1 = lm(Y~X+G+X:G,df)
  c(i,r,r2,anova(Z1)[1:3,5],coef(Z1)[-1])
})
names(df_dbl) = c("i","r","r2","px","pg","pi","bx","bg","bi")
df_dbl$si = df_dbl$pi<0.01

df_dbl %>% group_by(r,r2) %>% summarise(p = mean(si), b = sprintf("%0.8f",mean(bi))) %>% print(n=1000)

#App deployment
rsconnect::deployApp(appDir = "./Interaction_power_calculator")

# Archiving
save(d2,d3,d4,test,test3,test4,test5,new,cov,file = "./PhD Backup.RData",compress = "xz")
rm(d2,d3,d4,test,test3,test4,test5,new,cov)
load(file = "./PhD Backup.RData")
rm(d3b,d4b,test2,df2,cov2,cov3)
