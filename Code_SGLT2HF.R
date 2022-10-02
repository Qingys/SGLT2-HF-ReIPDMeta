
###IPD --------------------------------------------------------------------
# library(tidyverse)
# library(IPDfromKM)
# library(survminer)
# ipd <- extract_ipd('points_0.txt', 'points_1.txt',
#                    timeSeq = c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36),
#                    numRisk0 = c(2991, 2960, 2923, 2883, 2849, 2656, 2302, 2046, 1738, 1459, 1107, 781, 471),
#                    numRisk1 = c(2997, 2974, 2930, 2891, 2847, 2641, 2287, 2022, 1725, 1469, 1118, 794, 462),
#                    events0 = 244, events1 = 219)
# write.csv(ipd, 'EMPEROR-P_CVdeath.csv')
# 
# 
# ###save dt
# ipd$study <- 'EMPEROR-P'
# ipd$outcome <- 'CVdeath'
# ipd$EF <- 'HFpEF'
# 
# 
# #dt <- ipd
# dt <- bind_rows(dt, ipd)
# ###
# #dt[study=='EMPEROR-P'][outcome=='kidney']$time <- ipd$time
# #dt[study=='EMPEROR-P'][outcome=='kidney']$status <- ipd$status
# #dt[study=='EMPEROR-P'][outcome=='kidney']$treat <- ipd$treat
# 
# write.csv(dt, 'ipd.csv')

##########################two-stage meta-analysis#################
library(tidyverse)
library(broom)
library(survival)
library(boot)
library(flexsurv)

######individual study tRR
###function
run_study_tRR <- function(outcomei, studyi, timei){
  dt_outcome_study <- dt %>% filter(outcome==outcomei,
                                    study==studyi)
  boot_logtrr <- function(dt_b, index){
    dtsub <- dt_b[index,]
    ###flexible hazard regression
    fit_flex_HHF <- flexsurvspline(Surv(time, status) ~ treat + 
                                     gamma1(treat),
                                   data = dtsub, k=2)
    pred <- predict(fit_flex_HHF, 
                    type = 'survival',
                    newdata = tibble(treat=c(0,1)), 
                    times = c(6,12,24))
    pred_tseq <- tibble(
      treat=c(rep(0,3),rep(1,3)),
      time=rep(c(6,12,24),2),
      surv=c(pred$.pred[[1]][,2],pred$.pred[[2]][,2])
    )
    pred_tseq <- pred_tseq %>% 
      mutate(risk=1-surv) %>% 
      group_by(time) %>% 
      mutate(tRR=risk[2]/risk[1]) %>% 
      mutate(logtRR=log(tRR)) %>% 
      ungroup() %>% 
      filter(treat==1) %>% 
      select(time,logtRR)
    logtRR <- pred_tseq %>% filter(time==timei) %>% pull(logtRR)
    return(logtRR)
  }
  
  ###boot
  boot_i <- boot(dt_outcome_study, boot_logtrr, 
                 R=500, parallel = 'multicore',
                 ncpus = 4)
  sum_tb <- summary(boot_i)
  return(sum_tb)
}

###map data
#run_study_tRR('HHF', 'DAPA-HF', 24)

pred_list <- list(
  outcomei=c(rep('HHF',9),rep('ALLcause',9),rep('CVdeath',9),rep('kidney',9)),
  studyi=rep(rep(c('DAPA-HF','EMPEROR-P','EMPEROR-R'),3),4),
  timei=rep(c(rep(6,3),rep(12,3),rep(24,3)),4)
)

pred_tb <- pmap_dfr(pred_list, run_study_tRR)

pred_tb <- pred_tb %>% 
  mutate(outcome=pred_list$outcomei,
  study=pred_list$studyi,
  time=pred_list$timei) %>% 
  mutate(tRR=exp(original), 
         lb=exp(original-1.96*bootSE),
         ub=exp(original+1.96*bootSE))
sink('results.txt')
print(pred_tb)
sink()



###############################IPD meta-analysis
##########################one stage meta analysis##################
library(tidyverse)
library(survival)
library(survminer)
library(mexhaz)
library(boot)

###read data
dt <- read_csv('ipd.csv')

###1 HHF-----------------------------------------------------------------------
###fit HHF
dt_HHF <- dt %>% filter(outcome=='HHF') %>% as.data.frame()

###cox
#PH assumption
dt_HHF %>% coxph(Surv(time, status) ~ treat+strata(study), data =.) %>% cox.zph()
dt_HHF %>% coxph(Surv(time, status) ~ treat, data =.) %>% cox.zph()

###km
#naive
fit_km_HHF <- survfit(Surv(time, status) ~ treat, 
                          data = dt_HHF)

tiff('HHF.tiff', res = 300, width = 300*10, height = 300*7)
ggsurvplot(fit_km_HHF, 
           fun = 'event',
           palette = 'lancet',
           data = dt_HHF,
           risk.table = T,
           censor=F,
           conf.int = F,
           conf.int.style="step",
           xlim=c(0, 24),
           ylim=c(0, 0.2),
           break.x.by = 3,
           xlab='Months since randomization',
           ylab='Cumulative incidence',
           legend.labs = c('Placebo', 'SGLT2 inhibitors'),
           legend = c(0.1,0.9),
           legend.title='',
           surv.scale='percent',
           base_size = 11,
           ggtheme = theme_survminer(
             font.x = 12,
             font.y = 12,
             font.main = 12
           ),
           title=''
           )
dev.off()


###fit one-stage model
fit_frailty_HHF <- mexhaz(Surv(time, status) ~ treat + nph(treat), data = dt_HHF,
                          random = 'study', base = 'exp.ns', knots = median(dt_HHF$time))

###pred RR
fun_predRR <- function(time1){
  pred <- predict(fit_frailty_HHF, time.pts = time1, cluster = NULL,
                  data.val = data.frame(treat=c(0, 1)))
  tibble1 <- tibble(time=pred$results$time.pts[1], 
                    tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
  return(tibble1)
}
timeseq <- seq(0.1, 24, 0.1)
abs_HHF <- map_dfr(timeseq, fun_predRR)

###boot
fun_boot <- function(dat, index){
  library(mexhaz)
  library(purrr)
  datsub <- dat[index,]
  fit_frailty_HHF <- mexhaz(Surv(time, status) ~ treat + nph(treat), 
                            data = datsub,
                            random = 'study', base = 'exp.ns', knots = median(dt_HHF$time))
  fun_predRR <- function(time1){
    pred <- predict(fit_frailty_HHF, time.pts = time1, cluster = NULL,
                    data.val = data.frame(treat=c(0, 1)))
    tRRdf <- data.frame(time=pred$results$time.pts[1], 
                        tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
    return(tRRdf)
  }
  timeseq <- seq(0.1, 24, 0.1)
  abs_HHF <- map_dfr(timeseq, fun_predRR)
  tRR <- abs_HHF$tRR
  return(tRR)
}
###boot
boot_HHF <- boot(dt_HHF, fun_boot, R=500, parallel = 'multicore', ncpus = 4)
###boot CI
abs_HHF_CI <- apply(boot_HHF$t, 2, function(x) quantile(x, c(0.025, 0.975))) %>% 
  t() %>% as.data.frame()
names(abs_HHF_CI) <- c('lb', 'ub')
abs_HHF_CI$est <- boot_HHF$t0
abs_HHF_CI$time <- timeseq

####time-varying RR
p_tRR_HHF <- abs_HHF_CI %>% 
  ggplot(aes(time, est)) +
  geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.35, fill='#006d2c') +
  geom_line(size=1.5, color='#006d2c') +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24), limits = c(0, 24)) +
  scale_y_continuous(breaks = c(0.5, 0.6, 0.8, 1), trans = 'log',
                     limits = c(0.46, 1)) +
  geom_hline(yintercept = 1, color='#bdbdbd', size=1) +
  geom_hline(yintercept = 0.73, color='#a50f15', size=1, linetype='dashed') +
  labs(x='Months since randomization', y='Time-varying risk ratio') +
  theme_bw()
print(p_tRR_HHF)
ggsave('tRR_HHF.tiff', dpi = 300, width = 8, height = 6, plot = p_tRR_HHF)


###2 Allcause---------------------------------------------------------------------
###fit Allcause
dt_AC <- dt %>% filter(outcome=='ALLcause') %>% as.data.frame()

###cox
#PH assumption
coxph(Surv(time, status) ~ treat+strata(study), data = dt_AC) %>% cox.zph()

###km
#naive
fit_km_AC <- survfit(Surv(time, status) ~ treat, 
                  data = dt_AC)

tiff('ALLcause.tiff', res = 300, width = 300*10, height = 300*7)
ggsurvplot(fit_km_AC, 
           fun = 'event',
           palette = 'lancet',
           data = dt_AC,
           risk.table = T,
           censor=F,
           conf.int = F,
           conf.int.style="step",
           xlim=c(0, 24),
           ylim=c(0, 0.3),
           break.x.by = 3,
           xlab='Months since randomization',
           ylab='Cumulative incidence',
           legend.labs = c('Placebo', 'SGLT2 inhibitors'),
           legend = c(0.1,0.9),
           legend.title='',
           surv.scale='percent',
           base_size = 11,
           ggtheme = theme_survminer(
             font.x = 12,
             font.y = 12,
             font.main = 12
           ),
           title=''
)
dev.off()


###fit one-stage model
fit_frailty_AC <- mexhaz(Surv(time, status) ~ treat + nph(treat), data = dt_AC,
                          random = 'study', base = 'exp.ns', knots = median(dt_AC$time))

###pred RR
fun_predRR <- function(time1){
  pred <- predict(fit_frailty_AC, time.pts = time1, cluster = NULL,
                  data.val = data.frame(treat=c(0, 1)))
  tibble1 <- tibble(time=pred$results$time.pts[1], 
                    tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
  return(tibble1)
}
timeseq <- seq(0.1, 24, 0.1)
abs_AC <- map_dfr(timeseq, fun_predRR)

###boot
fun_boot <- function(dat, index){
  library(mexhaz)
  library(purrr)
  datsub <- dat[index,]
  fit_frailty_AC <- mexhaz(Surv(time, status) ~ treat + nph(treat), 
                            data = datsub,
                            random = 'study', base = 'exp.ns', knots = median(dt_AC$time))
  fun_predRR <- function(time1){
    pred <- predict(fit_frailty_AC, time.pts = time1, cluster = NULL,
                    data.val = data.frame(treat=c(0, 1)))
    tRRdf <- data.frame(time=pred$results$time.pts[1], 
                        tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
    return(tRRdf)
  }
  timeseq <- seq(0.1, 24, 0.1)
  abs_AC <- map_dfr(timeseq, fun_predRR)
  tRR <- abs_AC$tRR
  return(tRR)
}
###boot
boot_AC <- boot(dt_AC, fun_boot, R=500, parallel = 'multicore', ncpus = 4)
###boot CI
abs_AC_CI <- apply(boot_AC$t, 2, function(x) quantile(x, c(0.025, 0.975))) %>% 
  t() %>% as.data.frame()
names(abs_AC_CI) <- c('lb', 'ub')
abs_AC_CI$est <- boot_AC$t0
abs_AC_CI$time <- timeseq

####time-varying RR
p_tRR_AC <- abs_AC_CI %>% 
  ggplot(aes(time, est)) +
  geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.35, fill='#006d2c') +
  geom_line(size=1.5, color='#006d2c') +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24), limits = c(0, 24)) +
  scale_y_continuous(breaks = c(0.7, 0.8, 1), trans = 'log',
                     limits = c(0.6, 1.2)) +
  geom_hline(yintercept = 1, color='#bdbdbd', size=1) +
  geom_hline(yintercept = 0.92, color='#a50f15', size=1, linetype='dashed') +
  labs(x='Months since randomization', y='Time-varying risk ratio') +
  theme_bw()
print(p_tRR_AC)
ggsave('tRR_AC.tiff', dpi = 300, width = 8, height = 6, plot = p_tRR_AC)




###3 kidney----------------------------------------------------------------------
######################fit kidney
dt_kidney <- dt %>% filter(outcome=='kidney') %>% as.data.frame()

###cox
#PH assumption
coxph(Surv(time, status) ~ treat+strata(study), data=dt_kidney) %>% cox.zph()

###km
#naive
fit_km_kidney <- survfit(Surv(time, status) ~ treat, 
                      data = dt_kidney)

tiff('kidney.tiff', res = 300, width = 300*10, height = 300*7)
ggsurvplot(fit_km_kidney, 
           fun = 'event',
           palette = 'lancet',
           data = dt_kidney,
           risk.table = T,
           censor=F,
           conf.int = F,
           conf.int.style="step",
           xlim=c(0, 24),
           ylim=c(0, 0.1),
           break.x.by = 3,
           xlab='Months since randomization',
           ylab='Cumulative incidence',
           legend.labs = c('Placebo', 'SGLT2 inhibitors'),
           legend = c(0.1,0.9),
           legend.title='',
           surv.scale='percent',
           base_size = 11,
           ggtheme = theme_survminer(
             font.x = 12,
             font.y = 12,
             font.main = 12
           ),
           title=''
)
dev.off()


###fit one-stage model
fit_frailty_kidney <- mexhaz(Surv(time, status) ~ treat + nph(treat), data = dt_kidney,
                         random = 'study', base = 'exp.ns', knots = median(dt_kidney$time))

###pred RR
fun_predRR <- function(time1){
  pred <- predict(fit_frailty_kidney, time.pts = time1, cluster = NULL,
                  data.val = data.frame(treat=c(0, 1)))
  tibble1 <- tibble(time=pred$results$time.pts[1], 
                    tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
  return(tibble1)
}
timeseq <- seq(0.1, 24, 0.1)
abs_kidney <- map_dfr(timeseq, fun_predRR)

###boot
fun_boot <- function(dat, index){
  library(mexhaz)
  library(purrr)
  datsub <- dat[index,]
  fit_frailty_kidney <- mexhaz(Surv(time, status) ~ treat + nph(treat), 
                           data = datsub,
                           random = 'study', base = 'exp.ns', knots = median(dt_kidney$time))
  fun_predRR <- function(time1){
    pred <- predict(fit_frailty_kidney, time.pts = time1, cluster = NULL,
                    data.val = data.frame(treat=c(0, 1)))
    tRRdf <- data.frame(time=pred$results$time.pts[1], 
                        tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
    return(tRRdf)
  }
  timeseq <- seq(0.1, 24, 0.1)
  abs_kidney <- map_dfr(timeseq, fun_predRR)
  tRR <- abs_kidney$tRR
  return(tRR)
}
###boot
boot_kidney <- boot(dt_kidney, fun_boot, R=500, parallel = 'multicore', ncpus = 4)
###boot CI
abs_kidney_CI <- apply(boot_kidney$t, 2, function(x) quantile(x, c(0.025, 0.975))) %>% 
  t() %>% as.data.frame()
names(abs_kidney_CI) <- c('lb', 'ub')
abs_kidney_CI$est <- boot_kidney$t0
abs_kidney_CI$time <- timeseq

####time-varying RR
p_tRR_kidney <- abs_kidney_CI %>% 
  ggplot(aes(time, est)) +
  geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.35, fill='#006d2c') +
  geom_line(size=1.5, color='#006d2c') +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24), limits = c(0, 24)) +
  scale_y_continuous(breaks = c(0.5, 0.7, 1, 1.5, 2), trans = 'log',
                     limits = c(0.45, 2.02)) +
  geom_hline(yintercept = 1, color='#bdbdbd', size=1) +
  geom_hline(yintercept = 0.73, color='#a50f15', size=1, linetype='dashed') +
  labs(x='Months since randomization', y='Time-varying risk ratio') +
  theme_bw()
print(p_tRR_kidney)
ggsave('tRR_kidney.tiff', dpi = 300, width = 8, height = 6, plot = p_tRR_kidney)



###4 CV death-----------------------------------------------------------------------
######################fit CV death
dt_CVdeath <- dt %>% filter(outcome=='CVdeath') %>% as.data.frame()

###cox
#PH assumption
coxph(Surv(time, status) ~ treat+strata(study), data=dt_CVdeath) %>% cox.zph()

###km
#naive
fit_km_CVdeath <- survfit(Surv(time, status) ~ treat, 
                         data = dt_CVdeath)

tiff('CVdeath.tiff', res = 300, width = 300*10, height = 300*7)
ggsurvplot(fit_km_CVdeath, 
           fun = 'event',
           palette = 'lancet',
           data = dt_CVdeath,
           risk.table = T,
           censor=F,
           conf.int = F,
           conf.int.style="step",
           xlim=c(0, 24),
           ylim=c(0, 0.2),
           break.x.by = 3,
           xlab='Months since randomization',
           ylab='Cumulative incidence',
           legend.labs = c('Placebo', 'SGLT2 inhibitors'),
           legend = c(0.1,0.9),
           legend.title='',
           surv.scale='percent',
           base_size = 11,
           ggtheme = theme_survminer(
             font.x = 12,
             font.y = 12,
             font.main = 12
           ),
           title=''
)
dev.off()


###fit one-stage model
fit_frailty_CVdeath <- mexhaz(Surv(time, status) ~ treat + nph(treat), data = dt_CVdeath,
                             random = 'study', base = 'exp.ns', knots = median(dt_CVdeath$time))

###pred RR
fun_predRR <- function(time1){
  pred <- predict(fit_frailty_CVdeath, time.pts = time1, cluster = NULL,
                  data.val = data.frame(treat=c(0, 1)))
  tibble1 <- tibble(time=pred$results$time.pts[1], 
                    tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
  return(tibble1)
}
timeseq <- seq(0.1, 24, 0.1)
abs_CVdeath <- map_dfr(timeseq, fun_predRR)

###boot
fun_boot <- function(dat, index){
  library(mexhaz)
  library(purrr)
  datsub <- dat[index,]
  fit_frailty_CVdeath <- mexhaz(Surv(time, status) ~ treat + nph(treat), 
                               data = datsub,
                               random = 'study', base = 'exp.ns', knots = median(dt_CVdeath$time))
  fun_predRR <- function(time1){
    pred <- predict(fit_frailty_CVdeath, time.pts = time1, cluster = NULL,
                    data.val = data.frame(treat=c(0, 1)))
    tRRdf <- data.frame(time=pred$results$time.pts[1], 
                        tRR=(1-pred$results$surv[2])/(1-pred$results$surv[1]))
    return(tRRdf)
  }
  timeseq <- seq(0.1, 24, 0.1)
  abs_CVdeath <- map_dfr(timeseq, fun_predRR)
  tRR <- abs_CVdeath$tRR
  return(tRR)
}
###boot
boot_CVdeath <- boot(dt_CVdeath, fun_boot, R=500, parallel = 'multicore', ncpus = 4)
###boot CI
abs_CVdeath_CI <- apply(boot_CVdeath$t, 2, function(x) quantile(x, c(0.025, 0.975))) %>% 
  t() %>% as.data.frame()
names(abs_CVdeath_CI) <- c('lb', 'ub')
abs_CVdeath_CI$est <- boot_CVdeath$t0
abs_CVdeath_CI$time <- timeseq

####time-varying RR
p_tRR_CVdeath <- abs_CVdeath_CI %>% 
  ggplot(aes(time, est)) +
  geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.35, fill='#006d2c') +
  geom_line(size=1.5, color='#006d2c') +
  scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18, 21, 24), limits = c(0, 24)) +
  scale_y_continuous(breaks = c(0.7, 0.8, 1), trans = 'log',
                     limits = c(0.6, 1.2)) +
  geom_hline(yintercept = 1, color='#bdbdbd', size=1) +
  geom_hline(yintercept = 0.88, color='#a50f15', size=1, linetype='dashed') +
  labs(x='Months since randomization', y='Time-varying risk ratio') +
  theme_bw()
print(p_tRR_CVdeath)
ggsave('tRR_CVdeath.tiff', dpi = 300, width = 8, height = 6, plot = p_tRR_CVdeath)



###5 events ----------------------------------------------------------------------
######risk and event
library(tidyverse)
library(survival)
library(riskRegression)
library(broom)

fun_predictrisk <- function(outcomei, studyi, timei){
  dt_sub <- dt %>% filter(outcome==outcomei, study==studyi)
  ###model
  fit <- survfit(Surv(time, status) ~ treat, data = dt_sub)
  tb <- tidy(fit) %>% 
    filter(time<=timei) %>% 
    group_by(strata) %>% 
    summarise(time=last(time), event=last(cumsum(n.event))) %>% 
    ungroup()
  tb$outcome <- outcomei
  tb$study <- studyi
  return(tb)
}

###test
fun_predictrisk(outcomei = 'ALLcause', studyi = 'DAPA-HF', timei = 24)

###map
map_list <- list(outcomei=c(rep('HHF', 9), rep('ALLcause',9), rep('kidney',9), rep('CVdeath',9)),
                 studyi=rep(rep(c('DAPA-HF', 'EMPEROR-P', 'EMPEROR-R'),3),4),
                 timei=rep(c(rep(6,3), rep(12,3), rep(24,3)),4)
                 )
map_list

map_tb <- pmap_dfr(map_list, fun_predictrisk)

sink('events.txt')
map_tb %>% print(n=100)
sink()





###6 robust test for NPH ---------------------------------------------------
library(CauchyCP)

nph_HHF <- CauchyCP(time = dt_HHF$time, status = dt_HHF$status, x = dt_HHF$treat,
                    cutpoints = seq(3, 24, 3))

nph_Allcause <- CauchyCP(time = dt_AC$time, status = dt_AC$status, x = dt_AC$treat,
                    cutpoints = seq(3, 24, 3))

nph_kidney <- CauchyCP(time = dt_kidney$time, status = dt_kidney$status, x = dt_kidney$treat,
                         cutpoints = seq(3, 24, 3))

nph_CVdeath <- CauchyCP(time = dt_CVdeath$time, status = dt_CVdeath$status, x = dt_CVdeath$treat,
                       cutpoints = seq(3, 24, 3))

bind_rows(
tibble(outcome='HHF', time=seq(3, 24, 3), 
       HR_before=nph_HHF$hrs$before, HR_after=nph_HHF$hrs$after,
       Pvalue=nph_HHF$p.unadj, Poverall=nph_HHF$pval),
tibble(outcome='Allcause', time=seq(3, 24, 3), 
       HR_before=nph_Allcause$hrs$before, HR_after=nph_Allcause$hrs$after,
       Pvalue=nph_Allcause$p.unadj, Poverall=nph_Allcause$pval),
tibble(outcome='kidney', time=seq(3, 24, 3), 
       HR_before=nph_kidney$hrs$before, HR_after=nph_kidney$hrs$after,
       Pvalue=nph_kidney$p.unadj, Poverall=nph_kidney$pval),
tibble(outcome='CVdeath', time=seq(3, 24, 3), 
       HR_before=nph_CVdeath$hrs$before, HR_after=nph_CVdeath$hrs$after,
       Pvalue=nph_CVdeath$p.unadj, Poverall=nph_CVdeath$pval)) %>% 
  write_csv('nph_test.csv')


###7 print one-stage results -------------------------------------------------

bind_rows(
abs_HHF_CI %>% 
  filter(time==6|time==12|time==24) %>% 
  mutate(outcome='HHF'),
abs_AC_CI %>% 
  filter(time==6|time==12|time==24) %>% 
  mutate(outcome='Allcause'),
abs_CVdeath_CI %>% 
  filter(time==6|time==12|time==24) %>% 
  mutate(outcome='CVdeath'),
abs_kidney_CI %>% 
  filter(time==6|time==12|time==24) %>% 
  mutate(outcome='kidney')) %>% 
  write_csv('onestage_results.csv')



###8 change-point cox model ----------------------------------------------


fun_change_point <- function(timei){
survSplit(Surv(time, status) ~ treat, data = dt_kidney, cut = c(timei), 
          episode = 'tg', id='...1') %>% 
  coxph(Surv(time, status) ~ treat:strata(tg), data = .) %>% 
  broom::tidy(conf.int=T, exp=T) %>% 
    mutate(time=timei) %>% 
    return()
}
map_dfr(seq(6, 24, 3), fun_change_point) %>% 
  mutate(conf.low=round(conf.low, 3), conf.high=round(conf.high, 3)) %>% 
  mutate(ci=paste0('(', conf.low, ' to ', conf.high, ')')) %>% 
  write_csv('change_points_kidney.csv')














