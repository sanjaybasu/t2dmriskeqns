
rm(list=ls())
library(cvAUC)
library(Hmisc)
library(matrixStats)
load("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_dm_models.RData")
attach(accord_sets)
black = as.numeric(raceclass=='Black')
hisp =  as.numeric(raceclass=='Hispanic')
hdl = as.numeric(hdl)
tob = as.numeric(x4smoke==1)
bmi = wt_kg/(ht_cm/100)^2
bprx = as.numeric((diuretic==1)|(a2rb==1)|(acei==1)|(dhp_ccb==1)|(nondhp_ccb==1)|(alpha_blocker==1)|(central_agent==1)|(beta_blocker==1)|(vasodilator==1)|(reserpine==1)|(other_bpmed==1))
oraldmrx = as.numeric((sulfonylurea==1)|(biguanide==1)|(meglitinide==1)|(ag_inhibitor==1)|(tzd==1)|(other_diabmed==1))

intensivegly = (arm==3)|(arm==4)|(arm==7)|(arm==8)
intensivebp = (arm==1)|(arm==3)
fibratearm = (arm==7)|(arm==5)


# UKPDS nephro: renal failure
amphist=rep(0,length(baseline_age))
blindhist=rep(0,length(baseline_age))
african =black
agedx = baseline_age
sex = female
egfr =gfr
hgb = rep(12,length(baseline_age))
mmalb = as.numeric(ualb>50)
wbc = rep(7.75,length(baseline_age))
timeleft = rep(10,length(baseline_age))
ldluk = accord_sets$ldl * 0.02586
hdluk = accord_sets$hdl * 0.02586
ihd=cvd_hx_baseline
chf=rep(0,length(baseline_age))
hr=accord_sets$hr
a1c=hba1c



retin4 = as.numeric(accord_sets$Retin4==1)
t_censor = rowMaxs(cbind(accord_sets$Retin4Days))
t_retin4s = rowMaxs(cbind(accord_sets$Retin4Days*as.numeric(accord_sets$Retin4)))
t_retin4s[is.na(t_retin4s)]=0
t_retin4s[t_retin4s==0] = t_censor[t_retin4s==0]
t_retin4s[t_retin4s==0] = 'NA'
t_retin4s = as.numeric(t_retin4s)


neph3 = as.numeric(accord_sets$Neph3==1)
t_censor = rowMaxs(cbind(accord_sets$Neph3Days))
t_neph3s = rowMaxs(cbind(accord_sets$Neph3Days*as.numeric(accord_sets$Neph3)))
t_neph3s[is.na(t_neph3s)]=0
t_neph3s[t_neph3s==0] = t_censor[t_neph3s==0]
t_neph3s[t_neph3s==0] = 'NA'
t_neph3s = as.numeric(t_neph3s)



neuro4 = as.numeric(accord_sets$Neuro4==1)
t_censor = rowMaxs(cbind(accord_sets$Neuro4Days))
t_neuro4s = rowMaxs(cbind(accord_sets$Neuro4Days*as.numeric(accord_sets$Neuro4)))
t_neuro4s[is.na(t_neuro4s)]=0
t_neuro4s[t_neuro4s==0] = t_censor[t_neuro4s==0]
t_neuro4s[t_neuro4s==0] = 'NA'
t_neuro4s = as.numeric(t_neuro4s)

cvd = (accord_sets$censor_nmi==0)|(accord_sets$censor_nst==0)|(accord_sets$censor_cm==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25,accord_sets$fuyrs_nst*365.25,accord_sets$fuyrs_cm*365.25))
t_cvds = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25*(1-accord_sets$censor_nmi),accord_sets$fuyrs_nst*365.25*(1-accord_sets$censor_nst),accord_sets$fuyrs_cm*365.25*(1-accord_sets$censor_cm)))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
ascvd = cvd
t_ascvds=t_cvds

fatalmi = (censor_cm==0)&(censor_tst==1)
mi = (accord_sets$censor_nmi==0)|(fatalmi==1)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25,accord_sets$fuyrs_cm*365.25))
t_mis = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25*(1-accord_sets$censor_nmi),accord_sets$fuyrs_cm*365.25*(1-fatalmi)))
t_mis[t_mis==0] = t_censor[t_mis==0]
t_mis[t_mis==0] = 'NA'
t_mis = as.numeric(t_mis)

str = (accord_sets$censor_tst==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_tst*365.25))
t_strs = rowMaxs(cbind(accord_sets$fuyrs_tst*365.25*(1-accord_sets$censor_tst)))
t_strs[t_strs==0] = t_censor[t_strs==0]
t_strs[t_strs==0] = 'NA'
t_strs = as.numeric(t_strs)

chf = (accord_sets$censor_chf==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_chf*365.25))
t_chfs = rowMaxs(cbind(accord_sets$fuyrs_chf*365.25*(1-accord_sets$censor_chf)))
t_chfs[t_chfs==0] = t_censor[t_chfs==0]
t_chfs[t_chfs==0] = 'NA'
t_chfs = as.numeric(t_chfs)


neph = neph3
eye = retin4
t_nephs = t_neph3s
t_eyes = t_retin4s

dp<-data.frame(neuro4, t_neuro4s, retin4, t_retin4s, neph3, t_neph3s, neph,t_nephs,eye,t_eyes,mi,t_mis,str,t_strs,ascvd,t_ascvds,chf,t_chfs,cvd,t_cvds,
               intensivegly,intensivebp,fibratearm,
               baseline_age,female,black,hisp,tob,bmi,
               sbp,dbp,hr,gfr,
               bprx,oraldmrx,anti_coag,statin,
               cvd_hx_baseline,
               hba1c,chol,hdl,hdluk,screat,uacr,
               amphist,blindhist,african,agedx,sex,egfr,hgb,mmalb,wbc,timeleft,
               bmi,ldl,ldluk,sbp,timeleft,
               agedx,a1c,hr,sbp,wbc,chf,ihd)
dp=dp[complete.cases(dp),]
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_nephs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_nephs), 0, dp$neph)


# ukpds neph

q5=(3.549+0.686*dp$african+-0.029*dp$agedx+-0.869*dp$sex+-0.054*dp$bmi+
      -0.487*dp$egfr/10+-0.268*dp$hgb+
      0.027*dp$ldluk*10+1.373*dp$mmalb+(0.085/10)*dp$sbp+0.029*dp$wbc+
      1.108*dp$amphist+0.732*dp$blindhist)
riskrenal = 1-exp(exp(q5)*0-exp(q5)*10)
riskrenal[riskrenal<0]=0
riskrenal[riskrenal>1]=1
estinc_e=riskrenal
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$neph)

betax=(survcox_neph3$coefficients[1]*dp$baseline_age+
         survcox_neph3$coefficients[2]*dp$female+
         survcox_neph3$coefficients[3]*dp$black+
         survcox_neph3$coefficients[4]*dp$hisp+
         survcox_neph3$coefficients[5]*dp$tob+
         survcox_neph3$coefficients[9]*dp$sbp+
         survcox_neph3$coefficients[10]*dp$bprx+
         survcox_neph3$coefficients[11]*dp$oraldmrx+
         survcox_neph3$coefficients[12]*dp$anti_coag+
         survcox_neph3$coefficients[13]*dp$cvd_hx_baseline+
         survcox_neph3$coefficients[14]*dp$hba1c+
         survcox_neph3$coefficients[15]*dp$chol+
         survcox_neph3$coefficients[16]*dp$hdl+
         survcox_neph3$coefficients[17]*dp$screat+
         survcox_neph3$coefficients[18]*dp$uacr)
riskrenalnew = 1 - .96^exp(betax-mean(na.omit(betax)))
estinc_e=riskrenalnew
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$neph)




# ukpds retinopathy

adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_eyes, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_eyes), 0, dp$eye)

q3=(-11.607+0.047*dp$agedx+0.171*dp$a1c+(0.080/10)*dp$hr+(0.068/10)*dp$sbp+0.052*dp$wbc+
      0.841*dp$chf+0.610*dp$ihd) 
riskbld = 1-exp(exp(q3)*0-exp(q3)*10)
riskbld[riskbld<0]=0
riskbld[riskbld>1]=1
estinc_e=riskbld
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$eye)

betax=(survcox_retin4$coefficients[1]*dp$baseline_age+
         survcox_retin4$coefficients[2]*dp$female+
         survcox_retin4$coefficients[3]*dp$black+
         survcox_retin4$coefficients[4]*dp$sbp+
         survcox_retin4$coefficients[5]*dp$bprx+
         survcox_retin4$coefficients[6]*dp$oraldmrx+
         survcox_retin4$coefficients[7]*dp$cvd_hx_baseline+
         survcox_retin4$coefficients[8]*dp$hba1c+
         survcox_retin4$coefficients[9]*dp$chol+
         survcox_retin4$coefficients[10]*dp$hdl+
         survcox_retin4$coefficients[11]*dp$screat+
         survcox_retin4$coefficients[12]*dp$uacr)
riskbldnew = 1 - .93^exp(betax-mean(na.omit(betax)))
estinc_e=riskbldnew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$eye)

# neuropathy

adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_neuro4s, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_neuro4s), 0, dp$neuro4)


betax=(survcox_neuro4$coefficients[1]*dp$baseline_age+
         survcox_neuro4$coefficients[2]*dp$female+
         survcox_neuro4$coefficients[3]*dp$black+
         survcox_neuro4$coefficients[4]*dp$sbp+
         survcox_neuro4$coefficients[5]*dp$bprx+
         survcox_neuro4$coefficients[6]*dp$oraldmrx+
         survcox_neuro4$coefficients[7]*dp$cvd_hx_baseline+
         survcox_neuro4$coefficients[8]*dp$hba1c+
         survcox_neuro4$coefficients[9]*dp$chol+
         survcox_neuro4$coefficients[10]*dp$hdl+
         survcox_neuro4$coefficients[11]*dp$screat)
riskneuronew = 1 - .85^exp(betax-mean(na.omit(betax)))
estinc_e=riskneuronew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$neuro4)


# MI fatal/nonfat
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_mis, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_mis), 0, dp$mi)
riskmi = 1- exp(0-exp(-8.791+
                        0.045*dp$agedx+
                        0.108*dp$hba1c+
                        -0.049*dp$hdluk*10+
                        0.023*dp$ldluk*10+
                        0.203*dp$mmalb+
                        0.046*dp$sbp/10+
                        0.277*dp$tob+
                        0.026*7.75) * 10)
estinc_e=riskmi
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$mi)

betax=(survcox_mi$coefficients[1]*dp$baseline_age+
         survcox_mi$coefficients[2]*dp$female+
         survcox_mi$coefficients[3]*dp$black+
         survcox_mi$coefficients[4]*dp$tob+
         survcox_mi$coefficients[5]*dp$sbp+
         survcox_mi$coefficients[6]*dp$bprx+
         survcox_mi$coefficients[7]*dp$statin+
         survcox_mi$coefficients[8]*dp$anti_coag+
         survcox_mi$coefficients[9]*dp$cvd_hx_baseline+
         survcox_mi$coefficients[10]*dp$hba1c+
         survcox_mi$coefficients[11]*dp$chol+
         survcox_mi$coefficients[12]*dp$hdl+
         survcox_mi$coefficients[13]*dp$screat+
         survcox_mi$coefficients[14]*dp$uacr)
riskminew = 1 - .81^exp(betax-mean(na.omit(betax)))
estinc_e=riskminew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$mi)


# Stroke fatal/nonfat
dp$fu.time <- pmin(dp$t_strs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_strs), 0, dp$str)
riskstr = 1- exp(0-exp(-13.053 + 
                         .066*dp$agedx +
                         -.420*dp$female+
                         -.190*(dp$gfr/10)*(dp$gfr<60)+
                         0.092*dp$hba1c+
                         0.016*dp$ldluk*10+
                         0.420*dp$mmalb+
                         0.170*dp$sbp/10+
                         0.331*dp$tob+
                         0.040*7.75) *  (10)^1.466)
estinc_e=riskstr
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$str)

betax=(survcox_str$coefficients[1]*dp$baseline_age+
         survcox_str$coefficients[2]*dp$female+
         survcox_str$coefficients[3]*dp$black+
         survcox_str$coefficients[4]*dp$tob+
         survcox_str$coefficients[5]*dp$sbp+
         survcox_str$coefficients[6]*dp$bprx+
         survcox_str$coefficients[7]*dp$statin+
         survcox_str$coefficients[8]*dp$anti_coag+
         survcox_str$coefficients[9]*dp$cvd_hx_baseline+
         survcox_str$coefficients[10]*dp$hba1c+
         survcox_str$coefficients[11]*dp$chol+
         survcox_str$coefficients[12]*dp$hdl+
         survcox_str$coefficients[13]*dp$screat+
         survcox_str$coefficients[14]*dp$uacr)
riskstrnew = 1 - .974^exp(betax-mean(na.omit(betax)))
estinc_e=riskstrnew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=7))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$str)




# hard ASCVD
dp$fu.time <- pmin(dp$t_ascvds, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_ascvds), 0, dp$ascvd)


estinc_e=riskmi+riskstr
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$ascvd)

betax=(survcox_ascvd$coefficients[1]*dp$baseline_age+
         survcox_ascvd$coefficients[2]*dp$female+
         survcox_ascvd$coefficients[3]*dp$black+
         survcox_ascvd$coefficients[4]*dp$tob+
         survcox_ascvd$coefficients[5]*dp$sbp+
         survcox_ascvd$coefficients[6]*dp$bprx+
         survcox_ascvd$coefficients[7]*dp$statin+
         survcox_ascvd$coefficients[8]*dp$anti_coag+
         survcox_ascvd$coefficients[9]*dp$cvd_hx_baseline+
         survcox_ascvd$coefficients[10]*dp$hba1c+
         survcox_ascvd$coefficients[11]*dp$chol+
         survcox_ascvd$coefficients[12]*dp$hdl+
         survcox_ascvd$coefficients[13]*dp$screat+
         survcox_ascvd$coefficients[14]*dp$uacr)
riskascvdnew = 1 - .91^exp(betax-mean(na.omit(betax)))
estinc_e=riskascvdnew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$ascvd)

# CHF
dp$fu.time <- pmin(dp$t_chfs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_chfs), 0, dp$chf)
riskchf = 1- exp(0-exp(-12.332 + 
                         0.068*dp$agedx +
                         0.072*dp$bmi+
                         -0.220*(dp$gfr/10)*(dp$gfr<60)+
                         0.012*dp$ldluk*10+
                         0.771*dp$mmalb) * ( 10)^1.514)

estinc_e=riskchf
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$chf)

betax=(survcox_chf$coefficients[1]*dp$baseline_age+
         survcox_chf$coefficients[2]*dp$female+
         survcox_chf$coefficients[3]*dp$black+
         survcox_chf$coefficients[4]*dp$tob+
         survcox_chf$coefficients[5]*dp$sbp+
         survcox_chf$coefficients[6]*dp$bprx+
         survcox_chf$coefficients[7]*dp$statin+
         survcox_chf$coefficients[8]*dp$anti_coag+
         survcox_chf$coefficients[9]*dp$cvd_hx_baseline+
         survcox_chf$coefficients[10]*dp$hba1c+
         survcox_chf$coefficients[11]*dp$chol+
         survcox_chf$coefficients[12]*dp$hdl+
         survcox_chf$coefficients[13]*dp$screat+
         survcox_chf$coefficients[14]*dp$uacr)
riskchfnew = 1 - .963^exp(betax-mean(na.omit(betax)))
estinc_e=riskchfnew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=7))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$chf)




hiukpdsrenalrisk = riskrenal>.1
hinewrenalrisk = riskrenalnew>.1

hiukpdsbldrisk = riskbld>.1
hinewbldrisk = riskbldnew>.1

table(hiukpdsrenalrisk,hinewrenalrisk,dp$neph)
table(hiukpdsbldrisk,hinewbldrisk,dp$eye)



hiukpdsascvdrisk= (riskmi+riskstr)>.1
hiukpdsmirisk = riskmi>.05
hiukpdsstrrisk = riskstr>.05
hiukpdschfrisk = riskchf>.1

hinewascvdrisk= riskascvdnew>.1
hinewmirisk = riskminew>.05
hinewstrrisk = riskstrnew>.05
hinewchfrisk = riskchfnew>.1

table(hiukpdsascvdrisk,hinewascvdrisk,dp$ascvd)
table(hiukpdsmirisk,hinewmirisk,dp$mi)
table(hiukpdsstrrisk,hinewstrrisk,dp$str)
table(hiukpdschfrisk,hinewchfrisk,dp$chf)

ASCVD = riskascvdnew
CHF = riskchfnew
Retinopathy = riskbldnew
Neuropathy = riskneuronew
Nephropathy = riskrenalnew

pairs(~ASCVD+CHF+Retinopathy+Neuropathy+Nephropathy,xlim=c(0,1),ylim=c(0,1))



dp$fu.time <- pmin(dp$t_ascvds, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_ascvds), 0, dp$ascvd)

xbaic = (-29.799*log(dp$baseline_age)+
           4.884*log(dp$baseline_age)^2+
           13.540*log(dp$chol)-
           3.114*log(dp$baseline_age)*log(dp$chol)-
           13.578*log(dp$hdl)+
           3.149*log(dp$baseline_age)*log(dp$hdl)+
           2.019*log(dp$sbp)*dp$bprx+
           1.957*log(dp$sbp)*(1-dp$bprx)+
           7.574*dp$tob-
           1.665*log(dp$baseline_age)*dp$tob+
           0.661*1)*(dp$female==1)*(dp$black==0)+
  (12.344*log(dp$baseline_age)+
     11.853*log(dp$chol)-
     2.664*log(dp$baseline_age)*log(dp$chol)-
     7.990*log(dp$hdl)+
     1.769*log(dp$baseline_age)*log(dp$hdl)+
     1.797*log(dp$sbp)*dp$bprx+
     1.764*(1-dp$bprx)*log(dp$sbp)+
     7.837*dp$tob-
     1.795*log(dp$baseline_age)*dp$tob+
     0.658*1)*(dp$female==1)*(dp$black==1)+
  (17.114*log(dp$baseline_age)+
     .940*log(dp$chol)-
     18.920*log(dp$hdl)+
     4.475*log(dp$baseline_age)*log(dp$hdl)+
     29.291*log(dp$sbp)*dp$bprx-
     6.432*log(dp$baseline_age)*log(dp$sbp)*dp$bprx+
     27.820*log(dp$sbp)*(1-dp$bprx)-
     6.087*(1-dp$bprx)*log(dp$baseline_age)*log(dp$sbp)+
     0.691*dp$tob+0.874*1)*(dp$female==0)*(dp$black==0)+
  (2.469*log(dp$baseline_age)+
     0.302*log(dp$chol)-
     0.307*log(dp$hdl)+
     1.916*dp$bprx*log(dp$sbp)+
     1.809*(1-dp$bprx)*log(dp$sbp)+
     0.549*dp$tob+
     0.645*1)*(dp$female==0)*(dp$black==1)
yhat = xbaic 
meanyhat = -29.180*(dp$female==1)*(dp$black==0)+
  61.180*(dp$female==1)*(dp$black==1)+
  86.610*(dp$female==0)*(dp$black==0)+
  19.540*(dp$female==0)*(dp$black==1)
riskslope = 0.967*(dp$female==1)*(dp$black==0)+
  0.953*(dp$female==1)*(dp$black==1)+
  0.914*(dp$female==0)*(dp$black==0)+
  0.895*(dp$female==0)*(dp$black==1)
riskascvdacc = 1-riskslope^exp((yhat)-(meanyhat))
estinc_e=riskascvdacc
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$ascvd)

hiaccascvdrisk= riskascvdacc>.1
table(hiaccascvdrisk,hinewascvdrisk,dp$ascvd)







# #### Look Ahead data ####
library(matrixStats)
library(Hmisc)
library(cvAUC)
rm(list=ls())
load("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_dm_models.RData")
load("~/Data/Look_AHEAD_V4/Data/End_of_Intervention/Key_Data/lookahead_sbasu_cut.RData")
detach(accord_sets)
detach(dppos_sets)
attach(lookahead_sets)
baseline_age = lookahead_sets$age
female = (lookahead_sets$FEMALE=='Yes')
black = as.numeric(lookahead_sets$RACEVAR)==1
hisp = as.numeric(lookahead_sets$RACEVAR)==2
chol = lookahead_sets$CholmgDL
vldl = 'NA'
ldl = lookahead_sets$ldlchlmgDL
trig = 'NA'
hdl = lookahead_sets$hdlchlmgDL
dbp = lookahead_sets$dbp
sbp = lookahead_sets$sbp
oraldmrx = (lookahead_sets$DIABDRUG=="Yes")&(lookahead_sets$INSULINS=="No")
cvd_hx_baseline = (lookahead_sets$CVDHIS=="Yes")
hba1c = lookahead_sets$hba1cpct
fpg = lookahead_sets$glucosemgDL
screat = lookahead_sets$screatmgDL
hr = lookahead_sets$pulse
bprx = (lookahead_sets$HTNDRUG=="Yes")
insulinrx = (lookahead_sets$INSULINS=="Yes")
statin = (lookahead_sets$STATINS=="Yes")
fibrate = (lookahead_sets$FIBRATES=="Yes")
anti_coag = rep(0,length(baseline_age))
anti_inflam = rep(0,length(baseline_age))
platelet_agi = rep(0,length(baseline_age))
aspirin = (lookahead_sets$ASPIRIN=='Every day')
cpk = 'NA'
mincr = screat/(0.7*female+0.9*(1-female))
mincr[mincr>1] = 1
maxcr = screat/(0.7*female+0.9*(1-female))
maxcr[maxcr<1] = 1
gfr =  141*mincr^(-0.329*female+-0.411*(1-female))*(maxcr^-1.209)*(0.993^baseline_age)*(1.018*female+1*(1-female))*(1.159*black+1*(1-black))
ucreat = lookahead_sets$ucreatmgDL
ualb = lookahead_sets$ualbmgDL
uacr = ualb/ucreat*1000
alt = rep(0,length(baseline_age))
potassium.y = rep(0,length(baseline_age))
tob = as.numeric(lookahead_sets$SMOKING=="Present")
hdl = lookahead_sets$hdlchlmgDL 
ldl = lookahead_sets$ldlchlmgDL 
hdluk = lookahead_sets$hdlchlmgDL * 0.02586
ldluk = lookahead_sets$ldlchlmgDL * 0.02586
mmalb = as.numeric(ualb>50)
agedx = baseline_age
bmi = lookahead_sets$bmi

intensivegly=(rep(FALSE,length(baseline_age)))
intensivebp=(rep(FALSE,length(baseline_age)))
fibratearm=(rep(FALSE,length(baseline_age)))


mi = (lookahead_sets$AllMI==1)
mi[is.na(mi)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_AllMI,lookahead_sets$max))
t_mis = rowMaxs(cbind(lookahead_sets$t_AllMI*(lookahead_sets$AllMI)))
t_mis[is.na(t_mis)]=0
t_mis[t_mis==0] = t_censor[t_mis==0]
t_mis[t_mis==0] = 'NA'
t_mis = as.numeric(t_mis)

str = (lookahead_sets$AllStroke==1)
str[is.na(str)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_AllStroke,lookahead_sets$max))
t_strs = rowMaxs(cbind(lookahead_sets$t_AllStroke*(lookahead_sets$AllStroke)))
t_strs[is.na(t_strs)]=0
t_strs[t_strs==0] = t_censor[t_strs==0]
t_strs[t_strs==0] = 'NA'
t_strs = as.numeric(t_strs)

ascvd = (mi==1) | (str==1)
ascvd[is.na(ascvd)]=0
t_ascvds = rowMins(cbind(t_mis,t_strs))

cvd=ascvd
t_cvds=t_ascvds

chf = (lookahead_sets$CHF==1)
chf[is.na(chf)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_CHF,lookahead_sets$max))
t_chfs = rowMaxs(cbind(lookahead_sets$t_CHF*(lookahead_sets$CHF)))
t_chfs[is.na(t_chfs)]=0
t_chfs[t_chfs==0] = t_censor[t_chfs==0]
t_chfs[t_chfs==0] = 'NA'
t_chfs = as.numeric(t_chfs)






dp<-data.frame(mi,t_mis,str,t_strs,ascvd,t_ascvds,chf,t_chfs,cvd,t_cvds,
               agedx,hdluk, ldluk, mmalb, tob, bmi,
               baseline_age,female,black,hisp,
               sbp,dbp,hr,
               bprx,oraldmrx,statin,anti_coag,
               cvd_hx_baseline,
               hba1c,chol,trig,ldl,hdl,screat,gfr,uacr)
dp=dp[complete.cases(dp),]


# MI fatal/nonfat
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_mis, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_mis), 0, dp$mi)
riskmi = 1- exp(0-exp(-8.791+
                        0.045*dp$agedx+
                        0.108*dp$hba1c+
                      -0.049*dp$hdluk*10+
                        0.023*dp$ldluk*10+
                        0.203*dp$mmalb+
                        0.046*dp$sbp/10+
                        0.277*dp$tob+
                        0.026*7.75) * 10)
estinc_e=riskmi
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=9))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$mi)

betax=(survcox_mi$coefficients[1]*dp$baseline_age+
         survcox_mi$coefficients[2]*dp$female+
         survcox_mi$coefficients[3]*dp$black+
         survcox_mi$coefficients[4]*dp$tob+
         survcox_mi$coefficients[5]*dp$sbp+
         survcox_mi$coefficients[6]*dp$bprx+
         survcox_mi$coefficients[7]*dp$statin+
         survcox_mi$coefficients[8]*dp$anti_coag+
         survcox_mi$coefficients[9]*dp$cvd_hx_baseline+
         survcox_mi$coefficients[10]*dp$hba1c+
         survcox_mi$coefficients[11]*dp$chol+
         survcox_mi$coefficients[12]*dp$hdl+
         survcox_mi$coefficients[13]*dp$screat+
         survcox_mi$coefficients[14]*dp$uacr)
riskminew = 1 - .93^exp(betax-mean(na.omit(betax)))
estinc_e=riskminew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$mi)


# Stroke fatal/nonfat
dp$fu.time <- pmin(dp$t_strs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_strs), 0, dp$str)
riskstr = 1- exp(0-exp(-13.053 + 
                         .066*dp$agedx +
                         -.420*dp$female+
                         -.190*(dp$gfr/10)*(dp$gfr<60)+
                         0.092*dp$hba1c+
                         0.016*dp$ldluk*10+
                         0.420*dp$mmalb+
                         0.170*dp$sbp/10+
                         0.331*dp$tob+
                         0.040*7.75) *  (10)^1.466)
estinc_e=riskstr
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$str)

betax=(survcox_str$coefficients[1]*dp$baseline_age+
         survcox_str$coefficients[2]*dp$female+
         survcox_str$coefficients[3]*dp$black+
         survcox_str$coefficients[4]*dp$tob+
         survcox_str$coefficients[5]*dp$sbp+
         survcox_str$coefficients[6]*dp$bprx+
         survcox_str$coefficients[7]*dp$statin+
         survcox_str$coefficients[8]*dp$anti_coag+
         survcox_str$coefficients[9]*dp$cvd_hx_baseline+
         survcox_str$coefficients[10]*dp$hba1c+
         survcox_str$coefficients[11]*dp$chol+
         survcox_str$coefficients[12]*dp$hdl+
         survcox_str$coefficients[13]*dp$screat+
         survcox_str$coefficients[14]*dp$uacr)
riskstrnew = 1 - .974^exp(betax-mean(na.omit(betax)))
estinc_e=riskstrnew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=7))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$str)


# hard ASCVD
dp$fu.time <- pmin(dp$t_ascvds, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_ascvds), 0, dp$ascvd)


estinc_e=riskmi+riskstr
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$ascvd)

betax=(survcox_ascvd$coefficients[1]*dp$baseline_age+
         survcox_ascvd$coefficients[2]*dp$female+
         survcox_ascvd$coefficients[3]*dp$black+
         survcox_ascvd$coefficients[4]*dp$tob+
         survcox_ascvd$coefficients[5]*dp$sbp+
         survcox_ascvd$coefficients[6]*dp$bprx+
         survcox_ascvd$coefficients[7]*dp$statin+
         survcox_ascvd$coefficients[8]*dp$anti_coag+
         survcox_ascvd$coefficients[9]*dp$cvd_hx_baseline+
         survcox_ascvd$coefficients[10]*dp$hba1c+
         survcox_ascvd$coefficients[11]*dp$chol+
         survcox_ascvd$coefficients[12]*dp$hdl+
         survcox_ascvd$coefficients[13]*dp$screat+
         survcox_ascvd$coefficients[14]*dp$uacr)
riskascvdnew = 1 - .91^exp(betax-mean(na.omit(betax)))
estinc_e=riskascvdnew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$ascvd)

# CHF
dp$fu.time <- pmin(dp$t_chfs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_chfs), 0, dp$chf)
riskchf = 1- exp(0-exp(-12.332 + 
                         0.068*dp$agedx +
                         0.072*dp$bmi+
                         -0.220*(dp$gfr/10)*(dp$gfr<60)+
                         0.012*dp$ldluk*10+
                         0.771*dp$mmalb) * ( 10)^1.514)

estinc_e=riskchf
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$chf)

betax=(survcox_chf$coefficients[1]*dp$baseline_age+
         survcox_chf$coefficients[2]*dp$female+
         survcox_chf$coefficients[3]*dp$black+
         survcox_chf$coefficients[4]*dp$tob+
         survcox_chf$coefficients[5]*dp$sbp+
         survcox_chf$coefficients[6]*dp$bprx+
         survcox_chf$coefficients[7]*dp$statin+
         survcox_chf$coefficients[8]*dp$anti_coag+
         survcox_chf$coefficients[9]*dp$cvd_hx_baseline+
         survcox_chf$coefficients[10]*dp$hba1c+
         survcox_chf$coefficients[11]*dp$chol+
         survcox_chf$coefficients[12]*dp$hdl+
         survcox_chf$coefficients[13]*dp$screat+
         survcox_chf$coefficients[14]*dp$uacr)
riskchfnew = 1 - .963^exp(betax-mean(na.omit(betax)))
estinc_e=riskchfnew
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=7))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$chf)

hiukpdsascvdrisk= (riskmi+riskstr)>.1
hiukpdsmirisk = riskmi>.1
hiukpdsstrrisk = riskstr>.1
hiukpdschfrisk = riskchf>.1

hinewascvdrisk= riskascvdnew>.1
hinewmirisk = riskminew>.1
hinewstrrisk = riskstrnew>.1
hinewchfrisk = riskchfnew>.1

table(hiukpdsascvdrisk,hinewascvdrisk,dp$ascvd)
table(hiukpdsmirisk,hinewmirisk,dp$mi)
table(hiukpdsstrrisk,hinewstrrisk,dp$str)
table(hiukpdschfrisk,hinewchfrisk,dp$chf)





# ACC AHA PCEs for ASCVD


dp$fu.time <- pmin(dp$t_ascvds, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_ascvds), 0, dp$ascvd)

xbaic = (-29.799*log(dp$baseline_age)+
           4.884*log(dp$baseline_age)^2+
           13.540*log(dp$chol)-
           3.114*log(dp$baseline_age)*log(dp$chol)-
           13.578*log(dp$hdl)+
           3.149*log(dp$baseline_age)*log(dp$hdl)+
           2.019*log(dp$sbp)*dp$bprx+
           1.957*log(dp$sbp)*(1-dp$bprx)+
           7.574*dp$tob-
           1.665*log(dp$baseline_age)*dp$tob+
           0.661*1)*(dp$female==1)*(dp$black==0)+
  (12.344*log(dp$baseline_age)+
     11.853*log(dp$chol)-
     2.664*log(dp$baseline_age)*log(dp$chol)-
     7.990*log(dp$hdl)+
     1.769*log(dp$baseline_age)*log(dp$hdl)+
     1.797*log(dp$sbp)*dp$bprx+
     1.764*(1-dp$bprx)*log(dp$sbp)+
     7.837*dp$tob-
     1.795*log(dp$baseline_age)*dp$tob+
     0.658*1)*(dp$female==1)*(dp$black==1)+
  (17.114*log(dp$baseline_age)+
     .940*log(dp$chol)-
     18.920*log(dp$hdl)+
     4.475*log(dp$baseline_age)*log(dp$hdl)+
     29.291*log(dp$sbp)*dp$bprx-
     6.432*log(dp$baseline_age)*log(dp$sbp)*dp$bprx+
     27.820*log(dp$sbp)*(1-dp$bprx)-
     6.087*(1-dp$bprx)*log(dp$baseline_age)*log(dp$sbp)+
     0.691*dp$tob+0.874*1)*(dp$female==0)*(dp$black==0)+
  (2.469*log(dp$baseline_age)+
     0.302*log(dp$chol)-
     0.307*log(dp$hdl)+
     1.916*dp$bprx*log(dp$sbp)+
     1.809*(1-dp$bprx)*log(dp$sbp)+
     0.549*dp$tob+
     0.645*1)*(dp$female==0)*(dp$black==1)
yhat = xbaic 
meanyhat = -29.180*(dp$female==1)*(dp$black==0)+
  61.180*(dp$female==1)*(dp$black==1)+
  86.610*(dp$female==0)*(dp$black==0)+
  19.540*(dp$female==0)*(dp$black==1)
riskslope = 0.967*(dp$female==1)*(dp$black==0)+
  0.953*(dp$female==1)*(dp$black==1)+
  0.914*(dp$female==0)*(dp$black==0)+
  0.895*(dp$female==0)*(dp$black==1)
riskascvdacc = 1-riskslope^exp((yhat)-(meanyhat))
estinc_e=riskascvdacc
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$ascvd)

hiaccascvdrisk= riskascvdacc>.1
table(hiaccascvdrisk,hinewascvdrisk,dp$ascvd)



