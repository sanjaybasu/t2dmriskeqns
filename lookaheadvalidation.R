
# #### Look Ahead data ####
# rm(list=ls())
# library(sas7bdat)
# setwd("~/Data/Look_AHEAD_V4/Data/Baseline/Analysis_Data")
# baseline_combined = read.sas7bdat("baseline_combined.sas7bdat")
# setwd("~/Data/Look_AHEAD_V4/Data/End_of_Intervention/Key_Data")
# la4_outcomes1 = read.sas7bdat("la4_outcomes1.sas7bdat")
# la3_activitystatus = read.sas7bdat("la3_activitystatus.sas7bdat")
# save.image("lookahead_sbasu.RData")
#
# rm(list=ls())
# library(plyr)
# load("~/Data/Look_AHEAD_V4/Data/End_of_Intervention/Key_Data/lookahead_sbasu.RData")
# baseline_combined$MaskID = baseline_combined$P_ID
# lookahead_sets = merge(baseline_combined,la4_outcomes1,by="MaskID",all.x=TRUE,all.y=TRUE)
# la3_max = data.frame(MaskID = la3_activitystatus$MaskID, daysAvg = la3_activitystatus$daysAvg, active = la3_activitystatus$active)
# la3_max = la3_max[(la3_max$active==2),]
# la3_max = ddply(la3_max, "MaskID", summarize, max = max(daysAvg))
# lookahead_sets = merge(lookahead_sets,la3_max,by="MaskID",all.x=TRUE,all.y=TRUE)
# save.image("~/Data/Look_AHEAD_V4/Data/End_of_Intervention/Key_Data/lookahead_sbasu_cut.RData")


#### ext validation: lookahead ####
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
bmi=as.numeric(lookahead_sets$bmi)

intensivegly=(rep(FALSE,length(baseline_age)))
intensivebp=(rep(FALSE,length(baseline_age)))
fibratearm=(rep(FALSE,length(baseline_age)))


# sample = data.frame(AllMI,AllStroke,CHF,CVD_death,Death,baseline_age,female,black,hisp,tob,bmi,hr,
#                     sbp,dbp,
#                     bprx,oraldmrx,anti_coag,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
#                     cvd_hx_baseline,
#                     hba1c,chol,hdl,screat,ucreat,ualb,uacr,gfr,fpg)
# sample=sample[complete.cases(sample),]



# CVD mort
cvdmort = (lookahead_sets$CVD_death==1)
cvdmort[is.na(cvdmort)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_CVD_death,lookahead_sets$max))
t_cvdmorts = rowMaxs(cbind(lookahead_sets$t_CVD_death*(lookahead_sets$CVD_death)))
t_cvdmorts[is.na(t_cvdmorts)]=0
t_cvdmorts[t_cvdmorts==0] = t_censor[t_cvdmorts==0]
t_cvdmorts[t_cvdmorts==0] = 'NA'
t_cvdmorts = as.numeric(t_cvdmorts)

dp<-data.frame(cvdmort,t_cvdmorts,
               baseline_age,female,black,tob,
               sbp,
               bprx,statin,anti_coag,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_cvdmorts, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_cvdmorts), 0, dp$cvdmort)
betax=(survcox_cvdmort$coefficients[1]*dp$baseline_age+
         survcox_cvdmort$coefficients[2]*dp$female+
         survcox_cvdmort$coefficients[3]*dp$black+
         survcox_cvdmort$coefficients[4]*dp$tob+
         survcox_cvdmort$coefficients[5]*dp$sbp+
         survcox_cvdmort$coefficients[6]*dp$bprx+
         survcox_cvdmort$coefficients[7]*dp$statin+
         survcox_cvdmort$coefficients[8]*dp$anti_coag+
         survcox_cvdmort$coefficients[9]*dp$cvd_hx_baseline+
         survcox_cvdmort$coefficients[10]*dp$hba1c+
         survcox_cvdmort$coefficients[11]*dp$chol+
         survcox_cvdmort$coefficients[12]*dp$hdl+
         survcox_cvdmort$coefficients[13]*dp$screat+
         survcox_cvdmort$coefficients[14]*dp$uacr)

    risk = 1 - .974^exp(betax-mean(na.omit(betax)))
    estinc_e=risk
    #estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
    dp$dec=as.numeric(cut2(estinc_e, g=4))
    GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                         cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
    GND.result
    ci.cvAUC(estinc_e,dp$cvdmort)



# MI fatal/nonfat
mi = (lookahead_sets$AllMI==1)
mi[is.na(mi)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_AllMI,lookahead_sets$max))
t_mis = rowMaxs(cbind(lookahead_sets$t_AllMI*(lookahead_sets$AllMI)))
t_mis[is.na(t_mis)]=0
t_mis[t_mis==0] = t_censor[t_mis==0]
t_mis[t_mis==0] = 'NA'
t_mis = as.numeric(t_mis)

dp<-data.frame(mi,t_mis,
               baseline_age,female,black,tob,
               sbp,
               bprx,statin,anti_coag,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_mis, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_mis), 0, dp$mi)
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
risk = 1 - .93^exp(betax-mean(na.omit(betax)))
estinc_e=risk
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$mi)


# Stroke fatal/nonfat
str = (lookahead_sets$AllStroke==1)
str[is.na(str)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_AllStroke,lookahead_sets$max))
t_strs = rowMaxs(cbind(lookahead_sets$t_AllStroke*(lookahead_sets$AllStroke)))
t_strs[is.na(t_strs)]=0
t_strs[t_strs==0] = t_censor[t_strs==0]
t_strs[t_strs==0] = 'NA'
t_strs = as.numeric(t_strs)

dp<-data.frame(str,t_strs,
               baseline_age,female,black,tob,
               sbp,
               bprx,statin,anti_coag,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_strs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_strs), 0, dp$str)
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
risk = 1 - .976^exp(betax-mean(na.omit(betax)))
estinc_e=risk
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=7))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$str)


# hard ASCVD

ascvd = (lookahead_sets$AllStroke==1)|(lookahead_sets$AllMI==1)|(lookahead_sets$CVD_death==1)
ascvd[is.na(ascvd)]=0
t_censor= rowMaxs(cbind(lookahead_sets$t_AllMI,lookahead_sets$t_AllStroke,lookahead_sets$t_CVD_death,lookahead_sets$max))
t_ascvds = rowMaxs(cbind(lookahead_sets$t_AllStroke*(lookahead_sets$AllStroke),lookahead_sets$t_AllMI*(lookahead_sets$AllMI),lookahead_sets$t_CVD_death*(lookahead_sets$CVD_death)))
t_ascvds[is.na(t_ascvds)]=0
t_ascvds[t_ascvds==0] = t_censor[t_ascvds==0]
t_ascvds[t_ascvds==0] = 'NA'
t_ascvds = as.numeric(t_ascvds)

dp<-data.frame(ascvd,t_ascvds,intensivegly,intensivebp,fibratearm,
               baseline_age,female,black,tob,hisp,bmi,
               sbp,
               bprx,statin,anti_coag,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=10*365.25 
dp$fu.time <- pmin(dp$t_ascvds, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_ascvds), 0, dp$ascvd)
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

    risk = 1 - .852^exp(betax-mean(na.omit(betax)))
    estinc_e=risk
    #estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
    dp$dec=as.numeric(cut2(estinc_e, g=10))
    GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status,
                         cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
    GND.result
    ci.cvAUC(estinc_e,dp$ascvd)





 
# CHF
chf = (lookahead_sets$CHF==1)
chf[is.na(chf)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_CHF,lookahead_sets$max))
t_chfs = rowMaxs(cbind(lookahead_sets$t_CHF*(lookahead_sets$CHF)))
t_chfs[is.na(t_chfs)]=0
t_chfs[t_chfs==0] = t_censor[t_chfs==0]
t_chfs[t_chfs==0] = 'NA'
t_chfs = as.numeric(t_chfs)

dp<-data.frame(chf,t_chfs,
               baseline_age,female,black,tob,
               sbp,
               bprx,statin,anti_coag,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_chfs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_chfs), 0, dp$chf)
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
risk = 1 - .962^exp(betax-mean(na.omit(betax)))
estinc_e=risk
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=7))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$chf)


# allmort
allmort = (lookahead_sets$Death==1)
allmort[is.na(allmort)]=0
t_censor = rowMaxs(cbind(lookahead_sets$t_Death,lookahead_sets$max))
t_allmorts = rowMaxs(cbind(lookahead_sets$t_Death*(lookahead_sets$Death)))
t_allmorts[is.na(t_allmorts)]=0
t_allmorts[t_allmorts==0] = t_censor[t_allmorts==0]
t_allmorts[t_allmorts==0] = 'NA'
t_allmorts = as.numeric(t_allmorts)

dp<-data.frame(allmort,t_allmorts,
               baseline_age,female,black,tob,
               sbp,
               bprx,statin,anti_coag,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=10*365.25
dp$fu.time <- pmin(dp$t_allmorts, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_allmorts), 0, dp$allmort)
betax=(survcox_allmort$coefficients[1]*dp$baseline_age+
         survcox_allmort$coefficients[2]*dp$female+
         survcox_allmort$coefficients[3]*dp$black+
         survcox_allmort$coefficients[4]*dp$tob+
         survcox_allmort$coefficients[5]*dp$sbp+
         survcox_allmort$coefficients[6]*dp$bprx+
         survcox_allmort$coefficients[7]*dp$statin+
         survcox_allmort$coefficients[8]*dp$anti_coag+
         survcox_allmort$coefficients[9]*dp$cvd_hx_baseline+
         survcox_allmort$coefficients[10]*dp$hba1c+
         survcox_allmort$coefficients[11]*dp$chol+
         survcox_allmort$coefficients[12]*dp$hdl+
         survcox_allmort$coefficients[13]*dp$screat+
         survcox_allmort$coefficients[14]*dp$uacr)
risk = 1 - .933^exp(betax-mean(na.omit(betax)))
estinc_e=risk
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$allmort)



