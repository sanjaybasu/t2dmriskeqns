# #### DPPOS data ####
# rm(list=ls())
# library(sas7bdat)
# setwd("~/Data/dppos/Data/DPPOS_Phase2/Non_Form_Based")
# dppos_events = read.sas7bdat("events.sas7bdat")
# dppos_demographic = read.sas7bdat("demographic.sas7bdat")
# dppos_microvascular = read.sas7bdat("microvascular.sas7bdat")
# dppos_lab = read.sas7bdat("lab.sas7bdat")
# setwd("~/Data/dppos/Data/DPPOS_Phase1/Non_Form_Based")
# dppos_laboratory = read.sas7bdat("laboratory.sas7bdat")
# setwd("~/Data/dppos/Data/DPPOS_Phase1/Form_Based")
# dppos_r04 = read.sas7bdat("r04.sas7bdat")
# dppos_f03 = read.sas7bdat("f03.sas7bdat")
# dppos_f01 = read.sas7bdat("f01.sas7bdat")
# dppos_f02 = read.sas7bdat("f02.sas7bdat")
# dppos_f04 = read.sas7bdat("f04.sas7bdat")
# setwd("~/Data/dppos/Data/DPPOS_Phase2/Form_Based")
# dppos_r042 = read.sas7bdat("r04.sas7bdat")
# dppos_f032 = read.sas7bdat("f03.sas7bdat")
# dppos_f012 = read.sas7bdat("f01.sas7bdat")
# dppos_f022 = read.sas7bdat("f02.sas7bdat")
# dppos_f042 = read.sas7bdat("f04.sas7bdat")
# setwd("~/Data/dpp/Data/DPP_Data_2008/Form_Data/Data")
# dpp_q08 = read.sas7bdat("q08.sas7bdat")
# 
# save.image("~/Data/dppos/Data/DPPOS_Phase2/Non_Form_Based/dppos_sbasu.RData")
# 
# rm(list=ls())
# load("~/Data/dppos/Data/DPPOS_Phase2/Non_Form_Based/dppos_sbasu.RData")
# dppos_events_cut = dppos_events[which(dppos_events$DIABF==1),]
# dppos_laboratory_cut = dppos_laboratory[(dppos_laboratory$VISIT=="01A"),]
# dppos_sets = merge(dppos_events_cut,dppos_laboratory_cut,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dppos_sets$UCRE = NA
# dppos_sets = dppos_sets[which(dppos_sets$DIABF==1),]
# dppos_lab_cut = dppos_lab[(dppos_lab$VISIT=="11A"),]
# dppos_lab_cut = dppos_lab_cut[c("RELEASE_ID","UCRE")]
# dppos_sets2 = merge(dppos_sets,dppos_lab_cut,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dppos_sets3=dppos_sets2[order(dppos_sets2$RELEASE_ID),]
# dppos_sets3 = dppos_sets3[which(dppos_sets3$DIABF==1),]
# dppos_sets3$yrdiff = abs(dppos_sets3$DAYSRAND/365.25-dppos_sets3$DIABT)
# dppos_sets_minyrdiff = aggregate(yrdiff~RELEASE_ID,min,data=dppos_sets3)
# dppos_sets = merge(dppos_sets3,dppos_sets_minyrdiff,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dppos_sets_cut = dppos_sets[which(dppos_sets$yrdiff.x-dppos_sets$yrdiff.y<=1),]
# library(doBy)
# dppos_sets_cut = summaryBy(. ~RELEASE_ID,data = dppos_sets_cut, na.rm=TRUE)
# dppos_sets = merge(dppos_sets_cut,dppos_demographic,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dppos_sets = merge(dppos_sets,dppos_microvascular,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dppos_f01_cut = summaryBy(QPSBP1 + QPSBP2 + QPDBP1 + QPDBP2 ~RELEASE_ID,data = dppos_f01, na.rm=TRUE)
# # dppos_f02_cut = summaryBy(APSBP1 + APSBP2 + APDBP1 + APDBP2 ~RELEASE_ID,data = dppos_f02, na.rm=TRUE)
# # dppos_f03_cut = summaryBy(JISBP1 + JISBP2 + JIDBP1 + JIDBP2 ~RELEASE_ID,data = dppos_f03, na.rm=TRUE)
# #dppos_f04_cut = summaryBy(JIHYPMG+JISBP1 + JISBP2 + JIDBP1 + JIDBP2 ~RELEASE_ID,data = dppos_f04, na.rm=TRUE)
# # dppos_f012_cut = summaryBy(QPSBP1 + QPSBP2 + QPDBP1 + QPDBP2 ~RELEASE_ID,data = dppos_f012, na.rm=TRUE)
# #dppos_f022_cut = summaryBy(JIHYPMG+JISBP1 + JISBP2 + JIDBP1 + JIDBP2 ~RELEASE_ID,data = dppos_f022, na.rm=TRUE)
# #dppos_f032_cut = summaryBy(JIHYPMG+JISBP1 + JISBP2 + JIDBP1 + JIDBP2 ~RELEASE_ID,data = dppos_f032, na.rm=TRUE)
# #dppos_f042_cut = summaryBy(JIHYPMG+JISBP1 + JISBP2 + JIDBP1 + JIDBP2 ~RELEASE_ID,data = dppos_f042, na.rm=TRUE)
# colnames(dppos_f01_cut)[2:5] = c("SBP1","SBP2","DBP1","DBP2")
# # colnames(dppos_f02_cut)[2:5] = c("SBP1","SBP2","DBP1","DBP2")
# # colnames(dppos_f03_cut)[2:5] = c("SBP1","SBP2","DBP1","DBP2")
# # colnames(dppos_f012_cut)[2:5] = c("SBP1","SBP2","DBP1","DBP2")
# dppos_sets4=rbind(dppos_f01_cut)#,dppos_f02_cut,dppos_f03_cut,dppos_f012_cut)
# dppos_sets4 = summaryBy(SBP1 + SBP2 + DBP1 + DBP2 ~RELEASE_ID,data = dppos_sets4, na.rm=TRUE)
# dppos_sets = merge(dppos_sets,dppos_sets4,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dpp_q08_cut = summaryBy(IHMI+IHSTRK~RELEASE_ID,data = dpp_q08, na.rm=TRUE,FUN=min)
# dppos_sets = merge(dppos_sets,dpp_q08_cut,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dppos_r042_cut = summaryBy(. ~RELEASE_ID,data = dppos_r042, na.rm=TRUE,FUN=max)
# dppos_sets = merge(dppos_sets,dppos_r042_cut,by="RELEASE_ID",all.x=TRUE,all.y=TRUE)
# dppos_sets = dppos_sets[which(dppos_sets$DIABF==1),]
# 
# 
# save.image("~/Data/dppos/Data/DPPOS_Phase2/Non_Form_Based/dppos_sbasu_cut.RData")
# 

#### ext validation: dppos ####
rm(list=ls())
load("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_dm_models.RData")
load("~/Data/dppos/Data/DPPOS_Phase2/Non_Form_Based/dppos_sbasu_cut.RData")
baseline_age = 40*(AGEGROUP==1)+42.5*(AGEGROUP==2)+47.5*(AGEGROUP==3)+52.5*(AGEGROUP==4)+57.5*(AGEGROUP==5)+62.5*(AGEGROUP==6)+65*(AGEGROUP==7)
female = (SEX-1)
black = as.numeric(RACE_ETH==2)
hisp = (RACE_ETH==3)
chol = CHOL.mean
vldl = VLDL.mean
ldl = LDLC.mean
trig = TRIG.mean
hdl = CHDL.mean
dbp = rowMeans(cbind(DBP1.mean,DBP2.mean))
sbp = rowMeans(cbind(SBP1.mean,SBP2.mean))
oraldmrx = as.numeric(ASSIGN=="Metformin")   
cvd_hx_baseline = as.numeric((IHMI.min==1)|(IHSTRK.min==1))
hba1c = HBA1.mean
fpg = G000.mean
screat = CREA.mean
hr = rep(mean(na.omit(accord_sets$hr)),length(baseline_age))
bprx = as.numeric(CHAHMED.max==1)
insulinrx = rep(0,length(baseline_age))
statin = CHDRUG.max==1
fibrate = rep(0,length(baseline_age))
anti_coag = rep(0,length(baseline_age))
anti_inflam = rep(0,length(baseline_age))
platelet_agi = rep(0,length(baseline_age))
aspirin = rep(0,length(baseline_age))
cpk = rep(147.7329,length(baseline_age))
mincr = screat/(0.7*female+0.9*(1-female))
mincr[mincr>1] = 1
maxcr = screat/(0.7*female+0.9*(1-female))
maxcr[maxcr<1] = 1
gfr = 141*mincr^(-0.329*female+-0.411*(1-female))*(maxcr^-1.209)*(0.993^baseline_age)*(1.018*female+1*(1-female))*(1.159*black+1*(1-black))
ucreat = UCRE.y.mean
ualb = 10.91567/124.601*ucreat
uacr = ualb/ucreat*1000
alt = rep(26.27626,length(baseline_age))
potassium.y = rep(4.490411,length(baseline_age))

neph = (evtnep==1)#&(11-DIABT.mean<=5)&(11-DIABT.mean>=0)
eye= (evtret==1)#&(11-DIABT.mean<=5)&(11-DIABT.mean>=0)
neuro= (evtneu==1)#&(11-DIABT.mean<=5)&(11-DIABT.mean>=0)

intensivegly=(rep(FALSE,length(baseline_age)))
intensivebp=(rep(FALSE,length(baseline_age)))
fibratearm=(rep(FALSE,length(baseline_age)))


###### Nephropathy ######
#: micro- or macro-albuminuria (≥30 mg/gram creatinine, confirmed), (ACCORD NEPH 5 or NEPH 2)
# or renal dysfunction (end-stage renal disease, dialysis or renal transplant)  (ACCORD NEPH 3)
# or GFR < 45 ml per min based on serum creatinine, using the CKD-EPI equation or another validated algorithm; the qualifying criteria confirmed) (possibly add ACCORD NEPH 1)
t_censor = (11-DIABT.mean)*365.25
t_censor[t_censor<0]=0
t_dppnephs = (11-DIABT.mean)*365.25*neph
t_dppnephs[t_dppnephs<0]=0
t_dppnephs[is.na(t_dppnephs)]=0
t_dppnephs[t_dppnephs==0] = t_censor[t_dppnephs==0]
t_dppnephs[t_dppnephs==0] = 'NA'
t_dppnephs = as.numeric(t_dppnephs)
t_nephs = t_dppnephs

dp<-data.frame(neph,t_nephs,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=5*365.25
dp$fu.time <- pmin(dp$t_nephs, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_nephs), 0, dp$neph)
betax=(survcox_neph235$coefficients[1]*dp$baseline_age+survcox_neph235$coefficients[2]*dp$female+survcox_neph235$coefficients[3]*dp$black+survcox_neph235$coefficients[4]*dp$sbp+survcox_neph235$coefficients[5]*dp$bprx+survcox_neph235$coefficients[6]*dp$oraldmrx+survcox_neph235$coefficients[7]*dp$cvd_hx_baseline+survcox_neph235$coefficients[8]*dp$hba1c+survcox_neph235$coefficients[9]*dp$chol+survcox_neph235$coefficients[10]*dp$hdl+survcox_neph235$coefficients[11]*dp$screat+survcox_neph235$coefficients[12]*dp$uacr)
risk = 1 - .97^exp(betax-mean(na.omit(betax)))
estinc_e=risk
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=10))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$neph)



##### Retinopathy #####
#: retinopathy by fundus photography (ETDRS grade of 20 or greater) (possibly ACCORD Retin4)
# or adjudicated history of laser or other treatment for retinopathy  (Retin 1)  
t_censor = (11-DIABT.mean)*365.25
t_censor[t_censor<0]=0
t_dppeyes = (11-DIABT.mean)*365.25*eye
t_dppeyes[t_dppeyes<0]=0
t_dppeyes[is.na(t_dppeyes)]=0
t_dppeyes[t_dppeyes==0] = t_censor[t_dppeyes==0]
t_dppeyes[t_dppeyes==0] = 'NA'
t_dppeyes = as.numeric(t_dppeyes)
t_eyes = t_dppeyes

dp<-data.frame(eye,t_eyes,
               baseline_age,female,black,
               sbp,
               bprx,oraldmrx,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=5*365.25
dp$fu.time <- pmin(dp$t_eyes, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_eyes), 0, dp$eye)
betax=(survcox_retin14$coefficients[1]*dp$baseline_age+survcox_retin14$coefficients[2]*dp$female+survcox_retin14$coefficients[3]*dp$black+survcox_retin14$coefficients[4]*dp$sbp+survcox_retin14$coefficients[5]*dp$bprx+survcox_retin14$coefficients[6]*dp$oraldmrx+survcox_retin14$coefficients[7]*dp$cvd_hx_baseline+survcox_retin14$coefficients[8]*dp$hba1c+survcox_retin14$coefficients[9]*dp$chol+survcox_retin14$coefficients[10]*dp$hdl+survcox_retin14$coefficients[11]*dp$screat+survcox_retin14$coefficients[12]*dp$uacr)
risk = 1 - .97^exp(betax-mean(na.omit(betax)))
estinc_e=risk
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=8))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$eye)


##### Neuropathy #####
#: reduction or absence of light touch sensation to monofilament 
# (Semmes- Weinstein 10 gram) in either foot (< 8 of 10 applications detected). (neuro4)
t_censor = (11-DIABT.mean)*365.25
t_censor[t_censor<0]=0
t_dppneuro = (11-DIABT.mean)*365.25*neuro
t_dppneuro[t_dppneuro<0]=0
t_dppneuro[is.na(t_dppneuro)]=0
t_dppneuro[t_dppneuro==0] = t_censor[t_dppneuro==0]
t_dppneuro[t_dppneuro==0] = 'NA'
t_dppneuro = as.numeric(t_dppneuro)
t_neuro = t_dppneuro

dp<-data.frame(neuro,t_neuro,
               baseline_age,female,black,
               sbp,
               bprx,oraldmrx,
               cvd_hx_baseline,
               hba1c,chol,hdl,screat,uacr)
dp=dp[complete.cases(dp),]
adm.cens=5*365.25
dp$fu.time <- pmin(dp$t_neuro, adm.cens)
dp$status <- ifelse(as.numeric(adm.cens < dp$t_neuro), 0, dp$neuro)
betax=(survcox_neuro4$coefficients[1]*dp$baseline_age+survcox_neuro4$coefficients[2]*dp$female+survcox_neuro4$coefficients[3]*dp$black+survcox_neuro4$coefficients[4]*dp$sbp+survcox_neuro4$coefficients[5]*dp$bprx+survcox_neuro4$coefficients[6]*dp$oraldmrx+survcox_neuro4$coefficients[7]*dp$cvd_hx_baseline+survcox_neuro4$coefficients[8]*dp$hba1c+survcox_neuro4$coefficients[9]*dp$chol+survcox_neuro4$coefficients[10]*dp$hdl+survcox_neuro4$coefficients[11]*dp$screat+survcox_neuro4$coefficients[12]*dp$uacr)
risk = 1 - .97^exp(betax-mean(na.omit(betax)))
estinc_e=risk
#estinc_e=1-survfit_e$surv[dim(survfit_e$surv)[1],]
dp$dec=as.numeric(cut2(estinc_e, g=5))
GND.result=GND.calib(pred=estinc_e, tvar=dp$fu.time, out=dp$status, 
                     cens.t=adm.cens, groups=dp$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_e,dp$neuro)

