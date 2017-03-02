# 
# #### accord data ####
# 
# rm(list=ls())
# library(sas7bdat)
# setwd("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets")
# accord_key = read.sas7bdat("accord_key.sas7bdat")
# activitystatus = read.sas7bdat("activitystatus.sas7bdat")
# bloodpressure = read.sas7bdat("bloodpressure.sas7bdat")
# concomitantmeds = read.sas7bdat("concomitantmeds.sas7bdat")
# cvdoutcomes = read.sas7bdat("cvdoutcomes.sas7bdat")
# eye = read.sas7bdat("eye.sas7bdat")
# hba1c = read.sas7bdat("hba1c.sas7bdat")
# hypoglycemiaevents = read.sas7bdat("hypoglycemiaevents.sas7bdat")
# hypoglycemiatime1st = read.sas7bdat("hypoglycemiatime1st.sas7bdat")
# lipids = read.sas7bdat("lipids.sas7bdat")
# microvascularoutcomes = read.sas7bdat("microvascularoutcomes.sas7bdat")
# otherlabs = read.sas7bdat("otherlabs.sas7bdat")
# sae = read.sas7bdat("sae.sas7bdat")
# setwd("~/Data/accord/4-Data_Sets-CRFs/4a-CRF_Data_Sets")
# tobacco = read.sas7bdat("f01_inclusionexclusionsummary.sas7bdat")
# bmi = read.sas7bdat("f07_baselinehistoryphysicalexam.sas7bdat")
# save.image("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_sbasu.RData")
# 
# rm(list=ls())
# load("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_sbasu.RData")
# activitystatus_cut = activitystatus[which(activitystatus$Visit=='BLR'),]
# bloodpressure_cut = bloodpressure[which(bloodpressure$Visit=='BLR'),]
# concomitantmeds_cut = concomitantmeds[which(concomitantmeds$Visit=='BLR'),]
# hba1c_cut = hba1c[which(hba1c$Visit=='BLR'),]
# lipids_cut = lipids[which(lipids$Visit=='BLR'),]
# otherlabs_cut = otherlabs[which(otherlabs$Visit=='BLR'),]
# bmi_cut = bmi[which(bmi$Visit=='BLR'),]
# accord_sets = merge(accord_key,bloodpressure_cut,by="MaskID",all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,activitystatus_cut,by=c("MaskID"),all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,concomitantmeds_cut,by=c("MaskID"),all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,cvdoutcomes,by="MaskID",all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,eye,by="MaskID",all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,hba1c_cut,by=c("MaskID","Visit"),all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,lipids_cut,by=c("MaskID","Visit"),all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,otherlabs_cut,by=c("MaskID","Visit"),all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,microvascularoutcomes,by="MaskID",all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,sae,by="MaskID",all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,tobacco,by="MaskID",all.x=TRUE,all.y=TRUE)
# accord_sets = merge(accord_sets,bmi_cut,by="MaskID",all.x=TRUE,all.y=TRUE)
# save.image("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_sbasu_cut.RData")

##### model fitting ####

rm(list=ls())
load("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_sbasu_cut.RData")
library(gdata)
library(survival)
library(matrixStats)
library(Hmisc)
library(glmnet) # https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html#log
library(cvAUC)
library(doMC)
registerDoMC(cores=4)
keep(accord_sets, sure=TRUE)
attach(accord_sets)
kmdec=function(dec.num,dec.name, datain, adm.cens){
  stopped=0
  data.sub=datain[datain[,dec.name]==dec.num,]
  if (sum(data.sub$out)>1){
    avsurv=survfit(Surv(tvar,out) ~ 1, data=datain[datain[,dec.name]==dec.num,], error="g")
    avsurv.est=ifelse(min(avsurv$time)<=adm.cens,avsurv$surv[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],1)
    avsurv.stderr=ifelse(min(avsurv$time)<=adm.cens,avsurv$std.err[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
    avsurv.stderr=avsurv.stderr*avsurv.est
    avsurv.num=ifelse(min(avsurv$time)<=adm.cens,avsurv$n.risk[avsurv$time==max(avsurv$time[avsurv$time<=adm.cens])],0)
  } else {
    return(c(0,0,0,0,stopped=-1))
  }
  if (sum(data.sub$out)<5) stopped=1
  c(avsurv.est, avsurv.stderr, avsurv.num, dec.num, stopped) 
}
GND.calib = function(pred, tvar, out, cens.t, groups, adm.cens){
  tvar.t=ifelse(tvar>adm.cens, adm.cens, tvar)
  out.t=ifelse(tvar>adm.cens, 0, out)
  datause=data.frame(pred=pred, tvar=tvar.t, out=out.t, count=1, cens.t=cens.t, dec=groups)
  numcat=length(unique(datause$dec))
  groups=sort(unique(datause$dec))
  kmtab=matrix(unlist(lapply(groups,kmdec,"dec",datain=datause, adm.cens)),ncol=5, byrow=TRUE)
  # if (any(kmtab[,5] == -1)) stop("Stopped because at least one of the groups contains <2 events. Consider collapsing some groups.")
  # else if (any(kmtab[,5] == 1)) warning("At least one of the groups contains < 5 events. GND can become unstable.\ 
  #                                       (see Demler, Paynter, Cook 'Tests of Calibration and Goodness of Fit in the Survival Setting' DOI: 10.1002/sim.6428) \
  #                                       Consider collapsing some groups to avoid this problem.")
  hltab=data.frame(group=kmtab[,4],
                   totaln=tapply(datause$count,datause$dec,sum),
                   censn=tapply(datause$cens.t,datause$dec,sum),
                   numevents=tapply(datause$out,datause$dec,sum),
                   expected=tapply(datause$pred,datause$dec,sum),
                   kmperc=1-kmtab[,1], 
                   kmvar=kmtab[,2]^2, 
                   kmnrisk=kmtab[,3],
                   expectedperc=tapply(datause$pred,datause$dec,mean))
  hltab$kmnum=hltab$kmperc*hltab$totaln
  hltab$GND_component=ifelse(hltab$kmvar==0, 0,(hltab$kmperc-hltab$expectedperc)^2/(hltab$kmvar))
  print(hltab[c(1,2,3,4,10,5,6,9,7,11)], digits=4)
  plot(tapply(datause$pred,datause$dec,mean),1-kmtab[,1],xlab="Expected K-M rate",ylab="Observed K-M rate",xlim=c(0,1),ylim=c(0,1))
  abline(a=0,b=1, col = "gray60")
  calline = lm(hltab$kmperc~hltab$expectedperc)
  c(df=numcat-1, chi2gw=sum(hltab$GND_component),pvalgw=1-pchisq(sum(hltab$GND_component),numcat-1),slope=calline$coefficients[2],intercept = calline$coefficients[1])
}
diuretic = as.numeric((loop==1)|(thiazide==1)|(ksparing==1))
bprx = as.numeric((diuretic==1)|(a2rb==1)|(acei==1)|(dhp_ccb==1)|(nondhp_ccb==1)|(alpha_blocker==1)|(central_agent==1)|(beta_blocker==1)|(vasodilator==1)|(reserpine==1)|(other_bpmed==1))
oraldmrx = as.numeric((sulfonylurea==1)|(biguanide==1)|(meglitinide==1)|(ag_inhibitor==1)|(tzd==1)|(other_diabmed==1))
insulinrx = as.numeric((nphl_insulin==1)|(reg_insulin==1)|(la_insulin==1)|(othbol_insulin==1)|(premix_insulin==1))
black = as.numeric(raceclass=='Black')
hisp =  as.numeric(raceclass=='Hispanic')
hdl = as.numeric(hdl)
tob = as.numeric(x4smoke==1)
bmi = wt_kg/(ht_cm/100)^2

intensivegly = (arm==3)|(arm==4)|(arm==7)|(arm==8)
intensivebp = (arm==1)|(arm==3)
fibratearm = (arm==7)|(arm==5)

sample = data.frame(intensivegly,intensivebp,fibratearm,
                    baseline_age,female,black,hisp,tob,bmi,
                    sbp,dbp,hr,
                    bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                    cvd_hx_baseline,
                    hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,gfr,ualb,ucreat,uacr)
sample=sample[complete.cases(sample),]


##### hard ASCVD: MI or stroke, fatal or nonfatal  #####
cvd = (accord_sets$censor_nmi==0)|(accord_sets$censor_nst==0)|(accord_sets$censor_cm==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25,accord_sets$fuyrs_nst*365.25,accord_sets$fuyrs_cm*365.25))
t_cvds = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25*(1-accord_sets$censor_nmi),accord_sets$fuyrs_nst*365.25*(1-accord_sets$censor_nst),accord_sets$fuyrs_cm*365.25*(1-accord_sets$censor_cm)))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
cOutcome = Surv(time=t_cvds, event = cvd)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,gfr,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
cvdmodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(cvdmodel.cv)
coef.cv = coef(cvdmodel.cv, s = 'lambda.min')
coef.cv
c<-data.frame(cvd,t_cvds,
              baseline_age,female,black,tob,
              sbp,
              bprx,statin,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
c=c[complete.cases(c),]
adm.cens=5*365.25
c$fu.time <- pmin(c$t_cvds, adm.cens)
c$status <- ifelse(as.numeric(adm.cens < c$t_cvds), 0, c$cvd)
survcox_ascvd<-coxph(data=c, Surv(fu.time, status)~baseline_age+female+black+tob+
                   sbp+
                   bprx+statin+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_ascvd)
survfit_c=survfit(survcox_ascvd, newdata=c, se.fit=FALSE)
estinc_c=1-survfit_c$surv[dim(survfit_c$surv)[1],]
c$dec=as.numeric(cut2(estinc_c, g=10))
GND.result=GND.calib(pred=estinc_c, tvar=c$fu.time, out=c$status, 
                     cens.t=adm.cens, groups=c$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_c,c$cvd)


##### CVD mortality #####
cvmort = (accord_sets$censor_cm==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_cm*365.25))
t_cvmorts = rowMaxs(cbind(accord_sets$fuyrs_cm*365.25*(1-accord_sets$censor_cm)))
t_cvmorts[t_cvmorts==0] = t_censor[t_cvmorts==0]
t_cvmorts[t_cvmorts==0] = 'NA'
t_cvmorts = as.numeric(t_cvmorts)
cOutcome = Surv(time=t_cvmorts, event = cvmort)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,tob,bmi,
                        baseline_age,female,black,hisp,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,gfr,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
cvmortmodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(cvmortmodel.cv)
coef.cv = coef(cvmortmodel.cv, s = 'lambda.min')
coef.cv
d<-data.frame(cvmort,t_cvmorts,
              baseline_age,female,black,tob,
              sbp,
              bprx,statin,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_cvmorts, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_cvmorts), 0, d$cvmort)
survcox_cvdmort<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+tob+
                   sbp+
                   bprx+statin+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_cvdmort)
survfit_d=survfit(survcox_cvdmort, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$cvmort)


##### MI: fatal or nonfatal #####
fatalmi = (censor_cm==0)&(censor_tst==1)
mi = (accord_sets$censor_nmi==0)|(fatalmi==1)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25,accord_sets$fuyrs_cm*365.25))
t_mis = rowMaxs(cbind(accord_sets$fuyrs_nmi*365.25*(1-accord_sets$censor_nmi),accord_sets$fuyrs_cm*365.25*(1-fatalmi)))
t_mis[t_mis==0] = t_censor[t_mis==0]
t_mis[t_mis==0] = 'NA'
t_mis = as.numeric(t_mis)
cOutcome = Surv(time=t_mis, event = mi)
testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,tob,bmi,
                        baseline_age,female,black,hisp,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,gfr,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
mimodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(mimodel.cv)
coef.cv = coef(mimodel.cv, s = 'lambda.min')
coef.cv
d<-data.frame(mi,t_mis,
              baseline_age,female,black,tob,
              sbp,
              bprx,statin,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_mis, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_mis), 0, d$mi)
survcox_mi<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+tob+
                   sbp+
                   bprx+statin+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_mi)
survfit_d=survfit(survcox_mi, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$mi)

##### stroke: fatal or nonfatal #####
str = (accord_sets$censor_tst==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_tst*365.25))
t_strs = rowMaxs(cbind(accord_sets$fuyrs_tst*365.25*(1-accord_sets$censor_tst)))
t_strs[t_strs==0] = t_censor[t_strs==0]
t_strs[t_strs==0] = 'NA'
t_strs = as.numeric(t_strs)
cOutcome = Surv(time=t_strs, event = str)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,tob,bmi,
                        baseline_age,female,black,hisp,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,gfr,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
strmodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(strmodel.cv)
coef.cv = coef(strmodel.cv, s = 'lambda.min')
coef.cv
d<-data.frame(str,t_strs,
              baseline_age,female,black,tob,
              sbp,
              bprx,statin,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_strs, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_strs), 0, d$str)
survcox_str<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+tob+
                   sbp+
                   bprx+statin+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_str)
survfit_d=survfit(survcox_str, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$str)


##### CHF #####
chf = (accord_sets$censor_chf==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_chf*365.25))
t_chfs = rowMaxs(cbind(accord_sets$fuyrs_chf*365.25*(1-accord_sets$censor_chf)))
t_chfs[t_chfs==0] = t_censor[t_chfs==0]
t_chfs[t_chfs==0] = 'NA'
t_chfs = as.numeric(t_chfs)
cOutcome = Surv(time=t_chfs, event = chf)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,tob,bmi,
                        baseline_age,female,black,hisp,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,gfr,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
chfmodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(chfmodel.cv)
coef.cv = coef(chfmodel.cv, s = 'lambda.min')
coef.cv
d<-data.frame(chf,t_chfs,
              baseline_age,female,black,tob,
              sbp,
              bprx,statin,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_chfs, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_chfs), 0, d$chf)
survcox_chf<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+tob+
                   sbp+
                   bprx+statin+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_chf)
survfit_d=survfit(survcox_chf, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$chf)


##### all cause mort #####
allmort = (accord_sets$censor_tm==0)
t_censor = rowMaxs(cbind(accord_sets$fuyrs_tm*365.25))
t_allmorts = rowMaxs(cbind(accord_sets$fuyrs_tm*365.25*(1-accord_sets$censor_tm)))
t_allmorts[t_allmorts==0] = t_censor[t_allmorts==0]
t_allmorts[t_allmorts==0] = 'NA'
t_allmorts = as.numeric(t_allmorts)
cOutcome = Surv(time=t_allmorts, event = allmort)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,tob,bmi,
                        baseline_age,female,black,hisp,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,gfr,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
allmortmodel.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(allmortmodel.cv)
coef.cv = coef(allmortmodel.cv, s = 'lambda.min')
coef.cv
d<-data.frame(allmort,t_allmorts,
              baseline_age,female,black,tob,
              sbp,
              bprx,statin,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_allmorts, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_allmorts), 0, d$allmort)
survcox_allmort<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+tob+
                   sbp+
                   bprx+statin+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_allmort)
survfit_d=survfit(survcox_allmort, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$allmort)


##### nephropathy 1: SCr doubling or >20 mL/min decr in eGFR #####

neph1 = as.numeric(accord_sets$Neph1==1)
t_censor = rowMaxs(cbind(accord_sets$Neph1Days))
t_neph1s = rowMaxs(cbind(accord_sets$Neph1Days*as.numeric(accord_sets$Neph1)))
t_neph1s[is.na(t_neph1s)]=0
t_neph1s[t_neph1s==0] = t_censor[t_neph1s==0]
t_neph1s[t_neph1s==0] = 'NA'
t_neph1s = as.numeric(t_neph1s)
cOutcome = Surv(time=t_neph1s, event = neph1)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neph1model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neph1model.cv)
coef.cv = coef(neph1model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neph1,t_neph1s,
              baseline_age,female,black,hisp,tob,intensivegly,intensivebp,fibratearm,
              sbp,
              bprx,oraldmrx,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neph1s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neph1s), 0, d$neph1)
survcox_neph1<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+hisp+tob+intensivegly+intensivebp+fibratearm+
                   sbp+
                   bprx+oraldmrx+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neph1)
survfit_d=survfit(survcox_neph1, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neph1)


##### nephropathy 2: development of macro-albuminuria (UAlb>=300) #####

neph2 = as.numeric(accord_sets$Neph2==1)
t_censor = rowMaxs(cbind(accord_sets$Neph2Days))
t_neph2s = rowMaxs(cbind(accord_sets$Neph2Days*as.numeric(accord_sets$Neph2)))
t_neph2s[is.na(t_neph2s)]=0
t_neph2s[t_neph2s==0] = t_censor[t_neph2s==0]
t_neph2s[t_neph2s==0] = 'NA'
t_neph2s = as.numeric(t_neph2s)
cOutcome = Surv(time=t_neph2s, event = neph2)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neph2model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neph2model.cv)
coef.cv = coef(neph2model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neph2,t_neph2s,
              baseline_age,female,black,hisp,tob,intensivegly,intensivebp,fibratearm,
              sbp,
              bprx,oraldmrx,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neph2s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neph2s), 0, d$neph2)
survcox_neph2<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+hisp+tob+intensivegly+intensivebp+fibratearm+
                   sbp+
                   bprx+oraldmrx+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neph2)
survfit_d=survfit(survcox_neph2, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neph2)

##### nephropathy 3: renal failure OR ESRD (dialysis) OR SCr>3.3 #####

neph3 = as.numeric(accord_sets$Neph3==1)
t_censor = rowMaxs(cbind(accord_sets$Neph3Days))
t_neph3s = rowMaxs(cbind(accord_sets$Neph3Days*as.numeric(accord_sets$Neph3)))
t_neph3s[is.na(t_neph3s)]=0
t_neph3s[t_neph3s==0] = t_censor[t_neph3s==0]
t_neph3s[t_neph3s==0] = 'NA'
t_neph3s = as.numeric(t_neph3s)
cOutcome = Surv(time=t_neph3s, event = neph3)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neph3model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neph3model.cv)
coef.cv = coef(neph3model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neph3,t_neph3s,intensivegly,intensivebp,fibratearm,
              baseline_age,female,black,hisp,tob,
              sbp,
              bprx,oraldmrx,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neph3s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neph3s), 0, d$neph3)
survcox_neph3<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+hisp+tob+intensivegly+intensivebp+fibratearm+
                   sbp+
                   bprx+oraldmrx+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neph3)
survfit_d=survfit(survcox_neph3, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neph3)


##### nephropathy 4: any of Nephropathy Outcomes 1-3 #####

neph4 = as.numeric(accord_sets$Neph4==1)
t_censor = rowMaxs(cbind(accord_sets$Neph4Days))
t_neph4s = rowMaxs(cbind(accord_sets$Neph4Days*as.numeric(accord_sets$Neph4)))
t_neph4s[is.na(t_neph4s)]=0
t_neph4s[t_neph4s==0] = t_censor[t_neph4s==0]
t_neph4s[t_neph4s==0] = 'NA'
t_neph4s = as.numeric(t_neph4s)
cOutcome = Surv(time=t_neph4s, event = neph4)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neph4model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neph4model.cv)
coef.cv = coef(neph4model.cv, s = 'lambda.min')
coef.cv
d<-data.frame(neph4,t_neph4s,intensivegly,intensivebp,fibratearm,
              baseline_age,female,black,hisp,tob,
              sbp,
              bprx,oraldmrx,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neph4s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neph4s), 0, d$neph4)
survcox_neph4<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+hisp+tob+intensivegly+intensivebp+fibratearm+
                   sbp+
                   bprx+oraldmrx+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neph4)
survfit_d=survfit(survcox_neph4, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neph4)



##### nephropathy 5: development of micro-albuminuria (UAlb>=30)  #####

neph5 = as.numeric(accord_sets$Neph5==1)
t_censor = rowMaxs(cbind(accord_sets$Neph5Days))
t_neph5s = rowMaxs(cbind(accord_sets$Neph5Days*as.numeric(accord_sets$Neph5)))
t_neph5s[is.na(t_neph5s)]=0
t_neph5s[t_neph5s==0] = t_censor[t_neph5s==0]
t_neph5s[t_neph5s==0] = 'NA'
t_neph5s = as.numeric(t_neph5s)
cOutcome = Surv(time=t_neph5s, event = neph5)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neph5model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neph5model.cv)
coef.cv = coef(neph5model.cv, s = 'lambda.min')
coef.cv
d<-data.frame(neph5,t_neph5s,intensivegly,intensivebp,fibratearm,
              baseline_age,female,black,hisp,tob,
              sbp,
              bprx,oraldmrx,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neph5s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neph5s), 0, d$neph5)
survcox_neph5<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+hisp+tob+intensivegly+intensivebp+fibratearm+
                   sbp+
                   bprx+oraldmrx+anti_coag+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neph5)
survfit_d=survfit(survcox_neph5, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neph5)



##### nephropathy 2/3/5: for DPPOS validation #####

neph235 = (accord_sets$Neph2==1)|(accord_sets$Neph3==1)|(accord_sets$Neph5==1)
t_censor = rowMaxs(cbind(accord_sets$Neph2Days,accord_sets$Neph3Days,accord_sets$Neph5Days))
t_neph235s = rowMaxs(cbind(accord_sets$Neph2Days*(accord_sets$Neph2),accord_sets$Neph3Days*(accord_sets$Neph3),accord_sets$Neph5Days*(accord_sets$Neph5)))
t_neph235s[is.na(t_neph235s)]=0
t_neph235s[t_neph235s==0] = t_censor[t_neph235s==0]
t_neph235s[t_neph235s==0] = 'NA'
t_neph235s = as.numeric(t_neph235s)
cOutcome = Surv(time=t_neph235s, event = neph235)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neph235model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neph235model.cv)
coef.cv = coef(neph235model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neph235,t_neph235s,intensivegly,intensivebp,fibratearm,
              baseline_age,female,black,hisp,tob,
              sbp,
              bprx,oraldmrx,anti_coag,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neph235s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neph235s), 0, d$neph235)
survcox_neph235<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+hisp+tob+intensivegly+intensivebp+fibratearm+
                         sbp+
                         bprx+oraldmrx+anti_coag+
                         cvd_hx_baseline+
                         hba1c+chol+hdl+screat+uacr)
summary(survcox_neph235)
survfit_d=survfit(survcox_neph235, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neph235)



##### Retinopathy Outcome 1 : photocoagulation or vitrectomy  #####

retin1 = as.numeric(accord_sets$Retin1==1)
t_censor = rowMaxs(cbind(accord_sets$Retin1Days))
t_retin1s = rowMaxs(cbind(accord_sets$Retin1Days*as.numeric(accord_sets$Retin1)))
t_retin1s[is.na(t_retin1s)]=0
t_retin1s[t_retin1s==0] = t_censor[t_retin1s==0]
t_retin1s[t_retin1s==0] = 'NA'
t_retin1s = as.numeric(t_retin1s)
cOutcome = Surv(time=t_retin1s, event = retin1)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
retin1model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(retin1model.cv)
coef.cv = coef(retin1model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(retin1,t_retin1s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_retin1s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_retin1s), 0, d$retin1)
survcox_retin1<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                        sbp+
                        bprx+oraldmrx+
                        cvd_hx_baseline+
                        hba1c+chol+hdl+screat+uacr)
summary(survcox_retin1)
survfit_d=survfit(survcox_retin1, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$retin1)


##### Retinopathy Outcome 2: cataract extraction   #####

retin2 = as.numeric(accord_sets$Retin2==1)
t_censor = rowMaxs(cbind(accord_sets$Retin2Days))
t_retin2s = rowMaxs(cbind(accord_sets$Retin2Days*as.numeric(accord_sets$Retin2)))
t_retin2s[is.na(t_retin2s)]=0
t_retin2s[t_retin2s==0] = t_censor[t_retin2s==0]
t_retin2s[t_retin2s==0] = 'NA'
t_retin2s = as.numeric(t_retin2s)
cOutcome = Surv(time=t_retin2s, event = retin2)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
retin2model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(retin2model.cv)
coef.cv = coef(retin2model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(retin2,t_retin2s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_retin2s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_retin2s), 0, d$retin2)
survcox_retin2<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                   sbp+
                   bprx+oraldmrx+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_retin2)
survfit_d=survfit(survcox_retin2, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$retin2)



##### Retinopathy Outcome 3: 3-line change in visual acuity #####

retin3 = as.numeric(accord_sets$Retin3==1)
t_censor = rowMaxs(cbind(accord_sets$Retin3Days))
t_retin3s = rowMaxs(cbind(accord_sets$Retin3Days*as.numeric(accord_sets$Retin3)))
t_retin3s[is.na(t_retin3s)]=0
t_retin3s[t_retin3s==0] = t_censor[t_retin3s==0]
t_retin3s[t_retin3s==0] = 'NA'
t_retin3s = as.numeric(t_retin3s)
cOutcome = Surv(time=t_retin3s, event = retin3)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
retin3model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(retin3model.cv)
coef.cv = coef(retin3model.cv, s = 'lambda.min')
coef.cv
d<-data.frame(retin3,t_retin3s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,insulinrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_retin3s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_retin3s), 0, d$retin3)
survcox_retin3<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                   sbp+
                   bprx+oraldmrx+insulinrx+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_retin3)
survfit_d=survfit(survcox_retin3, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$retin3)




##### Retinopathy Outcome 4: severe vision loss (SF<20/200)  #####

retin4 = as.numeric(accord_sets$Retin4==1)
t_censor = rowMaxs(cbind(accord_sets$Retin4Days))
t_retin4s = rowMaxs(cbind(accord_sets$Retin4Days*as.numeric(accord_sets$Retin4)))
t_retin4s[is.na(t_retin4s)]=0
t_retin4s[t_retin4s==0] = t_censor[t_retin4s==0]
t_retin4s[t_retin4s==0] = 'NA'
t_retin4s = as.numeric(t_retin4s)
cOutcome = Surv(time=t_retin4s, event = retin4)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
retin4model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(retin4model.cv)
coef.cv = coef(retin4model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(retin4,t_retin4s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_retin4s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_retin4s), 0, d$retin4)
survcox_retin4<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                        sbp+
                        bprx+oraldmrx+
                        cvd_hx_baseline+
                        hba1c+chol+hdl+screat+uacr)
summary(survcox_retin4)
survfit_d=survfit(survcox_retin4, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$retin4)


##### Retinopathy Outcome 1/4: for DPPOS validation  #####
retin14 =  (accord_sets$Retin1==1)|(accord_sets$Retin4==1)
t_censor = rowMaxs(cbind(accord_sets$Retin4Days,accord_sets$Retin1Days))
t_retin14s = rowMaxs(cbind(accord_sets$Retin4Days*as.numeric(accord_sets$Retin4),accord_sets$Retin1Days*as.numeric(accord_sets$Retin1)))
t_retin14s[is.na(t_retin14s)]=0
t_retin14s[t_retin14s==0] = t_censor[t_retin14s==0]
t_retin14s[t_retin14s==0] = 'NA'
t_retin14s = as.numeric(t_retin14s)
cOutcome = Surv(time=t_retin14s, event = retin14)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
retin14model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(retin14model.cv)
coef.cv = coef(retin14model.cv, s = 'lambda.min')
coef.cv
d<-data.frame(retin14,t_retin14s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_retin14s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_retin14s), 0, d$retin14)
survcox_retin14<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                         sbp+
                         bprx+oraldmrx+
                         cvd_hx_baseline+
                         hba1c+chol+hdl+screat+uacr)
summary(survcox_retin14)
survfit_d=survfit(survcox_retin14, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$retin14)



##### Neuropathy Outcome 1: MNSI score >2.0 #####

neuro1 = as.numeric(accord_sets$Neuro1==1)
t_censor = rowMaxs(cbind(accord_sets$Neuro1Days))
t_neuro1s = rowMaxs(cbind(accord_sets$Neuro1Days*as.numeric(accord_sets$Neuro1)))
t_neuro1s[is.na(t_neuro1s)]=0
t_neuro1s[t_neuro1s==0] = t_censor[t_neuro1s==0]
t_neuro1s[t_neuro1s==0] = 'NA'
t_neuro1s = as.numeric(t_neuro1s)
cOutcome = Surv(time=t_neuro1s, event = neuro1)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neuro1model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neuro1model.cv)
coef.cv = coef(neuro1model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neuro1,t_neuro1s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neuro1s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neuro1s), 0, d$neuro1)
survcox_neuro1<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                   sbp+
                   bprx+oraldmrx+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neuro1)
survfit_d=survfit(survcox_neuro1, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neuro1)



##### Neuropathy Outcome 2: vibratory sensation loss #####

neuro2 = as.numeric(accord_sets$Neuro2==1)
t_censor = rowMaxs(cbind(accord_sets$Neuro2Days))
t_neuro2s = rowMaxs(cbind(accord_sets$Neuro2Days*as.numeric(accord_sets$Neuro2)))
t_neuro2s[is.na(t_neuro2s)]=0
t_neuro2s[t_neuro2s==0] = t_censor[t_neuro2s==0]
t_neuro2s[t_neuro2s==0] = 'NA'
t_neuro2s = as.numeric(t_neuro2s)
cOutcome = Surv(time=t_neuro2s, event = neuro2)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neuro2model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neuro2model.cv)
coef.cv = coef(neuro2model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neuro2,t_neuro2s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neuro2s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neuro2s), 0, d$neuro2)
survcox_neuro2<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                   sbp+
                   bprx+oraldmrx+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neuro2)
survfit_d=survfit(survcox_neuro2, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neuro2)



##### Neuropathy Outcome 3: ankle jerk loss  #####

neuro3 = as.numeric(accord_sets$Neuro3==1)
t_censor = rowMaxs(cbind(accord_sets$Neuro3Days))
t_neuro3s = rowMaxs(cbind(accord_sets$Neuro3Days*as.numeric(accord_sets$Neuro3)))
t_neuro3s[is.na(t_neuro3s)]=0
t_neuro3s[t_neuro3s==0] = t_censor[t_neuro3s==0]
t_neuro3s[t_neuro3s==0] = 'NA'
t_neuro3s = as.numeric(t_neuro3s)
cOutcome = Surv(time=t_neuro3s, event = neuro3)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neuro3model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neuro3model.cv)
coef.cv = coef(neuro3model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neuro3,t_neuro3s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neuro3s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neuro3s), 0, d$neuro3)
survcox_neuro3<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                   sbp+
                   bprx+oraldmrx+
                   cvd_hx_baseline+
                   hba1c+chol+hdl+screat+uacr)
summary(survcox_neuro3)
survfit_d=survfit(survcox_neuro3, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neuro3)


##### Neuropathy Outcome 4: pressure sensation loss #####


neuro4 = as.numeric(accord_sets$Neuro4==1)
t_censor = rowMaxs(cbind(accord_sets$Neuro4Days))
t_neuro4s = rowMaxs(cbind(accord_sets$Neuro4Days*as.numeric(accord_sets$Neuro4)))
t_neuro4s[is.na(t_neuro4s)]=0
t_neuro4s[t_neuro4s==0] = t_censor[t_neuro4s==0]
t_neuro4s[t_neuro4s==0] = 'NA'
t_neuro4s = as.numeric(t_neuro4s)
cOutcome = Surv(time=t_neuro4s, event = neuro4)

testsubset = data.frame(cOutcome,intensivegly,intensivebp,fibratearm,
                        baseline_age,female,black,hisp,tob,bmi,
                        sbp,dbp,hr,
                        bprx,oraldmrx,insulinrx,statin,fibrate,anti_coag,anti_inflam,platelet_agi,aspirin,
                        cvd_hx_baseline,
                        hba1c,chol,trig,ldl,hdl,fpg,alt,cpk,potassium.y,screat,ualb,ucreat,uacr)
testsubset=testsubset[complete.cases(testsubset),]
neuro4model.cv =  cv.glmnet(as.matrix(testsubset[,-c(1)]),as.matrix(testsubset[,1]),family="cox",parallel=TRUE)
plot(neuro4model.cv)
coef.cv = coef(neuro4model.cv, s = 'lambda.1se')
coef.cv
d<-data.frame(neuro4,t_neuro4s,
              baseline_age,female,black,
              sbp,
              bprx,oraldmrx,
              cvd_hx_baseline,
              hba1c,chol,hdl,screat,uacr)
d=d[complete.cases(d),]
adm.cens=5*365.25
d$fu.time <- pmin(d$t_neuro4s, adm.cens)
d$status <- ifelse(as.numeric(adm.cens < d$t_neuro4s), 0, d$neuro4)
survcox_neuro4<-coxph(data=d, Surv(fu.time, status)~baseline_age+female+black+
                        sbp+
                        bprx+oraldmrx+
                        cvd_hx_baseline+
                        hba1c+chol+hdl+screat+uacr)
summary(survcox_neuro4)
survfit_d=survfit(survcox_neuro4, newdata=d, se.fit=FALSE)
estinc_d=1-survfit_d$surv[dim(survfit_d$surv)[1],]
d$dec=as.numeric(cut2(estinc_d, g=10))
GND.result=GND.calib(pred=estinc_d, tvar=d$fu.time, out=d$status, 
                     cens.t=adm.cens, groups=d$dec, adm.cens=adm.cens)
GND.result
ci.cvAUC(estinc_d,d$neuro4)

save.image("~/Data/accord/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/accord_dm_models.RData")


