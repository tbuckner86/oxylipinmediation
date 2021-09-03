library(regmedint)
library(boot)

set.seed(1111)

# Load data
#laptop
#setwd("C:/Users/tbuck/OneDrive - The University of Colorado Denver/daisy/OXYLIPINS/METHYLATION MEDIATION/dataset")
#desktop
setwd("/home/bucknete")

#laptop
#load("C:/Users/tbuck/OneDrive - The University of Colorado Denver/daisy/OXYLIPINS/METHYLATION MEDIATION/dataset/sv_pret1d_dataset.Rdata")
#desktop
load("sv_pret1d_dataset.Rdata")


# Outcome and adjustment variables
ia = as.numeric(factor(sv$IAgroup2)) - 1
t1d = as.numeric(factor(sv$T1Dgroup)) - 1
SEX = as.numeric(factor(sv$SEX)) - 1
dr34 = sv$dr34
age_delta = as.numeric(pret1d$clinage) - as.numeric(sv$clinage)
age = as.numeric(sv$clinage)
covariates = as.data.frame(cbind(t1d,SEX,dr34,age,age_delta))


# Pairs
#oxylipins 17, 23, 37, 40, 41, and 42 were not on PCs
load("probesFromPipeline.Rdata")
oxylipins = c(paste0("bc_oxylipin",c(1:16,18:22,24:36,38,39,43:48)),"Prin1","Prin2")

pairs = expand.grid(probesFromPipeline,oxylipins,stringsAsFactors = F)


# Bootstrap function A
regmed_boot = function(d,i){
  regmedint_obj = 
    regmedint(data = d[i,],
              ## Variables
              yvar = "t1d",
              avar = "metab",
              mvar = "methyl",
              cvar = c("SEX","dr34","age","age_delta"),
              ## Values at which effects are evaluated
              a0 = 0,
              a1 = 1,
              m_cde = 1,
              c_cond = c(1,1,1,1),
              ## Model types
              mreg = "linear",
              yreg = "logistic",
              ## Additional specification
              interaction = T,
              casecontrol = T)
  summary(regmedint_obj)$summary_myreg["tnie",1]
}


# Bootstrap options
#error here: 
# Error in boot(data = df, statistic = regmed_boot, R = boots, parallel = "multicore",  : 
#no data in call to 'boot' 
boot_cores = 16
boots = 1000
# Iterate through all
methyl_pt1d_results = apply(pairs,1,function(r){
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(pret1d[,methyl]))
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(sv[,metab]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,covariates))
  df = df[complete.cases(df),]
  # Bootstrap
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  ci = boot.ci(b,type = "perc")
  ci_table = c(ci$t0,ci$percent[,4:5])
  return(ci_table)
})
methyl_pt1d_results = t(methyl_pt1d_results)
rownames(methyl_pt1d_results) = apply(pairs,1,paste,collapse = " & ")
colnames(methyl_pt1d_results) = c("Est.","Lower","Upper")
# Save
save(methyl_pt1d_results,file = "methyl_pt1d_results.Rdata")


################################################################

################################################################


# Bootstrap function B
regmed_boot = function(d,i){
  regmedint_obj = 
    regmedint(data = d[i,],
              ## Variables
              yvar = "t1d",
              avar = "methyl",
              mvar = "metab",
              cvar = c("SEX","dr34","age","age_delta"),
              ## Values at which effects are evaluated
              a0 = 0,
              a1 = 1,
              m_cde = 1,
              c_cond = c(1,1,1,1),
              ## Model types
              mreg = "linear",
              yreg = "logistic",
              ## Additional specification
              interaction = T,
              casecontrol = T)
  summary(regmedint_obj)$summary_myreg["tnie",1]
}

# Iterate through all
metab_pt1d_results = apply(pairs,1,function(r){
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(sv[,methyl]))
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(pret1d[,metab]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,covariates))
  df = df[complete.cases(df),]
  # Bootstrap
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  ci = boot.ci(b,type = "perc")
  ci_table = c(ci$t0,ci$percent[,4:5])
  return(ci_table)
})
metab_pt1d_results = t(metab_pt1d_results)
rownames(metab_pt1d_results) = apply(pairs,1,paste,collapse = " & ")
colnames(metab_pt1d_results) = c("Est.","Lower","Upper")
# Save
save(metab_pt1d_results,file = "metab_pt1d_results.Rdata")
