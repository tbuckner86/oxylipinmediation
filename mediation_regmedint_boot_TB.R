library(regmedint)
library(boot)
set.seed(111)
# Load data
setwd("C:/Users/bucknete/OneDrive - The University of Colorado Denver/daisy/OXYLIPINS/METHYLATION MEDIATION")
#setwd("~/Dropbox/School/MS Thesis")


psv <- read.csv("C:/Users/bucknete/OneDrive - The University of Colorado Denver/daisy/OXYLIPINS/METHYLATION MEDIATION/scoredpsvage.csv", na.strings="")
sv <- read.csv ("C:/Users/bucknete/OneDrive - The University of Colorado Denver/daisy/OXYLIPINS/METHYLATION MEDIATION/scoredsvage.csv", na.strings="")

load("C:/Users/bucknete/OneDrive - The University of Colorado Denver/daisy/OXYLIPINS/METHYLATION MEDIATION/psv_sv_dataset.RData")

#load("./data/raw_data/probesFromPipeline.Rdata")
#load("./data/raw_data/psv_sv_dataset.Rdata")
#load("./data/mediation/methyl_psv_candidates_p_01.Rdata")
#colnames(methyl_psv_candidates) = c("Var1","Var2")

#load("./data/mediation/metab_psv_candidates_p_01.Rdata")
#colnames(metab_psv_candidates) = c("Var1","Var2")


# Outcome and adjustment variables
ia = as.numeric(factor(psv$IAgroup2)) - 1
sex = as.numeric(factor(psv$SEX)) - 1
dr34 = psv$dr34
age_delta = as.numeric(sv$clinage) - as.numeric(psv$clinage)
age = as.numeric(psv$clinage)
covariates = as.data.frame(cbind(ia,sex,dr34,age,age_delta))

# Pairs
load("C:/Users/bucknete/OneDrive - The University of Colorado Denver/daisy/OXYLIPINS/METHYLATION MEDIATION/probesFromPipeline.RData")
oxylipins = paste0("bc_oxylipin",1:48)

pairs = expand.grid(probesFromPipeline,oxylipins,stringsAsFactors = F)

# Bootstrap function
regmed_boot = function(d,i){
  regmedint_obj = 
    regmedint(data = d[i,],
              ## Variables
              yvar = "ia",
              avar = "methyl",
              mvar = "metab",
              cvar = c("sex","dr34","age","age_delta"),
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
  summary(regmedint_obj)$summary_myreg[,1]
}
# Bootstrap options
boot_cores = 16
boots = 1000
# Testing
pairs = pairs[1:5,]
# Iterate through all
methyl_psv_results = apply(pairs,1,function(r){
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(psv[,methyl]))
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(sv[,metab]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,covariates))
  df = df[complete.cases(df),]
  # Bootstrap
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  return(b)
})
names(methyl_psv_results) = apply(pairs,1,paste,collapse = " & ")
# Save
save(methyl_psv_results,file = "methyl_psv_results.Rdata")


















# Bootstrap function
regmed_boot = function(d,i){
  regmedint_obj = 
    regmedint(data = d[i,],
              ## Variables
              yvar = "ia",
              avar = "metab",
              mvar = "methyl",
              cvar = c("sex","dr34","age","age_delta"),
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
  summary(regmedint_obj)$summary_myreg[,1]
}

# Iterate through all
metab_psv_results = apply(pairs,1,function(r){
  methyl = as.character(r["Var1"])
  methyl = as.numeric(scale(sv[,methyl]))
  metab = as.character(r["Var2"])
  metab = as.numeric(scale(psv[,metab]))
  # Dataframe 
  df = as.data.frame(cbind(methyl,metab,covariates))
  df = df[complete.cases(df),]
  # Bootstrap
  b = boot(data = df, statistic = regmed_boot, R = boots,parallel = "multicore",
           ncpus = boot_cores)
  return(b)
})
names(metab_psv_results) = apply(pairs,1,paste,collapse = " & ")
# Save
save(metab_psv_results,file = "metab_psv_results.Rdata")
