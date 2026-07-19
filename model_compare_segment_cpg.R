#try it for T1, with one file
path_try="/scratch/project_2010912/ezel/mixed_models/meth_MixMod_T1/meth_MixMod_T1_160001to240000"

e <- new.env()
load(path_try, envir = e)
#1st, 4th, 7th etc have the chr and nucleotide num info


model_results <- e$meth_MixMod_T1_160001to240000

#mod1.remlT.list = if  ( length(mod1.remlT)==1  ) NA else 
#list ( ranef(mod1.remlT) , summary(mod1.remlT)$coef$cond , mod1.remlT$fit$message , 
#mod1.remlT$fit$convergence , AIC(mod1.remlT) , residuals(mod1.remlT) )

#has estimated, stds, p values etc
coefs_weight_interc=(model_results[[2]])[[2]]
stat_code=(model_results[[2]])[[4]]
info_crit=(model_results[[2]])[[5]]