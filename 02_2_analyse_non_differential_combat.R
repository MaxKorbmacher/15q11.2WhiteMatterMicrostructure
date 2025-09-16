# Analysis of white matter microstructure in CNVs [15q11.2 deletion and duplication carriers]
#
# Max Korbmacher, July 2025
# max.korbmacher@gmail.com
#
#
# ----------------------------------- #
# --------------Structure------------ #
# ----------------------------------- #
# 0. Data wrangling------------------ #
# 1. Demographics ------------------- #
## 1.1 Duplication ------------------ #
## 1.2 Deletion --------------------- #
# 2. PCA training ------------------- #
## 2.1 Tracts------------------------ #
## 2.2 Regions----------------------- #
# 3. PCA predictions ---------------- #
## 3.1 Tracts------------------------ #
## 3.2 Regions----------------------- #
# 4. Whole brain assessment---------- #
## 4.1 Duplication ------------------ #
## 4.2 Deletion --------------------- #
# 5. Tract-level assessment---------- #
## 5.1 Duplication ------------------ #
## 5.2 Deletion --------------------- #
# 6. Region-level assessment--------- #
## 6.1 Duplication ------------------ #
## 6.2 Deletion --------------------- #
# ----------------------------------- #
# ----------------------------------- #
#
# 0. Data wrangling------------------ 
# wash your hands before eating
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# define the savepath
savepath = "/cluster/projects/p33/users/maxk/UKB/CNV/results/all_ComBat_"
#
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,lme4,lmerTest,haven,reshape2,
               ggseg3d, MatchIt, pwr,sjPlot,marginaleffects,
               vtable, effectsize, fastICA,factoextra,ggpubr,
               ggridges,viridis
               )
# read data
cross = read.csv("/cluster/projects/p33/users/maxk/UKB/CNV/data/cross_Combat_together.csv")
#
# match data
#cross$binarized_groups = ifelse(cross$groups == "HC","Control","Case")
# duplication = cross %>% filter(!groups == "Deletion")
# duplication$groups = factor(duplication$groups , ordered = FALSE )
# duplication = matchit(groups ~ age + sex + site + EstimatedTotalIntraCranialVol,
#                    data = duplication,
#                    method = "optimal",
#                    distance = "glm",
#                    estimand = "ATC",
#                    #exact = c("site","sex"),
#                    ratio = 10)
# duplication = match.data(duplication)
# table(duplication$groups)
#
#
deletion = cross %>% filter(!groups == "Duplication")
deletion$groups = factor(deletion$groups , ordered = FALSE )
deletion = matchit(groups ~ age + sex + site + EstimatedTotalIntraCranialVol,
                   data = deletion,
                   method = "optimal",
                   distance = "glm",
                   estimand = "ATC",
                   #exact = c("site","sex"),
                   ratio = 10)
deletion = match.data(deletion)

# 1. Demographics ------------------- 
## 1.1 Duplication ------------------
# continuous
# st(duplication %>% filter(groups == "Duplication")%>% select(age,sex,site),out = 'csv', 
#    file = paste(savepath,"Dup_demo",sep=""))
# st(duplication %>% filter(groups == "HC")%>% select(age,sex,site),out = 'csv', 
#    file = paste(savepath,"Dup_demo_HC",sep=""))
#
## 1.2 Deletion ---------------------
st(deletion %>% filter(groups == "Deletion")%>% select(age,sex,site),out = 'csv', 
   file = paste(savepath,"Del_demo",sep=""))
st(deletion %>% filter(groups == "HC")%>% select(age,sex,site),out = 'csv', 
   file = paste(savepath,"Del_demo_HC",sep=""))


# 2. PCA training ------------------- 
# pca_dup = dMRI_cross[!dMRI_cross$eid %in% duplication$eid,]
pca_del = cross[!cross$eid %in% deletion$eid,]

## 2.1 Tracts------------------------- 
tracts = deletion %>% select(ends_with("ATRL"),ends_with("ATRR"),
                                ends_with("CSTL"),ends_with("CSTR"),
                                ends_with("CINGL"),ends_with("CINGR"),
                                ends_with("CGL"),ends_with("CGR"),
                                ends_with("IFOFL"),ends_with("IFOFR"),
                                ends_with("ILFL"),ends_with("ILFR"),
                                ends_with("SLFL"),ends_with("SLFR"),
                                ends_with("UFL"),ends_with("UFR"),
                                ends_with("SLFTL"),ends_with("SLTFR"),
                                ends_with("FMIN"),ends_with("FMAJ")
) %>% names
#
# Visualisation
# estimate principal components and visualise the first 10 components
pcplots = function(data){
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
  a = fviz_eig(dkipc, main = "DKI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
  b = fviz_eig(dtipc, main = "DTI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
  c = fviz_eig(briapc, main = "BRIA", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
  d = fviz_eig(smtpc, main = "SMT", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
  smtmcpc =  data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T)
  e = fviz_eig(smtmcpc, main = "SMTmc", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T)
  f = fviz_eig(wmtipc, main = "WMTI", addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
  plot4 = ggarrange(a,b,c,d,e,f)
  return(plot4)
}
plot4.1 = pcplots(pca_del[tracts])
plot4.2 = pcplots(pca_dup[tracts])
ggsave(paste(savepath,"scree_tracts_deletion.pdf",sep=""),plot4.1,width=14,height=10)
# ggsave(paste(savepath,"scree_tracts_duplication.pdf",sep=""),plot4.2,width=14,height=10)
# rm(plot4.1,plot4.2)
#
#
## 2.2 Regions------------------------ 
regions = names(deletion)[!names(deletion) %in% tracts]
regions = regions[!regions %in% (deletion %>% select(contains("Mean")) %>% names)]
regions = regions[!regions %in% c("eid","X","X.1","age","sex","site","groups", "distance",
                                  "weights","subclass","EstimatedTotalIntraCranialVol")]
plot4.1 = pcplots(pca_del[regions])
# plot4.2 = pcplots(pca_dup[regions])
ggsave(paste(savepath,"scree_regions_deletion.pdf",sep=""),plot4.1,width=14,height=10)
# ggsave(paste(savepath,"scree_regions_duplication.pdf",sep=""),plot4.2,width=14,height=10)
# rm(plot4.1,plot4.2)
# 3. PCA predictions ----------------
create_pc = function(data){
  dkipc = data %>% select(starts_with("rk"), starts_with("ak"), starts_with("mk")) %>% prcomp(scale. = T, center = T)
  dtipc =  data %>% select(starts_with("rd"), starts_with("ad"), starts_with("md"), starts_with("fa")) %>% prcomp(scale. = T, center = T)
  briapc =  data %>% select(starts_with("v_"), starts_with("micro"), starts_with("Drad"), starts_with("Dax")) %>% prcomp(scale. = T, center = T)
  smtpc =  data %>% select(starts_with("smt_md"), starts_with("smt_long")) %>% prcomp(scale. = T, center = T)
  smtmcpc =  data %>% select(starts_with("smt_mc")) %>% prcomp(scale. = T, center = T)
  wmtipc =  data %>% select(starts_with("axEAD"), starts_with("radEAD")) %>% prcomp(scale. = T, center = T)
  pc_list = list(dkipc,dtipc,briapc,smtpc,smtmcpc,wmtipc)
  return(pc_list)
}
## 3.1 Tracts------------------------- 
# train
# duplication_tract_pcs = create_pc(pca_dup[tracts])
deletion_tract_pcs = create_pc(pca_del[tracts])
# predict
# duplication$DKI_PC = data.frame(predict(duplication_tract_pcs[[1]],newdata = duplication))$PC1
# duplication$DTI_PC = data.frame(predict(duplication_tract_pcs[[2]],newdata = duplication))$PC1
# duplication$BRIA_PC = data.frame(predict(duplication_tract_pcs[[3]],newdata = duplication))$PC1
# duplication$SMT_PC = data.frame(predict(duplication_tract_pcs[[4]],newdata = duplication))$PC1
# duplication$SMTmc_PC = data.frame(predict(duplication_tract_pcs[[5]],newdata = duplication))$PC1
# duplication$WMTI_PC = data.frame(predict(duplication_tract_pcs[[6]],newdata = duplication))$PC1
deletion$DKI_PC = data.frame(predict(deletion_tract_pcs[[1]],newdata = deletion))$PC1
deletion$DTI_PC = data.frame(predict(deletion_tract_pcs[[2]],newdata = deletion))$PC1
deletion$BRIA_PC = data.frame(predict(deletion_tract_pcs[[3]],newdata = deletion))$PC1
deletion$SMT_PC = data.frame(predict(deletion_tract_pcs[[4]],newdata = deletion))$PC1
deletion$SMTmc_PC = data.frame(predict(deletion_tract_pcs[[5]],newdata = deletion))$PC1
deletion$WMTI_PC = data.frame(predict(deletion_tract_pcs[[6]],newdata = deletion))$PC1
PCs = deletion %>% select(ends_with("PC")) %>% names
# associate
est = ci_l = ci_h = p = c()
# for (i in 1:length(PCs)){
#   fit <- lm(paste(PCs[i], "~ age+age^2+sex+groups"), data = duplication,weights = weights)
#   print(paste("Model",i,"done"))
#   est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
#   ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
#   ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
#   p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
# }
# dupl_tract_res = data.frame(PCs,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
# dupl_tract_res %>% filter(pFDR < 0.05)
# write.csv(x = dupl_tract_res, file = paste(savepath,"dupl_PC_tract_res.csv",sep=""))
for (i in 1:length(PCs)){
  fit <- lm(paste(PCs[i], "~ age+age^2+sex+groups"), data = deletion,weights = weights)
  print(paste("Model",i,"done"))
  est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
  ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
  ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
  p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
}
del_tract_res = data.frame(PCs,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
del_tract_res %>% filter(pFDR < 0.05)
write.csv(x = del_tract_res, file = paste(savepath,"del_PC_tract_res.csv",sep=""))
#
#
## 3.2 Regions------------------------
# duplication_region_pcs = create_pc(duplication[regions])
deletion_region_pcs = create_pc(deletion[regions])
# predict
# duplication$DKI_PC = data.frame(predict(duplication_region_pcs[[1]],newdata = duplication))$PC1
# duplication$DTI_PC = data.frame(predict(duplication_region_pcs[[2]],newdata = duplication))$PC1
# duplication$BRIA_PC = data.frame(predict(duplication_region_pcs[[3]],newdata = duplication))$PC1
# duplication$SMT_PC = data.frame(predict(duplication_region_pcs[[4]],newdata = duplication))$PC1
# duplication$SMTmc_PC = data.frame(predict(duplication_region_pcs[[5]],newdata = duplication))$PC1
# duplication$WMTI_PC = data.frame(predict(duplication_region_pcs[[6]],newdata = duplication))$PC1
deletion$DKI_PC = data.frame(predict(deletion_region_pcs[[1]],newdata = deletion))$PC1
deletion$DTI_PC = data.frame(predict(deletion_region_pcs[[2]],newdata = deletion))$PC1
deletion$BRIA_PC = data.frame(predict(deletion_region_pcs[[3]],newdata = deletion))$PC1
deletion$SMT_PC = data.frame(predict(deletion_region_pcs[[4]],newdata = deletion))$PC1
deletion$SMTmc_PC = data.frame(predict(deletion_region_pcs[[5]],newdata = deletion))$PC1
deletion$WMTI_PC = data.frame(predict(deletion_region_pcs[[6]],newdata = deletion))$PC1
PCs = deletion %>% select(ends_with("PC")) %>% names
# associate
est = ci_l = ci_h = p = c()
# for (i in 1:length(PCs)){
#   fit <- lm(paste(PCs[i], "~ age+age^2+sex+groups"), data = duplication,weights = weights)
#   print(paste("Model",i,"done"))
#   est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
#   ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
#   ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
#   p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
# }
# dupl_tract_res = data.frame(PCs,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
# dupl_tract_res %>% filter(pFDR < 0.05)
# write.csv(x = dupl_tract_res, file = paste(savepath,"dupl_PC_region_res.csv",sep=""))
for (i in 1:length(PCs)){
  fit <- lm(paste(PCs[i], "~ age+age^2+sex+groups"), data = deletion,weights = weights)
  print(paste("Model",i,"done"))
  est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
  ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
  ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
  p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
}
del_tract_res = data.frame(PCs,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
del_tract_res %>% filter(pFDR < 0.05)
write.csv(x = del_tract_res, file = paste(savepath,"del_PC_region_res.csv",sep=""))
#
# 4. Whole brain assessment----------
skeletons = deletion %>% select(contains("Mean")) %>% names
## 4.1 Duplication ------------------
# duplication = duplication %>% mutate(groups = relevel(groups, "HC"))
# est = ci_l = ci_h = p = c()
# for (i in 1:length(skeletons)){
#   fit <- lm(paste(skeletons[i], "~ age+age^2+sex+groups"), data = duplication,weights = weights)
#   print(paste("Model",i,"done"))
#   est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
#   ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
#   ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
#   p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
# }
# dupl_skeleton_res = data.frame(skeletons,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
# dupl_skeleton_res %>% filter(pFDR < 0.05)
# write.csv(x = dupl_skeleton_res, file = paste(savepath,"dupl_skeleton_res.csv",sep=""))
#
#
## 4.2 Deletion ---------------------
deletion = deletion %>% mutate(groups = relevel(groups, "HC"))
est = ci_l = ci_h = p = c()
for (i in 1:length(skeletons)){
  fit <- lm(paste(skeletons[i], "~ age+age^2+sex+groups"), data = deletion,weights = weights)
  print(paste("Model",i,"done"))
  est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
  ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
  ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
  p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
}
del_skeleton_res = data.frame(skeletons,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
del_skeleton_res %>% filter(pFDR < 0.05)
write.csv(x = del_skeleton_res, file = paste(savepath,"del_skeleton_res.csv",sep=""))

# 5. Tract-level assessment---------- 
## 5.1 Duplication ------------------
tracts = deletion %>% select(ends_with("ATRL"),ends_with("ATRR"),
                                ends_with("CSTL"),ends_with("CSTR"),
                                ends_with("CINGL"),ends_with("CINGR"),
                                ends_with("CGL"),ends_with("CGR"),
                                ends_with("IFOFL"),ends_with("IFOFR"),
                                ends_with("ILFL"),ends_with("ILFR"),
                                ends_with("SLFL"),ends_with("SLFR"),
                                ends_with("UFL"),ends_with("UFR"),
                                ends_with("SLFTL"),ends_with("SLTFR"),
                                ends_with("FMIN"),ends_with("FMAJ")
                                ) %>% names
est = ci_l = ci_h = p = c()
# for (i in 1:length(tracts)){
#   fit <- lm(paste(tracts[i], "~ age+age^2+sex+groups"), data = duplication,weights = weights)
#   print(paste("Model",i,"done"))
#   est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
#   ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
#   ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
#   p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
# }
# dupl_tract_res = data.frame(tracts,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
# dupl_tract_res %>% filter(pFDR < 0.05)
# write.csv(x = dupl_tract_res, file = paste(savepath,"dupl_tract_res.csv",sep=""))

## 5.2 Deletion ---------------------
est = ci_l = ci_h = p = c()
for (i in 1:length(tracts)){
  fit <- lm(paste(tracts[i], "~ age+age^2+sex+groups"), data = deletion,weights = weights)
  print(paste("Model",i,"done"))
  est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
  ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
  ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
  p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
}
del_tract_res = data.frame(tracts,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
del_tract_res %>% filter(pFDR < 0.05)
write.csv(x = del_tract_res, file = paste(savepath,"del_tract_res.csv",sep=""))
#
# 6. Region-wise  ---------------------
## 6.1 Duplication ------------------ 
regions = names(deletion%>%select(!contains("_PC")))[!names(deletion) %in% tracts]
regions = regions[!regions %in% skeletons]
regions = regions[!regions %in% c("eid","X","X.1","age","sex","site","groups", "distance",
                                  "weights","subclass","EstimatedTotalIntraCranialVol")]
regions = regions[1:1248]
# est = ci_l = ci_h = p = c()
# for (i in 1:length(regions)){
#   fit <- lm(paste(regions[i], "~ age+age^2+sex+groups"), data = duplication,weights = weights)
#   print(paste("Model",i,"done"))
#   est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
#   ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
#   ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
#   p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
# }
# dupl_regions_res = data.frame(regions,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
# dupl_regions_res %>% filter(pFDR < 0.05)
# write.csv(x = dupl_regions_res, file = paste(savepath,"dupl_regions_res.csv",sep=""))

## 6.2 Deletion --------------------- 
est = ci_l = ci_h = p = c()
for (i in 1:length(regions)){
  fit <- lm(paste(regions[i], "~ age+age^2+sex+groups"), data = deletion,weights = weights)
  print(paste("Model",i,"done"))
  est[i] = effectsize::standardize_parameters(fit)$Std_Coefficient[4]
  ci_l[i] = effectsize::standardize_parameters(fit)$CI_low[4]
  ci_h[i] = effectsize::standardize_parameters(fit)$CI_high[4]
  p[i] = as.numeric(avg_comparisons(fit,variables = "groups",vcov = ~subclass,wts = "weights")[6])
}
del_regions_res = data.frame(regions,est, ci_l, ci_h, p, pFDR = p.adjust(p,method = "fdr"))
del_regions_res %>% filter(pFDR < 0.05)
write.csv(x = del_regions_res, file = paste(savepath,"del_regions_res.csv",sep=""))


p_frame = del_regions_res %>% filter(pFDR < 0.05)

p_frame$regions[p_frame$regions =="AD_BodyCC"] = "DTI-AD BodyCC"
p_frame$regions[p_frame$regions =="MD_BodyCC"] = "DTI-MD BodyCC"
p_frame$regions[p_frame$regions =="AD_GenuCC"] = "DTI-AD GenyCC"
p_frame$regions[p_frame$regions =="RD_BodyCC"] = "DTI-RD BodyCC"
p_frame$regions[p_frame$regions =="RD_UncinatefasciculusR"] = "DTI-RD rUF"
p_frame$regions[p_frame$regions =="FA_UncinatefasciculusR"] = "DTI-FA rUF"
p_frame$regions[p_frame$regions =="AD_BodyCC"] = "DKI-AK SpleniumCC"
p_frame$regions[p_frame$regions =="awf_BodyCC"] = "WMTI-AWF BodyCC"
p_frame$regions[p_frame$regions =="micro_ADC_BodyCC"] = "BRIA-microADC BodyCC"
p_frame$regions[p_frame$regions =="micro_Rd_BodyCC"] = "BRIA-microRD BodyCC"
p_frame$regions[p_frame$regions =="radEAD_BodyCC"] = "WMTI-radEAD BodyCC"

p_frame$regions[p_frame$regions =="smt_fa_BodyCC"] = "SMT-FA BodyCC"
p_frame$regions[p_frame$regions =="smt_mc_extramd_BodyCC"] = "SMTmc-extraMD BodyCC"
p_frame$regions[p_frame$regions =="smt_mc_extratrans_BodyCC"] = "SMTmc-extraTrans BodyCC"
p_frame$regions[p_frame$regions =="smt_mc_extratrans_CingulumhippocampusR"] = "SMTmc-extraTrans rCH"
p_frame$regions[p_frame$regions =="smt_mc_extratrans_UncinatefasciculusR"] = "SMTmc-extraTrans rUF"
p_frame$regions[p_frame$regions =="smt_mc_intra_BodyCC"] = "SMTmc-intra BodyCC"
p_frame$regions[p_frame$regions =="smt_md_BodyCC"] = "SMT-MD BodyCC"
p_frame$regions[p_frame$regions =="ak_SpleniumCC"] = "DKI-AK SpleniumCC"
p_frame$regions[p_frame$regions =="smt_trans_BodyCC"] = "SMT-Trans BodyCC"

plot = ggplot(p_frame, aes(y=est, x=regions)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=ci_l, ymax=ci_h), width=.4,position=position_dodge(.9)) +
  ylab("Standardized Estimate with 95% Confidence Interval")+xlab("Metric in Region") + 
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw()
ggsave(filename = paste(savepath,"plot.pdf",sep=""),plot = plot, width = 10,height = 6)

# Now, plot raw values of the observed significant group differences
regs = del_regions_res %>% filter(pFDR < 0.05) %>% pull(regions)

#regs = c("radEAD_BodyCC","smt_mc_extratrans_BodyCC") # for test purposes
selected_data = deletion %>% select(all_of(regs), groups)

# We need to reshape the data from 'wide' to 'long' format to use ggplot2 effectively.
# This makes it easy to plot both metrics on the same y-axis (facet-wrapped).
data_long <- tidyr::pivot_longer(selected_data, 
                                 cols = all_of(regs), 
                                 names_to = "metric", 
                                 values_to = "value")
# Reorder the groups to ensure a logical plotting order
#data_long$groups <- factor(data_long$groups, levels = rev(groups))

# --- 3. Create the Ridgeline Plot ---
ridge = ggplot(data_long, aes(x = value, y = groups, fill = metric)) +
  geom_density_ridges(alpha = 0.8, rel_min_height = 0.01) +
  scale_fill_viridis_d(option = "plasma") + # Use viridis discrete scale
  #scale_fill_manual(values = c("#FF6B6B", "#4ECDC4")) + # Beautiful custom colors
  facet_wrap(~metric, scales = "free_x", ncol = 2) + # Separate plots for each metric
  labs(
    #title = "Distribution of Metrics Across 20 Groups",
    #subtitle = "Ridgeline plot for visualizing group differences",
    x = "Metric Value",
    y = "Group"
  ) +
  theme_minimal() + # A clean, minimalist theme
  #scale_color_viridis(discrete = TRUE, option = "D")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(t = 10, b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 20)),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(size = 8),
    legend.position = "none", # Hide legend since facets are labeled
    panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(2, "lines") # Add extra space between facets
  )

ggsave("/tsd/p33/home/p33-maxk/ridge.pdf",ridge,width=10,height=15)
ggsave(paste(savepath,"ridge.pdf",sep=""),ridge,width=10,height=15)

