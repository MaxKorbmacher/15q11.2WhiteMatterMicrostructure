# UKB prep for white matter microstructure differences in rare copy number variants 
# July 2025
# Code by Max Korbmacher (max.korbmacher@gmail.com)
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# CNV
cnv = read.delim("/cluster/projects/p33/groups/imaging/ukbio/genetics/CNVs/returndata/ukbreturn1701/All_CNVs_for_UKBB.dat")
key = read.delim(sep = " ",header = F,file="/cluster/projects/p33/groups/imaging/ukbio/genetics/CNVs/returndata/ukb27412bridge14421.txt")
#pheno = read.csv("/cluster/projects/p33/users/maxk/UKB/environment/data/demo/environment.txt.csv")
dMRI_cross = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/cross_sectional.csv")
#dMRI_T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/T1.csv")
#dMRI_T2 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/brainage/T2.csv")
T1w = read.csv("/cluster/projects/p33/users/maxk/UKB/data/T1w_50k/merged/T1_data.csv")
# lib
library(dplyr)
#install.packages("sva") # in case of neuroCombat not working
library(sva) # for ComBat()

# fix eid labels
names(key) = c("eid","f.eid")
cnv0 = merge(key,cnv,by="f.eid")
cnv = cnv0 %>% filter(filter == "Selected")

# make grouping for case-control
table(cnv$Pathogenic_CNVs)
cnv$groups = ifelse(cnv$Pathogenic_CNVs == "15q11.2del", "Deletion","discard")
cnv$groups = ifelse(cnv$Pathogenic_CNVs == "15q11.2dup", "Duplication",cnv$groups)
patho_CNV = unique((cnv %>% filter(NumCNV>0))$eid)
cnv = cnv %>% filter(!groups == "discard")
# filter only CNVs
cross = merge(cnv %>% select(eid,groups),dMRI_cross, by = "eid")

# add intra-cranial volume (will be used as a covariate)
cross = merge(T1w%>%select(eid,EstimatedTotalIntraCranialVol),cross, by="eid")
dMRI_cross = merge(T1w%>%select(eid,EstimatedTotalIntraCranialVol),dMRI_cross, by="eid")

#
# harmonise each dataset
## Deletion
del = cross %>% filter(groups == "Deletion")
covars = del %>% dplyr::select(eid,sex,site,age,groups)
del = del %>% dplyr::select(-all_of(names(covars))) %>% dplyr::select(-X)
del = t(as.matrix(del))
del = t(ComBat(del,batch=as.numeric(factor(covars$site)),mod=model.matrix(~age+as.numeric(sex), data = covars)))
del = cbind(covars,del)

## Duplication
dup = cross %>% filter(groups == "Duplication")
covars = dup %>% dplyr::select(eid,sex,site,age,groups)
dup = dup %>% dplyr::select(-all_of(names(covars))) %>% dplyr::select(-X)
dup = t(as.matrix(dup))
dup = t(ComBat(dup,batch=as.numeric(factor(covars$site)),mod=model.matrix(~age+as.numeric(sex), data = covars)))
dup = cbind(covars,dup)

## Non-carriers
dMRI = dMRI_cross[!dMRI_cross$eid %in% patho_CNV,]
covars = dMRI %>% dplyr::select(eid,sex,site,age)
dMRI = dMRI %>% dplyr::select(-all_of(names(covars))) %>% dplyr::select(-X)
dMRI = t(as.matrix(dMRI))
dMRI = t(ComBat(dMRI,batch=as.numeric(factor(covars$site)),mod=model.matrix(~age+as.numeric(sex), data = covars)))
dMRI = cbind(covars,dMRI)
dMRI$groups = "HC"
cross = rbind(dMRI, dup, del)
#dMRI_cross[!dMRI_cross$eid %in% dMRI$eid,]

# save
write.csv(cross, file = "/cluster/projects/p33/users/maxk/UKB/CNV/data/cross.csv")
# unfortunately not enough duplication and deletion carriers in the longitudinal data (N=8 each)
#TP1 = merge(cnv,dMRI_T1, by = "eid")
#TP2 = merge(cnv,dMRI_T2, by = "eid")
