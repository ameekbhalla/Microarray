#both sexes, quantile normalized across 
group2 <- with(pData(eset_mini), 
               paste(replicate, 
                     day,
                     sep = "."))
pData(eset_mini)$group2 <- as.factor(group2)
corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$group2)
corfit_1$consensus
[1] 0.6202567 #12 warnings


group2 <- with(pData(eset_mini), 
               paste(sex, 
                     replicate,
                     sep = "."))
pData(eset_mini)$group2 <- as.factor(group2)
corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$group2)
corfit_1$consensus
[1] 0.6202567 #23 warnings


group2 <- with(pData(eset_mini), 
               paste(sex, 
                     day,
                     sep = "."))
pData(eset_mini)$group2 <- as.factor(group2)
corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$group2)
corfit_1$consensus
[1] 0.4487796


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$day)
corfit_1$consensus
[1] 0.4487796


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$sex)
corfit_1$consensus
[1] 0.4763534


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$replicate)
corfit_1$consensus
[1] -0.005883273




#only females, quantile normalized across 
group2 <- with(pData(eset_mini), 
               paste(replicate, 
                     day, 
                     sep = "."))
pData(eset_mini)$group2 <- as.factor(group2)
corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$group2)
corfit_1$consensus
[1] 0.3365879


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$replicate)
corfit_1$consensus
[1] 0.3365862
 #warning message for one gene

corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$day)
corfit_1$consensus
[1] 0.131445 #warning message for three genes




####both sexes, q-nomralized with only f0 data
group2 <- with(pData(eset),
               paste(replicate, 
                     day, 
                     sep = "."))
pData(eset)$group2 <- as.factor(group2)
corfit_1 <- duplicateCorrelation(eset, design_1, block =  pData(eset)$group2)
corfit_1$consensus
[1] 0.4771085


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$replicate)
corfit_1$consensus
[1] -0.08968102

corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$day)
corfit_1$consensus
[1] 0.4771085


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$sex)
corfit_1$consensus
[1] 0.4119819




####only females, q-nomralized with only f0 data
group2 <- with(pData(eset_mini), 
               paste(replicate, 
                     day,
                     sep = "."))
pData(eset_mini)$group2 <- as.factor(group2)
corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$group2)
corfit_1$consensus
[1] 0.3263993


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$replicate)
corfit_1$consensus
[1] 0.3263993 #warning message for one gene


corfit_1 <- duplicateCorrelation(eset_mini, design_1, block =  pData(eset_mini)$day)
corfit_1$consensus
[1] 0.3263993 #warning message for three genes