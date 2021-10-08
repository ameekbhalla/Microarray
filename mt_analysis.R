## Load Packages
library(affy) # Bioconductor package for analysis of older Affymetrix arrays
#library(affyio)
#library(affyPLM)
#library(annaffy)
#library(annotate)
#library(arrayQualityMetrics)
library(BiocManager)
library(drosophila2.db) # annotation package for the model organism
library(drosophila2cdf)
library(formatR)
library(gcrma)
library(htmltools)
library(knitr)
library(limma)
#library(markdown)
library(rmarkdown)
#library(simpleaffy)
library(tidyverse)

# untar the file downloaded from GEO
untar("GSE119927_RAW.tar")

# on untarring gz files are obtained, g-unzip them
here::here() %>% ##folder path
  list.files(pattern = "gz$") %>% 
  map(R.utils::gunzip)

#Load the package into memory to get the whole data from GEO
library("GEOquery")

# read the whole data into memory as an ExpressionSet object (the GEO data is rma-ed already)
mt_whole_data <- getGEO(filename="GSE119927_series_matrix.txt.gz", getGPL = FALSE)
# sample_names_m <- which(pData(eset)$characteristics_ch1.1 == "generation: F0" & pData(eset)$characteristics_ch1.5 == "gender: male")
# eset_m <- eset[ , sample_names_m] #ExpressionSet objects cannot be subset with filter() or filter()

# get annotation information
sample_data <- pData(mt_whole_data) %>% filter(characteristics_ch1.1 == "generation: F0", characteristics_ch1.5 == "gender: male")

# read CEL files
raw_data <- ReadAffy()
phenoData(raw_data) <- AnnotatedDataFrame(sample_data)

# annotate
sampleNames(raw_data) <- sample_data$geo_accession

#check validity
validObject(raw_data)

# gcrma gives 1999 genes; rma gives 3820 DE genes; GEO gives 3244 DE genes therefore rma is closer to GEO results
eset <- gcrma(raw_data)
# eset <- rma(raw_data)

raw_expression_raw_data <- exprs(raw_data) %>% #plotDensities() accepts only matrix/ExpressionSet objects
  log2()
norm_expression_eset <- exprs(eset) %>% #plotDensities() accepts only matrix/ExpressionSet objects
  log2()

plotDensities(raw_expression_raw_data, legend = FALSE)
plotDensities(norm_expression_eset, legend = FALSE)
# plotDensities(exprs(eset_m), legend = FALSE)
# plotDensities(exprs(eset), legend = FALSE)

#####
pData(eset) <- pData(eset) %>% 
  mutate(
    diet = `diet fed:ch1` %>% as_factor())


# Create design matrix with group-means parametrization (no intercept)
design_1 <- model.matrix(~ 0 
                         # + day 
                         + diet, data = pData(eset))

colnames(design_1) <- levels(pData(eset)$diet)
colnames(design_1)

# Count the number of samples modeled by each coefficient
colSums(design_1)

#estimate the correlation between measurements made on the same repliactes:
# corfit_1 <- duplicateCorrelation(eset_f1_f, design_1, block =  pData(eset_f1_f)$replicate)
# corfit_1$consensus

# corfit_1$consensus <- by_sex$both[[1]]$consensus.correlation
# corfit_1$cor <- by_sex$both[[1]]$cor
# corfit_1$atanh.correlations <- by_sex$both[[1]]$atanh.correlations

#fit the linear model fit_model_1 with blocking
# fit_model_1 <- lmFit(eset_mini, design_1, 
#                      block =  pData(eset)$replicate, 
#                      correlation = corfit_1$consensus)

# #fit the linear model fit_model_1 without blocking
fit_model_1 <- lmFit(eset, design_1)

## Build the contrasts matrix for discovering main effect of genotype in females

cm_eset <- makeContrasts(
  f0_diet = HSD-CD,
  levels = design_1)

cm_eset

#fit the contrast matrix
fit2_model_1 <- contrasts.fit(fit_model_1, cm_eset)

#set trend = TRUE  since the data do not show a uniform variance across all expression intensities (this trend in the variance at different mean intensities is discovered with `plotSA(fit_model)` ) 
results_eset <- eBayes(fit2_model_1, 
                       trend = TRUE, #accommodates a mean-variance trend; 
                       robust = TRUE #increases power for some kinds of studies
) 

#Results of all contrasts simultaneously
decideTests(results_eset) %>% summary()

# #Venn diagram showing numbers of genes significant in all contrasts, simultaneously
decideTests(results_eset) %>% vennDiagram()

#A list of top genes differentially expressed in a contarst can be obtained by specifying its coef
# topTable(fit2_model_1, coef=c(1:3), adjust="BH", number = 30)
f0_diet <- topTable(
  results_eset,
  adjust="BH", number = nrow(fit2_model_1), sort.by = "none")


#annotation returns mutliple matches for some of the probes so perform it after obtaining DE gene list 
anno_eset <- AnnotationDbi::select(drosophila2.db,
                                   keys = (featureNames(eset)), #indexes the rows
                                   columns = c("SYMBOL", "GENENAME"), #indexes the columns
                                   keytype = "PROBEID") %>% #arranges the output by this column
  filter(!is.na(SYMBOL)) %>% # this also removes quality control probes
  group_by(PROBEID) %>%
  filter(n() == 1) # remove  probes which have multiple genes associated with them
# check that all quality control probes, containing the pattern AFFX in their name, were removed
# str_detect(anno_eset$PROBEID, "AFFX") %>% sum()

library(fdrtool)


f0_diet <- (f0_diet$t - median(f0_diet$t)) %>% 
  fdrtool(statistic = "normal", plot = F) %>% 
  cbind(f0_diet, .) %>% 
  rownames_to_column() %>% 
  left_join(anno_eset, by = c("rowname" = "PROBEID")) %>% 
  dplyr::rename(probe_id = rowname, symbol = SYMBOL, gene_name = GENENAME) %>% 
  dplyr::select("probe_id", "symbol", "gene_name", "logFC", "pval", "qval") #  %>% drop_na()

f0_diet %>% filter(qval <= 0.05) %>% count()
