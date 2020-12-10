# restore_packages.R
#
# installs each package from the stored list of packages

# load("installed_packages.rda")
load("//bhalla-desktop/Downloads/R Projects/Microarray/installed_packages.rda")

for (count in 1:length(installed_packages)) install.packages(installed_packages[count])


tmp = installed.packages()

new = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
bio = setdiff(installedpackages,new)

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

for (count in 1:length(bio)) BiocManager::install(bio[count])
BiocManager::install("affy")
library(BiocManager)
