## input data: preprocessed by survival_mRNA_cBio.py
## /Users/yiwenbu/PycharmProjects/survival/src/survival_mRNA_cBio.py
## output data: KM plot, divided by median value.

#### SET WORKING DIRECTORY
setwd("/Users/yiwenbu/PycharmProjects/TCGA_survival_plot/results")
library(dplyr)
library(survival)


#### READ DATA
inputfile = "/Users/yiwenbu/PycharmProjects/TCGA_survival_plot/results/CMTM8_survival_leukemia.xls"
data <- read.csv(file=inputfile, header=TRUE, sep="\t")
head(data)
gene = colnames(data)[ncol(data)]
gene_median = median(data[,gene])

# add one column as label: gene high and gene low by median expression
data2 = data %>% mutate(gene_level = case_when(
  eval(as.name(gene)) >= gene_median ~ "high", eval(as.name(gene)) < gene_median ~ "low"))

# Overall Survival : fit and plot
OS_surv = Surv(data2$OS_MONTHS, data2$OS_STATUS)
OS_fit = survfit(OS_surv~data2$gene_level)
plot(OS_fit, col=c("red", "blue"), mark=3, xlab="Time(Days)", ylab="Survival Probability", main="Overall Survival")
legend("topright", c(paste0(gene," High"), paste0(gene," Low")), col=c("red","blue"), lty=1)
# DFS survival: fit and plot
DFS_surv = Surv(data2$DFS_MONTHS, data2$DFS_STATUS)
DFS_fit = survfit(DFS_surv~data2$gene_level)
plot(DFS_fit, col=c("red", "blue"), mark=3, xlab="Time(Days)", ylab="Survival Probability", main="Disease Free Survival")
legend("topright", c(paste0(gene," High"), paste0(gene," Low")), col=c("red","blue"), lty=1)

# Statistics
OS_diff = survdiff(OS_surv~data2$gene_level)
OS_p_value = 1-pchisq(OS_diff$chisq,df=1)       # p = 0.02
print (paste0("Overall Survival: ", OS_p_value))
DFS_diff = survdiff(DFS_surv~data2$gene_level)
DFS_p_value = 1-pchisq(DFS_diff$chisq,df=1)     # p = 0.06
print (paste0("Disease Free Survival: ", DFS_p_value))
details = table(data2$gene_level)
details


