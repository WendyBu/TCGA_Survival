## loop to look for the best seperation of survival KM plot
## ## input data: preprocessed by survival_mRNA_cBio.py
## /Users/yiwenbu/PycharmProjects/survival/src/survival_mRNA_cBio.py

#### SET WORKING DIRECTORY
#setwd("/Users/yiwenbu/PycharmProjects/survival/results")
#setwd("~/Desktop/P12_Shu_hxia/SURVIVAL")
setwd("/Users/yiwenbu/PycharmProjects/TCGA_survival_plot/results")
library(dplyr)
library(survival)


#### READ DATA, inputFile
inputFile = '/Users/yiwenbu/PycharmProjects/TCGA_survival_plot/results/CMTM8_survival_leukemia.xls'
data <- read.csv(file=inputFile, header=TRUE, sep="\t")
head(data)  # data format
"""
                X DFS_MONTHS DFS_STATUS OS_MONTHS OS_STATUS   LILRB3
1 TCGA-AB-2932-03       42.1          0      42.1         0 5352.921
2 TCGA-AB-2987-03       59.0          0      59.0         0 5263.821
3 TCGA-AB-2987-03        5.7          1       8.8         1 5263.821
4 TCGA-AB-2987-03       32.7          1      32.7         0 5263.821
5 TCGA-AB-2987-03        8.3          1      10.7         1 5263.821
6 TCGA-AB-2987-03       12.0          1      20.5         1 5263.821
"""
gene = colnames(data)[ncol(data)]  # LILRB3


## functions:
# 1. labelData: label data with high or low, based on top cutoff or low cutoff 
# 2. survival_pvalue: calculate p value, giving dataframe and survival type(OS or DFS)
# 3. plot_surval: plot survival, giving labeled data, survival_type (OS or DFS)

################## label high or low  based on geneName, cutoff_top, cutoff_bottom
labelData <- function(dataset, gene, cutt_top, cutt_bottom){
  labeled_data = cbind(dataset)
  labeled_data$label=NA
  for (row in 1:nrow(labeled_data)){
    if (labeled_data[row, gene] >= cutt_top){
      labeled_data[row,"label"]='high'
    } else if (labeled_data[row, gene] < cutt_bottom){
      labeled_data[row,"label"]='low'
    } 
  }
  newDataSet = labeled_data[complete.cases(labeled_data),]
  return(newDataSet)
}


##################   calculate p value 
survival_pvalue <- function(data, survival_type){
  if (survival_type == "OS"){
    mySurv<-Surv(time=data$OS_MONTHS, event = data$OS_STATUS)
  }else
    if (survival_type == "DFS"){
    mySurv<-Surv(time=data$DFS_MONTHS, event = data$DFS_STATUS)
    }
  myfit<-survfit(mySurv~data$label)
  mydiff = survdiff(mySurv~data$label)
  myPvalue = pchisq(mydiff$chisq, length(mydiff$n)-1, lower.tail = FALSE)
  return (myPvalue)
}


###################### plot survival plot. return sample numbers and p value
plot_surval <- function(Fdata, survival_type){
  if (survival_type == "DFS"){
  mySurv<-Surv(time=Fdata$DFS_MONTHS, event = Fdata$DFS_STATUS)
  myfit<-survfit(mySurv~Fdata$label)
  plot(myfit, col=c("red", "blue"), mark=3, xlab="Time (Months)", ylab="Survival Probability", main="Disease Free Survival")
  } else
    if (survival_type == "OS"){
    mySurv<-Surv(time=Fdata$OS_MONTHS, event = Fdata$OS_STATUS)
    myfit<-survfit(mySurv~Fdata$label)
    plot(myfit, col=c("red", "blue"), mark=3, xlab="Time (Months)", ylab="Survival Probability", main="Overall Survival")
    }
  legend("topright", c(paste0(gene," High"), paste0(gene," Low")), col=c("red","blue"), lty=1)
  detail = table(Fdata$label)
  return (detail)
}


##### main #####
# 1. look for cutoff by looping, could adjust by = , each step of the change 
# 2. call the survival_pvalue, record p value< 0.05 data. 
# 3. chose the p <0.05 and sample numbers are the largest cutoff automatically
# 4. using the corresponding cutoff to plot survival. 
# 5. alternatively, manually chose the pvalue, cutoffs, and plot the survival curve. 
##################  LOOKING FOR CUTTOFF
survival_type = "DFS"
pvalue_tables = data.frame()  
bottom = quantile(data[,gene], probs = seq(0.50, 0.01, by= -0.03))  # could change by=, loop steps
for (n in bottom) {
  cutt_bottom = n
  top = quantile(data[,gene], probs = seq(0.51, 0.99, by= 0.03))
  for (j in top) {
    cutt_top = j
    newdata = labelData(data, gene, cutt_top, cutt_bottom)
    if (length(unique(newdata$label)) < 2){
      next
    }
    p = survival_pvalue(newdata, survival_type)
    if (p < 0.05){
      num_high_label = sum(newdata$label == "high")
      num_low_label = sum(newdata$label == "low")
      values <- c(cutt_top, cutt_bottom, p, num_high_label, num_low_label)
      pvalue_tables <- rbind(pvalue_tables, values)   # generate all p<0.05 top and bottom
    }
  }
}
## add colname to the pvalue_table
x <- c("top", "bottom", "pvalue", "num_high_samples", "num_low_samples")
colnames(pvalue_tables) <- x
print (pvalue_tables)
## find the pvalue <0.05 and the sample numbers are the biggest
diff_list = c()
for (i in (1:nrow(pvalue_tables))){
  diff = pvalue_tables[i,'top'] - pvalue_tables[i, "bottom"]
  diff_list[i] = diff
}
minpos = which.min(diff_list[1:length(diff_list)])
finalTop = pvalue_tables[minpos, 'top']
finalBottom = pvalue_tables[minpos, 'bottom']
print (paste0("Final p value: ", pvalue_tables[minpos, 'pvalue']))
finalset = labelData(data, gene, finalTop, finalBottom)  # label final dataset
details = plot_surval(finalset, survival_type)   # plot final survival graph
print ("Sample Numbers:")
print (details)


### call specific top and bottom for survival curve
# survival_type = "DFS"
# finalTop = 689.1010
# finalBottom = 81.25309
# finalset = labelData(data, gene, finalTop, finalBottom)
# details = plot_surval(finalset, survival_type)
# print ("Sample Numbers:")
# print (details)



