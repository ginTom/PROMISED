library(biomaRt)
library(hablar)
library(readr)
library(tidyverse)
library(stringr)
library(openxlsx)
library(mixOmics)
library(janitor)
library(Biobase)
library(genefilter)
library(ggplot2)
library(gridExtra)
library(kableExtra)


rm(list = ls())
setwd("C:/Users/tom_m/Desktop/PhD/DGF-Kidney/MOGONET/DGF")

prep.fun = function(x){
  table(unlist(sapply(x,class)))
}

transpose.df = function(x, first_colnames = TRUE){
  x_t <- as.data.frame(t(x))
  if (first_colnames){
    colnames(x_t) <- x_t[1,]
    x_t <- x_t[-1,]
  }
  return(x_t)
}

mirna_data_train <- read.csv(file = "1_tr.csv", header = FALSE, sep=",")
mirna_features <- read.csv(file = "1_featname.csv", header = FALSE, sep=",")$V1
colnames(mirna_data_train) <- mirna_features

snp_data_train <- read.csv(file = "2_tr.csv", header = FALSE, sep=",")
snp_features <- read.csv(file = "2_featname.csv", header = FALSE, sep=",")$V1
colnames(snp_data_train) <- snp_features

label_train = read.csv(file = "labels_tr.csv", header = FALSE, sep=",")$V1
label_train[label_train==0] = "DGF"
label_train[label_train==1] = "EGF"
label_train = as.factor(label_train)
summary(label_train)

data = list(miRNA = mirna_data_train, 
            snp = snp_data_train)

design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

basic.diablo.model = block.splsda(X = data, Y = label_train, ncomp = 5, design = design) 

perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 10, nrepeat = 10)  

plot(perf.diablo)

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"]
test.keepX = list (miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)), 
                   snp = c(5:9, seq(10, 18, 2), seq(20,30,5)))

tune.TCGA = tune.block.splsda(X = data, Y = label_train, ncomp = 3, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              dist = "mahalanobis.dist")
list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
list.keepX
final.diablo.model = block.splsda(X = data, Y = label_train, ncomp = 3, 
                                  keepX = list.keepX, design = design)

final.diablo.model$design # design matrix for the final model

sig_mirna <- selectVar(final.diablo.model, block = 'miRNA', comp = 1)
sig_mirna <- as.data.frame(sig_mirna[["miRNA"]][["value"]])

ord <- as.data.frame(sig_mirna[order(-sig_mirna$value.var),])
rownames(ord) <- rownames(sig_mirna)[order(-sig_mirna$value.var)]
write.csv(file="mirna_biomarker_mixomics.csv", ord, row.names = TRUE)

sig_snp <- selectVar(final.diablo.model, block = 'snp', comp = 1)
sig_snp <- as.data.frame(sig_snp[["snp"]][["value"]])
ord <- as.data.frame(sig_snp[order(-sig_snp$value.var),])
rownames(ord) <- rownames(sig_snp)[order(-sig_snp$value.var)]
# write.csv(file="snp_biomarker_mixomics.csv", ord, row.names = TRUE)

auc.splsda = auroc(final.diablo.model, roc.block = "miRNA", 
                   roc.comp = 2, print = TRUE)

mirna_data_test <- read.csv(file = "1_te.csv", header = FALSE, sep=",")
colnames(mirna_data_test) <- mirna_features

snp_data_test <- read.csv(file = "2_te.csv", header = FALSE, sep=",")
colnames(snp_data_test) <- snp_features

label_test = read.csv(file = "labels_te.csv", header = FALSE, sep=",")$V1
label_test[label_test==0] = "DGF"
label_test[label_test==1] = "EGF"
label_test = as.factor(label_test)
summary(label_test)

data_test = list(miRNA = mirna_data_test, 
            snp = snp_data_test)
predict.diablo = predict(final.diablo.model, newdata = data_test)
confusion.mat = get.confusion_matrix(truth = label_test,
                                     predicted = predict.diablo$WeightedVote$centroids.dist[,2])
(confusion.mat[1,1]+confusion.mat[2,2])/sum(confusion.mat)
