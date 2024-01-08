# Loading Data ############################
library(batchelor)
library(Spectre)

setwd("/Users/ekarimi/PycharmProjects/FACS_benchmark/clean code and data/")

identity <- read.csv("/Users/ekarimi/PycharmProjects/FACS_benchmark/clean code and data/data/identity_6_CellTypes.csv")
myeloid <- read.csv("/Users/ekarimi/PycharmProjects/FACS_benchmark/clean code and data/data/myeloid_6_CellTypes.csv")

# remove ssc and fsc features
id_markers <- colnames(identity)[3:25]
my_markers <- colnames(myeloid)[3:25]
identity <- identity[3:27]
myeloid <- myeloid[3:27]

# Selecting shared markers
mut_markers <- intersect(id_markers,my_markers)
mut_markers # 7 shared features

################### Entropy of batch mixing #############
BatchEntropy <- function(dataset, batch0, L=100, M=100, k=100) {
  #entropy of batch mixing
  # L is the number bootstrapping times
  # M is the number of randomly picked cells    
  # k is the number of nearest neighbours of cell (from all batches) to check   
  
  require(RANN)  
  nbatches<-length(unique(batch0))
  
  entropy<-matrix(0,L,1)
  set.seed(0) 
  for (boot in 1:L) {
    bootsamples<-sample(1:nrow(dataset),M)
    W21<-nn2(dataset,query=dataset[bootsamples,],k)
    
    for (i in 1:length(bootsamples)){
      
      for (j in 1:nbatches) {
        xi<-max(1,sum(batch0[W21$nn.idx[i,]]==j))
        entropy[boot]<-entropy[boot]+xi*log(xi)
      }
    }
  }
  
  return( (-1)*entropy/length(bootsamples) )
}


################## Batch mixing entropy on PCAs
source('BatchMixingEntropy.R')

entropy.unc<-BatchEntropy(pca.unc$x[,1:2],batch)
entropy.mnn<-BatchEntropy(pca.mnn$x[,1:2],batch)
entropy.lm<-BatchEntropy(pca.lm$x[,1:2],batch)
entropy.combat<-BatchEntropy(pca.combat$x[,1:2],batch)

entropies<-cbind(entropy.unc,entropy.mnn,entropy.lm,entropy.combat)

png(file="results/entropy_batches_pcaspace.png",width=900,height=700) #sils_alltypes_fullspace.png
par(mfrow=c(1,1),mar=c(8,8,5,3),cex.axis=3,cex.main=2,cex.lab=3)
boxplot(entropies,main="",names=c("Uncorrected","MNN","limma","ComBat"),lwd=4,ylab="Batch mixing entropy")#,col="Yellow",ylab="Alpha dists")
dev.off()


## MNN method ######################################

mnn_identity <- as.data.frame(identity)
mnn_identity <- mnn_identity[,c(mut_markers,'Population')]
mnn_identity$batch <- rep('identity',nrow(mnn_identity))
mnn_identity

mnn_myeloid <- as.data.frame(myeloid)
mnn_myeloid <- mnn_myeloid[,c(mut_markers,'Population')]
mnn_myeloid$batch <- rep('myeloid',nrow(mnn_myeloid))
mnn_myeloid

names(mnn_identity[,c(1:9)])
identity1 <- mnn_identity[,c(1:7)]
myeloid1 <- mnn_myeloid[,c(1:7)]

mnn.out <- mnnCorrect(t(identity1),t(myeloid1))
output <- mnn.out@assays@data$corrected %>% t() %>% as.data.frame
colnames(output) <- paste(colnames(output),'MNN_corrected', sep = '_')
output$Population <- c(identity$Population,myeloid$Population)
output$Batch <- c(identity$batch,myeloid$batch)


write.files(output, file.prefix = 'data/imputed_MNN/mixed_merged', write.fcs = FALSE, write.csv = TRUE)

## CyCombine method ######################################
library(cyCombine)
library(tidyverse)


cc_identity <- identity[1:23]
cc_myeloid <- myeloid[1:23]


# Define the overlap (16 markers)
overlap_12 <- intersect(get_markers(cc_identity), get_markers(cc_myeloid))

# Define markers unique to each panel
missing_1 <- get_markers(cc_myeloid)[!(get_markers(cc_myeloid) %in%
                                         overlap_12)]
missing_2 <- get_markers(cc_identity)[!(get_markers(cc_identity) %in%
                                          overlap_12)]

# Perform imputations (and measure runtime)
start_time_cC <- Sys.time()
panel_12 <- impute_across_panels(dataset1 = cc_identity, dataset2 = cc_myeloid,
                                 overlap_channels = overlap_12, impute_channels1 = missing_1,
                                 impute_channels2 = missing_2)
end_time_cC <- Sys.time()

cc_imputed_identity <- panel_12$dataset1
cc_imputed_myeloid <- panel_12$dataset2
cc_imputed_identity["Population"] <- identity$Population
cc_imputed_identity["batch"] <- identity$batch

cc_imputed_myeloid["Population"] <- myeloid$Population
cc_imputed_myeloid["batch"] <- myeloid$batch

write.files(cc_imputed_identity, file.prefix = 'data/imputed_CyCombine/identity_6_CellTypes', write.fcs = FALSE, write.csv = TRUE)
write.files(cc_imputed_myeloid, file.prefix = 'data/imputed_CyCombine/myeloid_6_CellTypes', write.fcs = FALSE, write.csv = TRUE)



