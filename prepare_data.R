library("batchelor")
library(Spectre)
?do.subsample
library(batchelor)
library(ggplot2)

label_sample <- function(df, cell_type, size= 70){
  df2 <- df[df$Population==cell_type, ]
  return(df2[sample(1:nrow(df2), size=size), ])
}



############################################
#### SIMULATION DID NOT WORK!!!! THERE WAS AN ERROR IN THE IMPUTE ACROSS PANELS FUNCTION!!! DK :(
#panel_1 <- read.csv("/Users/ekarimi/PycharmProjects/Modi/simulations/sim_cycombine/df1.csv")
#panel_2 <- read.csv("/Users/ekarimi/PycharmProjects/Modi/simulations/sim_cycombine/df2.csv")

#panel_1 <- subset(panel_1, select = -c(X) )
#panel_2 <- subset(panel_2, select = -c(X) )
############################################



setwd("/Users/ekarimi/PycharmProjects/FACS_benchmark/clean code and data/")

identity <- read.csv("/Users/ekarimi/PycharmProjects/FACS_benchmark/clean code and data/data/01022022_Focal_Identity_Sample_for_Ehsan.csv")
myeloid <- read.csv("/Users/ekarimi/PycharmProjects/FACS_benchmark/clean code and data/data/01022022_Focal_Myeloid_Sample_for_Ehsan.csv")

hist(myeloid$FJComp.BV510.A_efluor506_asinh_aligned, breaks = 100)
hist(myeloid$FJComp.APC.H7.A_HLADR_asinh_aligned, breaks = 100)
myeloid <- myeloid[myeloid$FJComp.BV510.A_efluor506_asinh_aligned > -3, ]

colnames(identity)[which(names(identity) == "Refined_clustering")] <- "Population"

table(identity$Population)
table(myeloid$Population)
hist(myeloid$FJComp.BV510.A_efluor506_asinh_aligned, breaks = 100)
hist(myeloid$FJComp.APC.H7.A_HLADR_asinh_aligned, breaks = 100)

### generating random samples
library(dplyr)
set.seed(12345)
#identity <- identity %>%
#  group_by(Population) %>%
#  sample_n(100)

#myeloid <- myeloid %>%
#  group_by(Population) %>%
#  sample_n(100)



myeloid["Population"][myeloid["Population"] == "cDCs1"] <- "Dendritic cells"
myeloid["Population"][myeloid["Population"] == "cDCs2"] <- "Dendritic cells"
myeloid["Population"][myeloid["Population"] == "pDCs"] <- "Dendritic cells"
myeloid["Population"][myeloid["Population"] == "Classical monocytes"] <- "Classical Monocytes"
myeloid["Population"][myeloid["Population"] == "Endothelial cells"] <- "Endothelial cells"
myeloid["Population"][myeloid["Population"] == "Granulocytes"] <- "Granulocytic cells"
myeloid["Population"][myeloid["Population"] == "Classical monocytes"] <- "Classical Monocytes"
myeloid["Population"][myeloid["Population"] == "Non-classical monocytes"] <- "Non-classical Monocytes"

identity["Population"][identity["Population"] == "Myeloblasts and Promyelocytes"] <- "HSCs/MPPs/Myeloblasts/Promyelocytes"
identity["Population"][identity["Population"] == "HSCs and MPPs"] <- "HSCs/MPPs/Myeloblasts/Promyelocytes"


table(identity$Population)
table(myeloid$Population)

identity <- identity[identity$Population != "T cell", ]
identity <- identity[identity$Population != "Pro-B cells and Pre-B cells", ]
identity <- identity[identity$Population != "Pre-Pro B cells", ]
identity <- identity[identity$Population != "Plasma cells", ]
identity <- identity[identity$Population != "Non switched memory B cells", ]
identity <- identity[identity$Population != "Mesenchymal stromal cells", ]
identity <- identity[identity$Population != "Mature naive B cells", ]
identity <- identity[identity$Population != "Immature B cells", ]
identity <- identity[identity$Population != "Erythroid progenitors", ]
identity <- identity[identity$Population != "Class switched memory B cells", ]
identity <- identity[identity$Population != "CD56dim CD16+ NK cells", ]
identity <- identity[identity$Population != "CD56bright CD16- NK cells", ]


myeloid <- myeloid[myeloid$Population != "Basophils", ]
myeloid <- myeloid[myeloid$Population != "Eosinophils", ]

##### save both of the datasets after stratification and changing labels
#write.files(identity, file.prefix = 'data/identity_6_CellTypes', write.fcs = FALSE, write.csv = TRUE)
#write.files(myeloid, file.prefix = 'data/myeloid_6_CellTypes', write.fcs = FALSE, write.csv = TRUE)

#####TEST###### save both of the datasets 
#write.files(identity, file.prefix = 'identity_test', write.fcs = FALSE, write.csv = TRUE)
#write.files(myeloid, file.prefix = 'myeloid_test', write.fcs = FALSE, write.csv = TRUE)




table(identity$Population)
table(myeloid$Population)

id_markers <- colnames(identity)[59:83]

my_markers <- colnames(myeloid)[59:83]

# Selecting shared markers
mut_markers <- intersect(id_markers,my_markers)
mut_markers # 9 shared features

identity <- as.data.frame(identity)
identity <- identity[,c(id_markers,'Population')]
identity$batch <- rep('identity',nrow(identity))

myeloid <- as.data.frame(myeloid)
#myeloid <- myeloid[,c(mut_markers,'Population')]
myeloid <- myeloid[,c(my_markers,'Population')]
myeloid$batch <- rep('myeloid',nrow(myeloid))



#### MNNNNNNNNNNNNNNNNN part to check the plots #######

# remove ssc and fsc features
id_markers <- setdiff(colnames(identity), c("FSC.A_asinh_aligned", "SSC.A_asinh_aligned", "batch", "Population"))
my_markers <- setdiff(colnames(myeloid), c("FSC.A_asinh_aligned", "SSC.A_asinh_aligned", "batch", "Population"))
identity <- identity[3:27]
myeloid <- myeloid[3:27]

# Selecting shared markers
mut_markers <- intersect(id_markers,my_markers)
mut_markers # 7 shared features




# cc -> cycombine
# cbb -> cytobackbone
# ctm -> cytofmerge
mut_marker_label <- c(mut_markers, c("batch", "Population"))
mixed_mut_data <- rbind(identity[,mut_marker_label], myeloid[,mut_marker_label])

library(umap)

u_mixed <- run.umap(mixed_mut_data, mut_markers)
#par(mfrow=c(1,2))

plot1 <- ggplot(u_mixed, aes(x=UMAP_X, y=UMAP_Y, colour=Population)) + geom_point() + ggtitle("mixed_mut_before_mnn_umap_population")
plot1
ggsave("mut_umap_before_mnn_population.png")

plot2 <- ggplot(u_mixed, aes(x=UMAP_X, y=UMAP_Y, colour=batch)) + geom_point() + ggtitle("mixed_mut_before_mnn_umap_batch")
plot2
ggsave("mut_umap_before_mnn_batch.png")


## MNN method ######################################

mnn_identity <- as.data.frame(identity)
mnn_identity <- mnn_identity[,c(mut_markers,'Population')]
mnn_identity$batch <- rep('identity',nrow(mnn_identity))
# normalisation identity
mnn_identity[,mut_markers] <- mnn_identity[,mut_markers] / 
  matrix(sqrt(colSums(mnn_identity[,mut_markers]*mnn_identity[,mut_markers])),
         nrow=nrow(mnn_identity[,mut_markers]),
         ncol=ncol(mnn_identity[,mut_markers]), byrow=TRUE)
mnn_identity

mnn_myeloid <- as.data.frame(myeloid)
mnn_myeloid <- mnn_myeloid[,c(mut_markers,'Population')]
mnn_myeloid$batch <- rep('myeloid',nrow(mnn_myeloid))
# normalisation myeloid
mnn_myeloid[,mut_markers] <- mnn_myeloid[,mut_markers] / 
  matrix(sqrt(colSums(mnn_myeloid[,mut_markers]*mnn_myeloid[,mut_markers])),
         nrow=nrow(mnn_myeloid[,mut_markers]),
         ncol=ncol(mnn_myeloid[,mut_markers]), byrow=TRUE)
mnn_myeloid

names(mnn_identity[,c(1:9)])
identity1 <- mnn_identity[,c(1:7)]
myeloid1 <- mnn_myeloid[,c(1:7)]

mnn.out <- mnnCorrect(t(identity1),t(myeloid1))
output <- mnn.out@assays@data$corrected %>% t() %>% as.data.frame
colnames(output) <- paste(colnames(output),'MNN_corrected', sep = '_')
output$Population <- c(identity$Population,myeloid$Population)
output$Batch <- c(identity$batch,myeloid$batch)

# mnn_after_identity <- output[output$Batch=="identity",]
# mnn_after_myeloid <- output[output$Batch=="myeloid",]

mnn_mut_markers <- setdiff(colnames(output), c("Population", "Batch"))

u_mnn <- run.umap(output, mnn_mut_markers)
plot3 <- ggplot(u_mnn, aes(x=UMAP_X, y=UMAP_Y, colour=Population)) + geom_point() + ggtitle("mnn_mut_umap_population")
plot3
ggsave("mnn_mut_umap_population.png")

# u_myeloid_mnn <- run.umap(mnn_after_myeloid, mnn_mut_markers)
plot4 <- ggplot(u_mnn, aes(x=UMAP_X, y=UMAP_Y, colour=Batch)) + geom_point() + ggtitle("mnn_mut_umap_batch")
plot4
ggsave("mnn_mut_umap_batch.png")





write.files(output, file.prefix = 'data/imputed_MNN/mixed_merged', write.fcs = FALSE, write.csv = TRUE)


##### END of MNNNNNNNNNNNNNN #####
# save the data

write.files(identity, file.prefix = 'data/identity_6_CellTypes', write.fcs = FALSE, write.csv = TRUE)
write.files(myeloid, file.prefix = 'data/myeloid_6_CellTypes', write.fcs = FALSE, write.csv = TRUE)

