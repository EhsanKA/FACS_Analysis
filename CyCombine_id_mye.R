library("batchelor")
library(Spectre)
library(batchelor)

label_sample <- function(df, cell_type, size= 70){
  df2 <- df[df$Population==cell_type, ]
  return(df2[sample(1:nrow(df2), size=size), ])
}

setwd("/Users/ekarimi/PycharmProjects/FACS_benchmark/")

identity <- read.csv("01022022_Focal_Identity_Sample_for_Ehsan.csv")
myeloid <- read.csv("01022022_Focal_Myeloid_Sample_for_Ehsan.csv")

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
#  sample_n(200)

#myeloid <- myeloid %>%
#  group_by(Population) %>%
#  sample_n(200)



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
write.files(identity, file.prefix = 'identity', write.fcs = FALSE, write.csv = TRUE)
write.files(myeloid, file.prefix = 'myeloid', write.fcs = FALSE, write.csv = TRUE)

#####TEST###### save both of the datasets 
#write.files(identity, file.prefix = 'identity_test', write.fcs = FALSE, write.csv = TRUE)
#write.files(myeloid, file.prefix = 'myeloid_test', write.fcs = FALSE, write.csv = TRUE)




table(identity$Population)
table(myeloid$Population)

id_markers <- colnames(identity)[59:83]

# Run DR
my_markers <- colnames(myeloid)[59:83]


# Selecting shared markers
mut_markers <- intersect(id_markers,my_markers)
mut_markers # 9 shared features

identity_sf <- identity
#identity_sf <- identity_sf[,-c('UMAP_X','UMAP_Y')]

#identity_sf <- run.umap(identity_sf, use.cols = mut_markers)
#identity_sf


myeloid_sf <- myeloid
#myeloid_sf <- myeloid_sf[,-c('UMAP_X','UMAP_Y')]

#myeloid_sf <- run.umap(myeloid_sf, use.cols = mut_markers)
#myeloid_sf



# Aligning two panels using MNN
#identity <- read.csv("/fast/home/e/ekarimi/Modi/data_transfer/identity_ready_for_RMNN.csv")
#identity$batch <- rep('identity',nrow(identity))

#myeloid <- read.csv("/fast/home/e/ekarimi/Modi/data_transfer/myeloid_ready_for_RMNN.csv")
#myeloid$batch <- rep('myeloid',nrow(myeloid))

identity <- as.data.frame(identity)
identity <- identity[,c(id_markers,'Population')]
identity$batch <- rep('identity',nrow(identity))
identity

myeloid <- as.data.frame(myeloid)
#myeloid <- myeloid[,c(mut_markers,'Population')]
myeloid <- myeloid[,c(my_markers,'Population')]
myeloid$batch <- rep('myeloid',nrow(myeloid))
myeloid

# matrices with features corresponding to rows and cells corresponding to columns
names(identity[,c(1:9)])

library(cyCombine)
library(tidyverse)
library(Spectre)


# Set data directory
#data_dir <- "/Users/ekarimi/PycharmProjects/FACS_imputation/FR-FCM-ZYV2"

# Load five panels' data - here we assume no substantial
# batch effects between the five samples
#panel_A <- prepare_data(data_dir = data_dir, pattern = "HEA_profile#A_pregated.fcs",
#                        down_sample = F, batch_ids = "A", sample_ids = 1)
#panel_B <- prepare_data(data_dir = data_dir, pattern = "HEA_profile#B_pregated.fcs",
#                        down_sample = F, batch_ids = "B", sample_ids = 1)

#panel_A <- run.umap(panel_A, get_markers(panel_A))
#plot1 <- ggplot(panel_A, aes(x=UMAP_X, y=UMAP_Y,)) + ggtitle('Identity - all features')

# Define the overlap (16 markers)

panel_A <- subset(identity, select=-c(Population, batch))
panel_B <- subset(myeloid, select=-c(Population, batch))
#panel_A <- identity[,-c(Population, batch)]
#panel_B <- myeloid[,-c("Population", "batch")]
overlap_AB <- intersect(get_markers(panel_A), get_markers(panel_B))

# Define markers unique to each panel
missing_A <- get_markers(panel_B)[!(get_markers(panel_B) %in%
                                      overlap_AB)]
missing_B <- get_markers(panel_A)[!(get_markers(panel_A) %in%
                                      overlap_AB)]

# Perform imputations (and measure runtime)
start_time_cC <- Sys.time()
panel_AB <- impute_across_panels(dataset1 = panel_A, dataset2 = panel_B,
                                 overlap_channels = overlap_AB, impute_channels1 = missing_A,
                                 impute_channels2 = missing_B)
end_time_cC <- Sys.time()


panel_A_imp <- panel_AB$dataset1
panel_B_imp <- panel_AB$dataset2

identity$Population

panel_A_imp["Population"] <- identity$Population
panel_B_imp["Population"] <- myeloid$Population


write.files(panel_A_imp, file.prefix = 'identity_cycombine_imputed', write.fcs = FALSE, write.csv = TRUE)
write.files(panel_B_imp, file.prefix = 'myeloid_cycombine_imputed', write.fcs = FALSE, write.csv = TRUE)

# Combine the dataframes for each set
panel_AB <- rbind(panel_AB$dataset1, panel_AB$dataset2)



#### SIMULATION DID NOT WORK!!!! THERE WAS AN ERROR IN THE IMPUTE ACROSS PANELS FUNCTION!!! DK :(  ################


panel_1 <- read.csv("/Users/ekarimi/PycharmProjects/Modi/simulations/sim_cycombine/df1.csv")
panel_2 <- read.csv("/Users/ekarimi/PycharmProjects/Modi/simulations/sim_cycombine/df2.csv")

panel_1 <- subset(panel_1, select = -c(X) )
panel_2 <- subset(panel_2, select = -c(X) )






# Define the overlap (16 markers)
overlap_12 <- intersect(get_markers(panel_1), get_markers(panel_2))

# Define markers unique to each panel
missing_1 <- get_markers(panel_2)[!(get_markers(panel_2) %in%
                                      overlap_12)]
missing_2 <- get_markers(panel_1)[!(get_markers(panel_1) %in%
                                      overlap_12)]

# Perform imputations (and measure runtime)
start_time_cC <- Sys.time()
panel_12 <- impute_across_panels(dataset1 = panel_1, dataset2 = panel_2,
                                 overlap_channels = overlap_12, impute_channels1 = missing_1,
                                 impute_channels2 = missing_2)
end_time_cC <- Sys.time()




######## original code of the cycombine R markdown tutorial about merging different panels   ##################


# based on https://biosurf.org/cyCombine_panel_merging.html
# link to the dataset http://flowrepository.org/id/FR-FCM-ZYV2
library(cyCombine)
library(tidyverse)
library(Spectre)


# Set data directory
data_dir <- "/Users/ekarimi/PycharmProjects/FACS_imputation/FR-FCM-ZYV2"

# Load five panels' data - here we assume no substantial
# batch effects between the five samples
panel_A <- prepare_data(data_dir = data_dir, pattern = "HEA_profile#A_pregated.fcs",
                        down_sample = F, batch_ids = "A", sample_ids = 1)
panel_B <- prepare_data(data_dir = data_dir, pattern = "HEA_profile#B_pregated.fcs",
                        down_sample = F, batch_ids = "B", sample_ids = 1)
panel_C <- prepare_data(data_dir = data_dir, pattern = "HEA_profile#C_pregated.fcs",
                        down_sample = F, batch_ids = "C", sample_ids = 1)
panel_D <- prepare_data(data_dir = data_dir, pattern = "HEA_profile#D_pregated.fcs",
                        down_sample = F, batch_ids = "D", sample_ids = 1)
panel_E <- prepare_data(data_dir = data_dir, pattern = "HEA_profile#E_pregated.fcs",
                        down_sample = F, batch_ids = "E", sample_ids = 1)


#panel_A <- run.umap(panel_A, get_markers(panel_A))
#plot1 <- ggplot(panel_A, aes(x=UMAP_X, y=UMAP_Y,)) + ggtitle('Identity - all features')

# Define the overlap (16 markers)
overlap_AB <- intersect(get_markers(panel_A), get_markers(panel_B))

# Define markers unique to each panel
missing_A <- get_markers(panel_B)[!(get_markers(panel_B) %in%
                                      overlap_AB)]
missing_B <- get_markers(panel_A)[!(get_markers(panel_A) %in%
                                      overlap_AB)]

# Perform imputations (and measure runtime)
start_time_cC <- Sys.time()
panel_AB <- impute_across_panels(dataset1 = panel_A, dataset2 = panel_B,
                                 overlap_channels = overlap_AB, impute_channels1 = missing_A,
                                 impute_channels2 = missing_B)
end_time_cC <- Sys.time()

# Combine the dataframes for each set
panel_AB <- rbind(panel_AB$dataset1, panel_AB$dataset2)


# Define the overlap (16 markers)
overlap_CD <- intersect(get_markers(panel_C), get_markers(panel_D))

# Define markers unique to each panel
missing_C <- get_markers(panel_D)[!(get_markers(panel_D) %in%
                                      overlap_CD)]
missing_D <- get_markers(panel_C)[!(get_markers(panel_C) %in%
                                      overlap_CD)]

# Perform imputations
panel_CD <- impute_across_panels(dataset1 = panel_C, dataset2 = panel_D,
                                 overlap_channels = overlap_CD, impute_channels1 = missing_C,
                                 impute_channels2 = missing_D)


# Combine the dataframes for each set
panel_CD <- rbind(panel_CD$dataset1, panel_CD$dataset2)



# And let's look at the marker density of A+B and C+D
# relative to the full panel (E)
plot_density(uncorrected = panel_E, corrected = rbind(panel_AB[,
                                                               colnames(panel_E)], panel_CD[, colnames(panel_E)]), dataset_names = c("Panel E",
                                                                                                                                     "Panel A+B and C+D"), ncol = 5)
