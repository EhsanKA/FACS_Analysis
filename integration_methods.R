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




# cc -> cycombine
# cbb -> cytobackbone
# ctm -> cytofmerge

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




## CytoBackBone method ##################################################################
library(CytoBackBone, exclude = "plot")  # v. 1.0.0, CytoBackBone has a function called 'plot', which masks the base::plot, but we overwrite that behavior
library(flowCore)

DataList=list()
cbb_identity <- subset(identity, select=-c(Population, batch))
cbb_myeloid <- subset(myeloid, select=-c(Population, batch))

DataList[["identity"]] <- cbb_identity
DataList[["myeloid"]] <- cbb_myeloid
### USING THE GITHUB CODE FROM SUBARNA TO CONVERT THE CSV FILE TO FSC
AllSampleNames <- names(DataList)

## Chech data quality
head(DataList)

##### END USER INPUT #
library(data.table)
x <- Sys.time()
x <- gsub(":", "-", x)
x <- gsub(" ", "_", x)

FCSList=list()
for(i in c(1:length(AllSampleNames))){
  data_subset <- DataList[i]
  data_subset <- rbindlist(as.list(data_subset))
  dim(data_subset)
  a <- names(DataList)[i]
  
  metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=dimnames(data_subset)[[2]])
  
  ## Create FCS file metadata - ranges, min, and max settings
  #metadata$range <- apply(apply(data_subset,2,range),2,diff)
  #metadata$minRange <- apply(data_subset,2,min)
  #metadata$maxRange <- apply(data_subset,2,max)
  
  data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
  head(data_subset.ff)
  FCSList[[a]] <- data_subset.ff
  write.FCS(data_subset.ff, paste0(getwd(), "/data/imputed_CytoBackBone/", a, ".fcs"))
}


panel_A <- FCSList[['identity']]
panel_B <- FCSList[['myeloid']]

panel_A_CBB <-import.FCS(paste0(getwd(), "/data/imputed_CytoBackBone/identity.fcs"), trans = "none")
panel_B_CBB <- import.FCS(paste0(getwd(), "/data/imputed_CytoBackBone/myeloid.fcs"), trans = "none")


overlap_AB <- intersect(colnames(panel_A), colnames(panel_B))

# Define markers unique to each panel
missing_A <- colnames(panel_B)[!(colnames(panel_B) %in%
                                   overlap_AB)]
missing_B <- colnames(panel_A)[!(colnames(panel_A) %in%
                                   overlap_AB)]

start_time_CBB <- Sys.time()
panel_AB_CBB <- CytoBackBone::merge(FCS1 = panel_A_CBB, FCS2 = panel_B_CBB,
                                    BBmarkers = overlap_AB, th = 3,
                                    leftout = T, normalize = FALSE)

end_time_CBB <- Sys.time()

cbb_out <- data.frame(panel_AB_CBB$merged@intensities)
colnames(cbb_out) <- panel_AB_CBB$merged@markers

write.files(cbb_out, file.prefix = 'data/imputed_CytoBackBone/mixed_merged', write.fcs = FALSE, write.csv = TRUE)


## CyTOFMerge method ##################################################################


devtools::source_url("https://github.com/tabdelaal/CyTOFmerge/blob/master/CombineFCS.R?raw=TRUE")



DataList=list()
cm_identity <- subset(identity, select=-c(Population, batch))
cm_myeloid <- subset(myeloid, select=-c(Population, batch))

DataList[["identity"]] <- cm_identity
DataList[["myeloid"]] <- cm_myeloid
### USING THE GITHUB CODE FROM SUBARNA TO CONVERT THE CSV FILE TO FSC
AllSampleNames <- names(DataList)

## Chech data quality
head(DataList)

##### END USER INPUT #
library(data.table)
x <- Sys.time()
x <- gsub(":", "-", x)
x <- gsub(" ", "_", x)

FCSList=list()
for(i in c(1:length(AllSampleNames))){
  data_subset <- DataList[i]
  data_subset <- rbindlist(as.list(data_subset))
  dim(data_subset)
  a <- names(DataList)[i]
  
  metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=dimnames(data_subset)[[2]])
  
  ## Create FCS file metadata - ranges, min, and max settings
  #metadata$range <- apply(apply(data_subset,2,range),2,diff)
  #metadata$minRange <- apply(data_subset,2,min)
  #metadata$maxRange <- apply(data_subset,2,max)
  
  data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
  head(data_subset.ff)
  FCSList[[a]] <- data_subset.ff
  write.FCS(data_subset.ff, paste0(getwd(), "/data/imputed_CyTOFMerge/", a, ".fcs"))
}

panel_A <- FCSList[['identity']]
panel_B <- FCSList[['myeloid']]

panel_A_cm <- import.FCS(paste0(getwd(), "/data/imputed_CyTOFMerge/identity.fcs"), trans = "none")
panel_B_cm <- import.FCS(paste0(getwd(), "/data/imputed_CyTOFMerge/myeloid.fcs"), trans = "none")


overlap_AB <- intersect(colnames(panel_A), colnames(panel_B))

# Define markers unique to each panel
missing_A <- colnames(panel_B)[!(colnames(panel_B) %in%
                                   overlap_AB)]
missing_B <- colnames(panel_A)[!(colnames(panel_A) %in%
                                   overlap_AB)]

start_time_cm <- Sys.time()
panel_AB_cm <- CombineFCS(FCSfile1 = paste0(getwd(), "/data/imputed_CyTOFMerge/identity.fcs"),RelevantMarkers1 = 1:23,
                          FCSfile2 = paste0(getwd(), "/data/imputed_CyTOFMerge/myeloid.fcs"), RelevantMarkers2 = 1:23, 
                          arcsinhTrans = F)

end_time_cm <- Sys.time()

write.files(panel_AB_cm, file.prefix = 'data/imputed_CyTOFMerge/mixed_merged', write.fcs = FALSE, write.csv = TRUE)

