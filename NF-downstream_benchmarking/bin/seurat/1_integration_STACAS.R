#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
reticulate::use_python('/usr/bin/python3.7')
library(Seurat)
library(STACAS)
library(future)
library(tidyverse)

spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')

    plot_path = "./output/NF-downstream_analysis/integration_seurat/plots/"
    rds_path = "./output/NF-downstream_analysis/integration_seurat/rds_files/"
    data_path = "./output/NF-scRNAseq_alignment/cellranger/count/filtered_feature_bc_matrix"
    ncores = 8

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')

    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/NF-scRNAseq_alignment/cellranger/count/filtered_feature_bc_matrix"
    ncores = opt$cores

    # Multi-core when running from command line
    plan("multiprocess", workers = ncores)
    options(future.globals.maxSize = 12* 1024^3) # 32gb

  } else {
    stop("--runtype must be set to 'nextflow'")
  }

  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

# Make dataframe with stage and replicate info extracted from path
input <- list.dirs(data_path, recursive = FALSE, full.names = TRUE)
input <- data.frame(sample = sub('.*/', '', input), run = str_split(sub('.*/', '', input), pattern = "_", simplify = T)[,2], path = input)

# Init list of seurat objects then merge
seurat_list <- apply(input, 1, function(x) CreateSeuratObject(counts= Read10X(data.dir = x[["path"]]), project = x[["sample"]]))
names(seurat_list) <- input$sample
seurat_all <- merge(x = seurat_list[[1]], y=seurat_list[-1], add.cell.ids = names(seurat_list), project = "chick.10x")

# Add metadata col for seq run
seurat_all@meta.data[["run"]] <- gsub(".*_", "", as.character(seurat_all@meta.data$orig.ident))
seurat_all@meta.data[["stage"]] <- gsub("_.*", "", as.character(seurat_all@meta.data$orig.ident))

# Convert metadata character cols to factors
seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)] <- lapply(seurat_all@meta.data[sapply(seurat_all@meta.data, is.character)], as.factor)

# Remove genes expressed in fewer than 5 cells
seurat_all <- DietSeurat(seurat_all, features = names(which(Matrix::rowSums(GetAssayData(seurat_all) > 0) >=5)))

# Store mitochondrial percentage in object meta data
seurat_all <- PercentageFeatureSet(seurat_all, pattern = "^MT-", col.name = "percent.mt")

# Remove data which do not pass filter threshold
seurat_all <- subset(seurat_all, subset = c(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15))


#####################################################################################################
#                           Integrate data from different 10x runs                                  #
####################################################################################################

# Split object by run and find integration points
seurat_split <- SplitObject(seurat_all, split.by = "run")

# Log normalize data and find variable features
seurat_split <- lapply(seurat_split, function(x) {
  NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Get array of all genes across all datasets in order to integrate using all features
all_features <- lapply(seurat_split, row.names) %>% Reduce(intersect, .)
# Find anchors used for integration
intergration_data <- Run.STACAS(seurat_split, dims = 1:30, anchor.features = 2000)
# Integrate data
intergration_data <- IntegrateData(anchorset = intergration_data, dims = 1:30, features.to.integrate = all_features)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(intergration_data) <- "integrated"

intergration_data <- ScaleData(intergration_data, features = rownames(intergration_data), vars.to.regress = "percent.mt")

# Save RDS after integration
saveRDS(intergration_data, paste0(rds_path, "intergration_data.RDS"))
