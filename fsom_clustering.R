library(readxl)
library(CATALYST)
library(SingleCellExperiment)
library(diffcyt)
library(HDCytoData)
library(DT)
library(ggplot2)
library(cowplot)
library(FlowSOM)
library(ConsensusClusterPlus)

#load the data and construct the sce file
md <- read_excel("data/meta_11vs8.xlsx")
panel <- read_excel("data/panel_umap.xlsx")
files <- list.files(path = "files_na", pattern = "\\.fcs$", full.names = TRUE)
#exclude cvd004b
# very important that all loaded fcs files are in the meta data and vice versa
files <- files[-c(4,20)]

fs <- read.flowSet(files, transformation = FALSE, truncate_max_range = FALSE)
fs = fs[,panel$fcs_colname]
md$condition <- factor(md$condition, levels = c("healthy", "patient"))

md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])
ids1 <- fsApply(fs, identifier)

md = subset(md,file_name %in% ids1)
sce <- prepData(fs, panel, md, features = panel$fcs_colname, cofactor = 5)


get_features = function (x, fs) 
{
  if (is.null(fs)) {
    fs <- rownames(x)
  }
  else {
    stopifnot(is.character(fs))
    foo <- tryCatch(error = function(e) e, match.arg(fs, 
                                                     c("type", "state", "none")))
    if (!inherits(foo, "error")) {
      fs <- foo
      stopifnot(!is.null(marker_classes(x)))
      fs <- rownames(x)[marker_classes(x) == fs]
      if (length(fs) == 0) 
        stop("No features matched the specified marker class.")
    }
    else stopifnot(fs %in% rownames(x))
  }
  return(fs)
}

features = "type"
features <- get_features(sce, features)
if (is.null(marker_classes(sce))) {
  rowData(sce)$marker_class <- factor(c("state", "type")[as.numeric(rownames(scw) %in% 
                                                                    features) + 1], levels = c("type", "state", "none"))
}
rowData(sce)$used_for_clustering <- rownames(sce) %in% features
es = t(assay(sce, "exprs"))


all_data <- data.frame(colData(sce), es)

only_healthy = all_data[all_data$condition == "healthy",]
only_patients= all_data[all_data$condition == "patient",]



# all data representation
fsom <- ReadInput(flowFrame(t(assay(sce, "exprs"))))

xdim = 10
ydim = 10
maxK = 20
seed = 1234

som <- BuildSOM(fsom, colsToUse = features, silent = TRUE, 
                xdim = xdim, ydim = ydim)

fSOM <- BuildMST(som)

k <- xdim * ydim
mcs <- seq_len(maxK)[-1]
mc <- suppressWarnings(suppressMessages(ConsensusClusterPlus(t(som$map$codes), 
                                                             maxK = maxK, reps = 100, distance = "euclidean", seed = seed, 
                                                             plot = NULL)))
# number of clusters selected
n_clust = 12
code_clustering1 <- mc[[n_clust]]$consensusClass

cell_clustering_som <- som$map$mapping[,1]
cell_clustering1 <- code_clustering1[cell_clustering_som]
PlotStars(fSOM,backgroundValues = as.factor(code_clustering1))



# healthy vs patients separetly

fsom_healthy = ReadInput(flowFrame(data.matrix(only_healthy)),scale = FALSE)
fsom_patients = ReadInput(flowFrame(data.matrix(only_patients)),scale = FALSE)

som_healthy <- BuildSOM(fsom_healthy, colsToUse = features, silent = TRUE,
                xdim = xdim, ydim = ydim)
som_patients<- BuildSOM(fsom_patients, colsToUse = features, silent = TRUE,
                        xdim = xdim, ydim = ydim)

fSOM_healthy <- BuildMST(som_healthy, tSNE=FALSE)


fSOM_patientss <- BuildMST(som_patients, tSNE=FALSE)



#### building merged objects



fsom_mhealthy = fSOM_healthy

fsom_mhealthy$MST$graph <- fSOM$MST$graph
fsom_mhealthy$MST$l <- fSOM$MST$l

fsom_mpatients = fSOM_patientss

fsom_mpatients$MST$graph <- fSOM$MST$graph
fsom_mpatients$MST$l <- fSOM$MST$l

PlotStars(fsom_mhealthy,backgroundValues = as.factor(code_clustering1))
PlotStars(fsom_mpatients,backgroundValues = as.factor(code_clustering1))

# umap
set.seed(1234)
sce <- runDR(sce, "UMAP", cells =min(n_cells(sce)), features = "type")
sce$cluster_id <- factor(cell_clustering1)
metadata(sce)$cluster_codes <- cell_clustering1 
metadata(sce)$SOM_codes <- som$map$codes


dr = "UMAP"
assay = "exprs"
color_by = "condition"
dims = c(1,2)
xy <- reducedDim(sce, dr)[, dims]
clusters = sce$cluster_id
es <- as.matrix(assay(sce, assay))
colnames(xy) <- c("x", "y")
df <- data.frame(colData(sce), xy, t(es),clusters)
df <- df[!(is.na(df$x) | is.na(df$y)), ]

outDir = "data/umap_na.csv"
write.csv(df,outDir, row.names = FALSE)

pdf(file="full_trees_na.pdf")

# star plots
a = PlotStars(fsom_mhealthy,backgroundValues = as.factor(code_clustering1))
title("Not activated, healthy")



b = PlotStars(fsom_mpatients,backgroundValues = as.factor(code_clustering1))
title("Not activated, patients")
dev.off()

ms = c('CD63', 'CD107a', 'PAC1', 'CD62P', 'CD154')
pdf(file="trees_healthy_na.pdf")
for (m in ms){
  PlotMarker(fsom_mhealthy,m,backgroundValues = as.factor(code_clustering1))
}
dev.off()

ms = c('CD63', 'CD107a', 'PAC1', 'CD62P', 'CD154')
pdf(file="trees_patients_na.pdf")
for (m in ms){
  PlotMarker(fsom_mpatients,m,backgroundValues = as.factor(code_clustering1))
}
dev.off()

