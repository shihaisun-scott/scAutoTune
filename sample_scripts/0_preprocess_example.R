# preprocess data and create Seurat object
# this is just an example script for Seurat object construction and should not be treated as universally optimal
# preprocessing choices must be tailored to the specific dataset
# Note: this script does not produce any output, and is only treated as an example


library(Seurat)
library(Matrix)

exon_csv_file <- "enter_csv_exon_file_name.csv"
intron_csv_file <- "enter_csv_intron_file_name.csv"
metadata_csv_file <- "enter_metadata_file_name.csv"

exon <- read.csv(exon_csv_file, header = TRUE, row.names = 1)
intron <- read.csv(intron_csv_file, header = TRUE, row.names = 1)

data <- exon + intron
data <- sweep(data,2,colSums(data),"/")*10^6
data <- Matrix(as.matrix(data))
keepgen <- setdiff(rownames(data), grep("^LOC|^MT-", rownames(data), val = TRUE))
data <- data[keepgen,]
metadata <- read.csv(metadata_csv_file)
idx <- metadata$class_label != 'Low Quality'
keepcells1 = colnames(data[,idx])
metadata <- metadata[keepcells1,]
data <- data[,keepcells1]
object <- CreateSeuratObject(counts = data, meta.data = metadata)
saveRDS(object, file = "my_seurat_object.rds")
