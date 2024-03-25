library(Seurat)

dataset_dir <- '/export/home/s2992547/scratch/'
liver_cosmx_seurat_file <- file.path(dataset_dir,'LiverDataReleaseSeurat_newUMAP.RDS')
so <- readRDS(liver_cosmx_seurat_file )

# summary of samples
table(so$fov, so$orig.ident)

# save all metadata to disk
write.table(so@meta.data, file="LiverData_metadata.tsv", sep='\t')

# Take an assortment of fovs. Deliberately non-contigous.
keep_fovs <- c(40,60,80,100,120,140,160,
               340,360,380,400,420,440,460)

so <- so[,so$fov %in% keep_fovs]
dim(so)
saveRDS(so,file="LiverDataReleaseSeurat_subset.RDS")
