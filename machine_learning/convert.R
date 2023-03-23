library(SeuratDisk)

Convert("/workdir/cfrna/references/human/tabula-sapiens/figshare/TabulaSapiens.h5ad", "h5seurat")

counts <- Seurat::LoadH5Seurat(
  "/workdir/cfrna/references/human/tabula-sapiens/figshare/TabulaSapiens.h5seurat",
  use.names = F
)

save(counts, "ts.Rdata")
