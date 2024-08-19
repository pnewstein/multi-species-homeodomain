library(biomaRt)

drosophila <- "dmelanogaster_gene_ensembl"
mouse <- "mmusculus_gene_ensembl"

get_hdtf_status <- function(species) {
  ensembl <- useMart("ensembl")
  ensembl <- useEnsembl(biomart = "ensembl", dataset=species, mirror="useast")
  all_at <- listAttributes(mart = ensembl)
  at <- c("ensembl_gene_id", "external_gene_name", "smart")
  att_db <- getBM(attributes = at, mart = ensembl,)
  unique_genes = unique(att_db$ensembl_gene_id)
  out_df <- data.frame(
    ensembl_gene_id = rep(NA, length(unique_genes)),
    external_gene_name = rep(NA, length(unique_genes)),
    is_hdtf = rep(NA, length(unique_genes))
  )
  for (i in seq_along(unique_genes)) {
    this_gene <- att_db[att_db$ensembl_gene_id == unique_genes[i], ]
    out_df$ensembl_gene_id[i] <- this_gene$ensembl_gene_id[1]
    out_df$external_gene_name[i] <- this_gene$external_gene_name[1]
    out_df$is_hdtf[i] <- "SM00389" %in% this_gene$smart
  }
  return(out_df)
}

hdtfs <- get_hdtf_status(mouse)
write.csv(hdtfs, "data/mouse_hdtf.csv")
hdtfs <- get_hdtf_status(drosophila)
write.csv(hdtfs, "data/drosophila_hdtf.csv")
