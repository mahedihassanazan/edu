res <- read.csv("Retina/Retina_DEGs_Result.csv", row.names=1)
res <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]
# The rownames look like "ENSMUSG...|GeneName"
gene_names <- sapply(strsplit(rownames(res), "|", fixed=TRUE), function(x) x[2])
out <- data.frame(Gene_Name = gene_names, Log_Fold_Change = res$log2FoldChange, P_Value = res$pvalue)
write.csv(out, "Simplified_Retina_DEGs_No_NA.csv", row.names=FALSE)
cat("Done\n")
