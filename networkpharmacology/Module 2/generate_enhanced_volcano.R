options(repos = c(CRAN = "https://cloud.r-project.org"))
library(EnhancedVolcano)
library(ggplot2)

cat("Reading data...\n")
res <- read.csv("Retina/Retina_DEGs_Result.csv", row.names=1)

# Handle NA values
res$padj[is.na(res$padj)] <- 1
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0

# Extract just the gene symbol into a new column to avoid duplicate rowname errors
res$Symbol <- sapply(strsplit(rownames(res), "\\|"), function(x) x[2])

# Create custom colors: Left = Blue (Downregulated), Right = Red (Upregulated)
keyvals <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
      ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red2',
        'grey30'))
keyvals[is.na(keyvals)] <- 'grey30'
names(keyvals)[keyvals == 'red2'] <- 'Upregulated'
names(keyvals)[keyvals == 'grey30'] <- 'Not Significant'
names(keyvals)[keyvals == 'royalblue'] <- 'Downregulated'

cat("Generating EnhancedVolcano Plot...\n")
png("Enhanced_Volcano_Plot.png", width=1200, height=1000, res=150)
p <- EnhancedVolcano(res,
    lab = res$Symbol,
    x = 'log2FoldChange',
    y = 'padj',
    title = 'Glaucoma vs Normal',
    subtitle = 'EnhancedVolcano',
    pCutoff = 0.05,
    FCcutoff = 1,
    colCustom = keyvals,
    pointSize = 1.5,
    labSize = 4.0,
    legendPosition = 'top',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5) + 
    coord_cartesian(xlim = c(-12, 12))

print(p)
dev.off()

cat("Done! Plot saved as Enhanced_Volcano_Plot.png\n")
