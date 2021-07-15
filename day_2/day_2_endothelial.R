library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(Glimma)
library(EnsDb.Hsapiens.v79)
library(Trendy)
library(gplots)

library("pheatmap")
library(EnhancedVolcano)
library(DESeq2)
library(e1071)

count_data = read.csv("C:/drug discovery hackathon/day_2/counts_endothelial.csv", row.names = 1)
clin_data = data.frame("Label" = c(0,1,1,1,1,0,0,0,0,0), row.names = colnames(count_data) )
clin_data$Label = as.factor(clin_data$Label)
levels(clin_data$Label) = c("normal", "diabetes")
dds = DESeqDataSetFromMatrix(countData = count_data, colData = clin_data, design = ~Label)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Label","normal","diabetes"))
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$pvalue < 0.05)

#Volcano plot
EnhancedVolcano(res,lab = rownames(res), title = "Diabetes versus control",x = 'log2FoldChange',y = 'pvalue',xlim = c(-5, 8))

#MA plot
plotMA(resLFC1)

#pheatmap
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Label")])
test_df = assay(ntd)[select,]
colnames(test_df) = rownames(df)
pheatmap(test_df, cluster_rows=FALSE, show_rownames=T,
         cluster_cols=T, annotation_col=df, main = "Diabetes versus Normal")


genes = resLFC1[which((resLFC1$pvalue < 0.05) == TRUE),]
gene_list = rownames(genes)
write.csv(data.frame(gene_list), "C:/drug discovery hackathon/files/sign_genes.csv")

hypoxia_df = read.csv("C:/drug discovery hackathon/files/GO_term_summary_20200915_170513.csv")
hypoxia_genes = hypoxia_df$Symbol
hypoxia_genes = toupper(hypoxia_genes)
hypoxia_diabetes_genes = intersect(gene_list,hypoxia_genes)


