library(ExperimentHub)
eh = ExperimentHub()
query(eh , "GSE62944")
tcga_data <- eh[["EH1"]]
unique(phenoData(tcga_data)$CancerType)

lgg_data <- tcga_data[, which(phenoData(tcga_data)$CancerType=="LGG")]
phe <- lgg_data@phenoData
colnames(phe)

mut_idx <- which(phenoData(lgg_data)$idh1_mutation_found=="YES")
mut_data <- exprs(lgg_data)[, mut_idx]

wt_idx <- which(phenoData(lgg_data)$idh1_mutation_found=="NO")
wt_data <- exprs(lgg_data)[, wt_idx]

countData <- cbind(mut_data, wt_data)

samples= c(colnames(mut_data), colnames(wt_data))
group =c(rep("mut",length(mut_idx)), rep("wt", length(wt_idx)))
coldata <- cbind(samples, group)
colnames(coldata) <- c("sampleName", "Group")
coldata[,"Group"] <- factor(coldata[,"Group"], c("wt","mut"))

library(DESeq2)
ddsMat <- DESeqDataSetFromMatrix(countData = countData,
                                 colData = DataFrame(coldata),
                                 design = ~ Group)

dds <- ddsMat
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

res <- results(dds) 
summary(res)
