# TPM
TPM values calculation and display for RNA-seq experiments



```
library(BiocManager)

BiocManager::install("GenomicFeatures")

library(GenomicFeatures)

install.packages("pheatmap")

library(pheatmap)


calcTPM <- function(counts, lengths) {
  lengths_kb <- lengths / 1000
  rate <- counts / lengths_kb
  tpm <- rate / sum(rate) * 1e6
  return(tpm)
}

setwd("")

count_file <- read.table(file = "count_table_Arsenic_Grazie_4132023.txt",
                         row.names = 1,
                         sep = "\t", header = TRUE
                      )


gtf_file <- "Arabidopsis_thaliana.TAIR10.60.gff3.gz"

txdb <- makeTxDbFromGFF(gtf_file, format = "gff3")

exons_by_gene <- exonsBy(txdb, by = "gene")

gene_lengths <- vapply(exons_by_gene, function(gr) sum(width(reduce(gr))), numeric(1))



common_genes <- intersect(rownames(count_file), names(gene_lengths))

if (length(common_genes) == 0) {
  stop("No overlapping gene IDs between the count table and the GTF-derived gene lengths.")
}

counts_data <- count_file[common_genes, ]
gene_lengths <- gene_lengths[common_genes]


attr(count_file, "gene_lengths") <- gene_lengths

counts_data[] <- lapply(counts_data, function(x) as.numeric(as.character(x)))

tpm_matrix <- apply(counts_data, 2, calcTPM, lengths = gene_lengths)


tpm_matrix <- as.data.frame(tpm_matrix)

write.table(tpm_matrix, file = "TPM_values.txt", sep = "\t", quote = FALSE, col.names = NA)


genes_of_interest <- c("AT1G53430", "AT1G53440", "AT5G01950")



heatmap_data <- tpm_matrix[rownames(tpm_matrix) %in% genes_of_interest, ]

if (nrow(heatmap_data) == 0) {
  stop("None of the specified genes of interest were found in the TPM matrix")
}

log_tpm <- log2(heatmap_data + 1)


metadata <- read.csv("Sample_Metadata.csv", row.names = 1)

pheatmap(log_tpm,
         cluster_rows = TRUE, cluster_cols = FALSE, annotation_col =  metadata,
         main = "Heatmap of TPM Values for Selected Arabidopsis Genes")

selected_gene <- "AT1G63440"
if (selected_gene %in% rownames(tpm_matrix)) {
  gene_data <- data.frame(Sample = colnames(tpm_matrix), TPM = as.numeric(tpm_matrix[selected_gene, ]))
  
  # Merge with metadata to get group information
  gene_data <- merge(gene_data, metadata, by.x = "Sample", by.y = "row.names", all.x = TRUE)
  
ggplot(gene_data, aes(x = factor(group), y = TPM, fill = factor(group))) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_jitter(width = 0.1, alpha = 0.8, color = "black") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("Violin Plot of TPM Values for HMA5-", selected_gene), x = "Group", y = "TPM Value") +
    theme(legend.position = "none") + theme_bw()
} else {
  warning(paste("Selected gene", selected_gene, "not found in TPM matrix"))
}

ggplot(gene_data, aes(x = factor(group), y = TPM, fill = factor(group))) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.7, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = paste("Violin Plot of TPM Values for", selected_gene), x = "Group", y = "TPM Value") +
  theme(legend.position = "none") + theme_bw() +
  stat_compare_means(method = "kruskal.test")

summary_data <- gene_data %>%
  group_by(group) %>%
  summarise(mean_TPM = mean(TPM), 
            sd_TPM = sd(TPM), 
            n = n(),
            se_TPM = sd_TPM / sqrt(n))

ggplot(summary_data, aes(x = factor(group), y = mean_TPM, fill = factor(group))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_TPM - se_TPM, ymax = mean_TPM + se_TPM), width = 0.2) +
  labs(title = paste("Bar Plot of Mean TPM Values for", selected_gene), x = "Group", y = "Mean TPM Value") +
  theme_minimal() +
  ggpubr::stat_compare_means(method = "kruskal.test", label = "p.signif")

} else {
  warning(paste("Selected gene", selected_gene, "not found in TPM matrix"))
}


```


### Insights into HMA family

![image](https://github.com/user-attachments/assets/221fe00b-81b5-4fbc-bd61-3c54e4d78c63)

![image](https://github.com/user-attachments/assets/ef31d35e-117d-4506-99e3-2fe531b29990)



