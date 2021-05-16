



# Sets design
Hib_design <- data.frame(row.names = colnames(Hib_count_data), condition = str_sub(colnames(Hib_count_data), start = 1L, end = 3L))

group <- as.factor(Hib_design$condition)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- colnames(Hib_count_data)

count_data <- Hib_count_data

# count_data_reads
col.cell <- brewer.pal(6, "Set3")[as.factor(Hib_design$condition)]

sum_data <- data.frame(Counts = c(mean(colSums(count_data)), (apply(count_data, 2, sum))))
sum_data <- sum_data %>% mutate(Samples = c("Average", colnames(count_data)))
sum_data <- sum_data %>% mutate(Colors = c("yellow", col.cell))
sum_data <- sum_data %>% mutate(Conditions = c("Average", Hib_design$condition %>% as.character()))

read_count_plot <- sum_data %>%
  ggplot(aes(x = Samples, y = Counts, fill = Conditions)) +
  geom_bar(colour = "black", stat = "identity") +
  geom_hline(mapping = NULL, data = NULL, yintercept = 2e+06, na.rm = FALSE, show.legend = NA, colour = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = "Read counts") +
  theme(legend.position="right") +
  scale_fill_manual(values=sum_data$Colors %>% unique())

read_count_plot

# Filter
d <- count_data %>% DGEList(lib.size = colSums(count_data), group = group)
keep <- rowSums(d %>% cpm() > 1) >= 1
filtered <- d[keep, ]
keep %>% summary()
filtered %>% dim()

filtered_count_data <- filtered$counts

# Least significantly DE genes based on a first-pass DE analysis performed prior to RUVg normalization.
fit <- filtered_count_data %>% DGEList(lib.size = colSums(filtered_count_data), group = group) %>% calcNormFactors(method = "TMM") %>% estimateDisp(design = design, tagwise = TRUE, robust = TRUE) %>% glmFit(design)
lrt <- fit %>% glmLRT(coef = 2:6)
top <- topTags(lrt, n = nrow(d))$table
empirical <- rownames(filtered_count_data)[which(!(rownames(filtered_count_data) %in% rownames(top)[1:500]))]

# Here, we consider all but the top 500 genes as ranked by edgeR p-values
set2 <- fit$counts %>% as.matrix() %>% RUVg(empirical, k = 1)

# PlotRLE
library(cowplot)
p1 <- ~ {
  filtered_count_data %>%
    as.matrix() %>%
    plotRLE(outline = FALSE, ylim = c(-2.5, 2.5), col = col.cell, main = "", cex.main = 2, las = 2, cex.axis = 0.8, style = "full", outlier.alpha = 0.1, outlier.shape = 3, outlier.size = 0, legend = TRUE)
  legend("topright", inset = c(0.0, 0), legend = group %>% levels(), col = col.cell %>% unique(), ncol = 2, cex = 1, box.lwd = 0, border = "black", fill = col.cell %>% unique())
}

p2 <- ~ {
  set2$normalizedCounts %>% 
    plotRLE(outline = FALSE, ylim = c(-2.5, 2.5), col = col.cell, main = "", cex.main = 2, las = 2, cex.axis = 0.8, style = "full", outlier.alpha = 0.1, outlier.shape = 3, outlier.size = 0)
  legend("topright", inset = c(0.0, 0), legend = group %>% levels(), col = col.cell %>% unique(), ncol = 2, cex = 1, box.lwd = 0, border = "black", fill = col.cell %>% unique())
}

plot_RLE <- plot_grid(p1, p2, labels = c("A", "B"), label_size = 12) + draw_label(label = "Relative Log Expression (RLE)", x = 0.04, y = 0.5, angle = 90) + draw_label(label = "Samples", x = 0.5, y = 0.05)
plot_RLE


# DE analysis
THRESHOLD <- 0.05
LFC <- 1
nSets <- 15

comparisons <- c("S01_x_S02", "S01_x_S03", "S01_x_S04", "S01_x_S05", "S01_x_S06", "S02_x_S03", "S02_x_S04", "S02_x_S05", "S02_x_S06", "S03_x_S04", "S03_x_S05", "S03_x_S06", "S04_x_S05", "S04_x_S06", "S05_x_S06")

contr.matrix <- limma::makeContrasts(
  con1 = S02 - S01,
  con2 = S03 - S01,
  con3 = S04 - S01,
  con4 = S05 - S01,
  con5 = S06 - S01,
  con6 = S03 - S02,
  con7 = S04 - S02,
  con8 = S05 - S02,
  con9 = S06 - S02,
  con10 = S04 - S03,
  con11 = S05 - S03,
  con12 = S06 - S03,
  con13 = S05 - S04,
  con14 = S06 - S04,
  con15 = S06 - S05,
  levels = c("S01", "S02", "S03", "S04", "S05", "S06")
)

# edgeR pipeline
fit <- set2$normalizedCounts %>% DGEList(lib.size = colSums(set2$normalizedCounts), group = group) %>% calcNormFactors("TMM") %>% estimateDisp(design, robust = TRUE) %>% glmFit(design)
rownames(fit$counts) <- set2$normalizedCounts %>% rownames()
lrt1 <- purrr::map(1:nSets, function(i) glmLRT(fit, contrast = contr.matrix[, i]))
edgeR_results <- purrr::map(seq_along(lrt1), function(i) topTags(lrt1[[i]], n = nrow(lrt1[[i]]), sort.by = "p.value"))
edgeR_results <- purrr::map(1:nSets, function(i) edgeR_results[[i]]$table) %>% purrr::map(rownames_to_column) %>% purrr::map(mutate, pi_value = abs(logFC) * (-1) * log10(FDR)) %>% purrr::map(dplyr::arrange, FDR) %>% purrr::set_names(c(comparisons))
edgeR_results <- edgeR_results %>%
  purrr::map(mutate, Product = if_else(is.na(match(rowname, as.character(RefSeq_gff.genes$Name))),
                                       as.character(RefSeq_gff.genes$product[match(rowname, as.character(RefSeq_gff.genes$OldLocusTag))]),
                                       as.character(RefSeq_gff.genes$product[match(rowname, as.character(RefSeq_gff.genes$Name))])
  ))
edgeR_top_DE1 <- edgeR_results %>% purrr::map(dplyr::filter, FDR < THRESHOLD)
edgeR_top_DE2 <- edgeR_top_DE1 %>% purrr::map(dplyr::filter, abs(logFC) >= LFC)
edgeR_DE_names <- edgeR_top_DE2 %>% purrr::map(dplyr::select, rowname) %>% purrr::map(pull) %>% purrr::set_names(c(comparisons))
edgeR_n_DEGs <- edgeR_DE_names %>% purrr::map(length)