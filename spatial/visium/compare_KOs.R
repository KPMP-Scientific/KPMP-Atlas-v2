library(ggplot2)

oracle <- read.csv("../celloracle_atlas/RNA_aPT-S1-2_GRN_aPT-S1-2/SOX4/SOX4_deg_before_after_ko.csv")
head(oracle)
oracle[oracle$pvals_adj == 0, "pvals_adj"] <- 2.348379e-308

# ko_genes <- c('CXCL11', 'DACH1', 'TNFRSF18', 'CCL2', 'CDKN1A',
#            'VCAM1', 'RUNX2', 'MID1', 'KHDRBS1', 'GDF15', 
#            'PROM1')
ko_genes <- c('CXCL11', 'DACH1', 'TNFRSF18', 'CCL2', 'CDKN1A', 'HNF1A', 
          'THRB', 'HNF4G', 'VCAM1', 'CEBPD', 'JUNB','GDF15', 'PROM1')

# Ensure ko_genes are filtered in the oracle data
any(!ko_genes %in% oracle$names)
oracle_filtered <- oracle[oracle$names %in% ko_genes, ]

# Create the ggplot
ggplot(oracle_filtered, aes(x = names, y = subclass, color = logfoldchanges, size = -log10(pvals_adj))) +
  geom_point() +
  scale_color_gradient2(low='#67001f',mid='white',high='#053061', midpoint = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "KO Genes", y = "Subclass", color = "Log Fold Change", size = "-log10(p-value)")
ggsave("../figures/dotplot_oracle.pdf", width = 8, height = 6)

oracle_filtered$adjusted_size <- pmin(pmax(-log10(oracle_filtered$pvals_adj), 0), 30)
oracle_filtered$adjusted_size <- pmin(-log10(oracle_filtered$pvals_adj), 30)

ggplot(oracle_filtered, aes(x = names, y = subclass, color = logfoldchanges, size = -log10(pvals_adj))) +
  geom_point() +
  scale_color_gradient2(low='#67001f',mid='white',high='#053061', midpoint = 0,
                        limits=c(-0.1,0.1),oob=scales::squish) +
  scale_size_continuous(limits = c(0, 30)) +
  scale_x_discrete(limits = ko_genes) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "KO Genes", y = "Subclass", color = "Log Fold Change", size = "-log10(p-value)")
ggsave("../figures/dotplot_oracle_rescaled_abridged.pdf", width = 4.5, height = 1.7)

invitro  <- read.csv('../siRNA/condition_SOX4_vs_Negative.csv')
head(invitro)

any(!ko_genes %in% invitro$annot.symbol)
invitro_filtered <- invitro[invitro$annot.symbol %in% ko_genes, ]
invitro_filtered$siRNA <- 'siRNA'

ggplot(invitro_filtered, aes(x = annot.symbol, y = siRNA, color = log2FoldChange, size = -log10(padj))) +
  geom_point() +
  scale_color_gradient2(low='#67001f',mid='white',high='#053061', midpoint = 0) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "siRNA Genes", y = "Subclass", color = "Log Fold Change", size = "-log10(p-value)")
ggsave("../figures/dotplot_siRNA.pdf", width = 8, height = 1.5)

ggplot(invitro_filtered, aes(x = annot.symbol, y = siRNA, color = log2FoldChange, size = -log10(padj))) +
  geom_point() +
  scale_color_gradient2(low='#67001f',mid='white',high='#053061', midpoint = 0,
                        limits=c(-1,1),oob=scales::squish) +
  scale_x_discrete(limits = ko_genes) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "siRNA Genes", y = "Subclass", color = "Log Fold Change", size = "-log10(p-value)")
ggsave("../figures/dotplot_siRNA_rescaled_abridged.pdf", width = 4.2, height = 1.5)

invitro_filtered[order(invitro_filtered$log2FoldChange, decreasing = TRUE), 'annot.symbol']
invitro_filtered[order(invitro_filtered$log2FoldChange, decreasing = TRUE), 'log2FoldChange']


# Scale oracle$logfoldchanges to -1 to 1, preserving sign
maxFC <- max(abs(oracle$logfoldchanges), na.rm = TRUE)
oracle$scaled <- oracle$logfoldchanges / maxFC
max(oracle$scaled)
#scale invitro$log2FoldChange to -1 to 1, preserving sign
invitro_scaled <- invitro[!is.na(invitro$log2FoldChange), ]
maxFC <- max(abs(invitro_scaled$log2FoldChange), na.rm = TRUE)
invitro_scaled$scaled <- invitro_scaled$log2FoldChange / maxFC
max(invitro_scaled$scaled)

# keep only the genes that are in both oracle and invitro_scaled
oracle <- oracle[oracle$names %in% invitro_scaled$annot.symbol, ]
invitro_scaled <- invitro_scaled[invitro_scaled$annot.symbol %in% oracle$names, ]

#multiply oracle$scaled by invitro_scaled$scaled where they match and are not NA
oracle$agrement <- oracle$scaled * invitro_scaled$scaled[match(oracle$names, invitro_scaled$annot.symbol)]
write.csv(oracle, '../celloracle_atlas/RNA_aPT-S1-2_GRN_aPT-S1-2/SOX4/SOX4_deg_before_after_ko_matched.csv', row.names = FALSE)
