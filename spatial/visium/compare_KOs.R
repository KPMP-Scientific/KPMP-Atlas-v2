library(ggplot2)

oracle <- read.csv("../celloracle_atlas/RNA_aPT-S1-2_GRN_aPT-S1-2/SOX4/SOX4_deg_before_after_ko.csv")
head(oracle)
oracle[oracle$pvals_adj == 0, "pvals_adj"] <- 2.348379e-308

ko_genes <- c('CXCL11', 'DACH1', 'TNFRSF18', 'CCL2', 'CDKN1A', 'HNF1A', 
          'THRB', 'HNF4G', 'VCAM1', 'CEBPD', 'JUNB','GDF15', 'PROM1')

any(!ko_genes %in% oracle$names)
oracle_filtered <- oracle[oracle$names %in% ko_genes, ]

oracle_filtered$adjusted_size <- pmin(pmax(-log10(oracle_filtered$pvals_adj), 0), 30)
oracle_filtered$adjusted_size <- pmin(-log10(oracle_filtered$pvals_adj), 30)

ocplot <- ggplot(oracle_filtered, aes(x = names, y = subclass, color = logfoldchanges, size = -log10(pvals_adj))) +
            geom_point() +
            scale_color_gradient2(low='#67001f',mid='white',high='#053061', midpoint = 0,
                                  limits=c(-0.1,0.1),oob=scales::squish) +
            scale_size_continuous(limits = c(0, 30)) +
            scale_x_discrete(limits = ko_genes) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "KO Genes", y = "Subclass", color = "Log Fold Change", size = "-log10(p-value)")
ocplot            
ggsave("../figures/dotplot_oracle_rescaled_abridged.pdf", width = 4.5, height = 1.7)
write.csv(ocplot$data,'R_final/data_source/Fig_5e_insilico_dot.csv',row.names = F)

invitro  <- read.csv('../siRNA/condition_SOX4_vs_Negative.csv')
head(invitro)

any(!ko_genes %in% invitro$annot.symbol)
invitro_filtered <- invitro[invitro$annot.symbol %in% ko_genes, ]
invitro_filtered$siRNA <- 'siRNA'

ivplot <- ggplot(invitro_filtered, aes(x = annot.symbol, y = siRNA, color = log2FoldChange, size = -log10(padj))) +
            geom_point() +
            scale_color_gradient2(low='#67001f',mid='white',high='#053061', midpoint = 0,
                                  limits=c(-1,1),oob=scales::squish) +
            scale_x_discrete(limits = ko_genes) +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = "siRNA Genes", y = "Subclass", color = "Log Fold Change", size = "-log10(p-value)")
ivplot            
ggsave("../figures/dotplot_siRNA_rescaled_abridged.pdf", width = 4.2, height = 1.5)
write.csv(ivplot$data,'R_final/data_source/Fig_5e_invitro_dot.csv',row.names = F)
