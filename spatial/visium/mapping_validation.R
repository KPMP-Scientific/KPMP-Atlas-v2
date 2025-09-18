library(Seurat)
library(ggplot2)

spatial <- readRDS('../data/new_samples/spatial_niches_new_group_FTU.RDS')

annot <- rio::import('../data/Supplementary Tables Version 1 05-2024.xlsx',
                    sheet= 'Table S12 - sn markers',range='B5:N133')
rownames(annot) <- annot$v2.subclass.l3
head(annot)

histo_180006 <- read.csv('../data/mapping_validation/A1_S1_18-0006_JCI.csv')
rownames(histo_180006) <- histo_180006$Barcode
head(histo_180006)
histo_180006$annot <- stringr::str_split(histo_180006$Graph.based,' ',simplify = TRUE)[,1]
histo_180006$annot <- ifelse(histo_180006$annot %in% c('DCT-CNT','CD'), 'CNT-CD', histo_180006$annot)

unique(spatial$orig.ident)
st <- subset(spatial, orig.ident == 'V19S25-016_XY01_18-0006')

data <- reshape2::melt(as.matrix(st@assays$predv2subclassl3@data))
head(data)
data <- data[data$Var1 != 'max',]
data$Var2 <- stringr::str_split(data$Var2,pattern = '_',simplify = TRUE)[,1]
data$annot.symbol <- histo_180006[match(data$Var2,rownames(histo_180006)),'annot']
head(data)

tempdf <- aggregate(value ~ Var1 + annot.symbol, data = data, FUN = sum)
head(tempdf)


library(dplyr)

tempdf <- tempdf %>%
  group_by(Var1) %>%
  mutate(scaled_value = value / sum(value)) %>%
  ungroup()

tempdf$Var1 <- factor(tempdf$Var1,levels = annot$v2.subclass.l3)
unique(tempdf$annot.symbol)
tempdf$annot.symbol <- factor(tempdf$annot.symbol,levels = unique(tempdf$annot.symbol)[c(4,1,2,6,5,3)])


ggplot(filter(tempdf, !is.na(scaled_value)), aes(x = Var1, y = annot.symbol, fill = scaled_value)) +
  geom_tile() +
  scale_fill_gradient(low='white',high = 'firebrick') +
  coord_equal() +
  theme(axis.text.x = element_text(angle = 90,vjust = .1,hjust=1))

ggsave('../data/mapping_validation/heatmap.pdf',width=16,height = 4)
write.csv(tempdf[,c('Var1','annot.symbol','scaled_value')],'R_final/data_source/SD4_a_validation_heatmap.csv',row.names = F)
