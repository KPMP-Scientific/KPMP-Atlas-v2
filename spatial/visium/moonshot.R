library(ggplot2)

mstable <- rio::import_list('../data/new_samples/moonshot.xlsx')

mstable
df <- reshape2::melt(mstable[[1]])
df <- df[!is.na(df$value),]

ggplot(df,aes(x=variable,y=value,color=variable)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=.4,size=1) +
    theme_bw() +
    coord_flip()
ggsave('../figures/new_samples_abridged/moonshot_ckd.pdf',width=4,height=1.5)
write.csv(df,'R_final/data_source/Fig_ED5c_moonshot_CKD.csv',row.names = F)

df <- reshape2::melt(mstable[[2]])
df <- df[!is.na(df$value),]

ggplot(df,aes(x=variable,y=value,color=variable)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=.4,size=1) +
    theme_bw() +
    coord_flip()
ggsave('../figures/new_samples_abridged/moonshot_aki.pdf',width=4,height=1.5)
write.csv(df,'R_final/data_source/Fig_ED5c_moonshot_AKI.csv',row.names = F)
