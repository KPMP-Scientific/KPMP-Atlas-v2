library(Seurat)
library(ggplot2)
library(foreach)
library(doParallel)
library(abind)

spatial <- readRDS('../data/new_samples/spatial_niches_new_group_FTU.RDS')

target_niches <- paste0('PT_',c(3,34,22,10,11,4,30,5,38,27,28,48,14,41))


names(spatial@images)
head(spatial)

neighbors <- as.data.frame(matrix(NA,nrow = nrow(spatial@meta.data),ncol=6))
rownames(neighbors) <- rownames(spatial@meta.data)
for (smp in unique(spatial$orig.ident)){
  #smp <- unique(spatial$orig.ident)[1]
  suffix <- unique(stringr::str_split(rownames(spatial@meta.data[spatial$orig.ident == smp,]),'_',simplify = T)[,2])
  coords <- read.csv(paste0('../positions/',smp,'_positions.csv'))
  rownames(coords) <- paste0(coords$barcode,'_',suffix)
  for (sp in 1:nrow(coords)){
    spot <- rownames(coords)[sp]
    if(spot %in% rownames(spatial@meta.data)){
      i <- coords[sp,"array_row"]
      j <- coords[sp,"array_col"]
      neigh <- rownames(coords[coords$array_row==i-1 & coords$array_col == j-1,])
      neigh <- c(neigh,rownames(coords[coords$array_row==i-1 & coords$array_col == j+1,]))
      neigh <- c(neigh,rownames(coords[coords$array_row==i & coords$array_col == j-2,]))
      neigh <- c(neigh,rownames(coords[coords$array_row==i & coords$array_col == j+2,]))
      neigh <- c(neigh,rownames(coords[coords$array_row==i+1 & coords$array_col == j-1,]))
      neigh <- c(neigh,rownames(coords[coords$array_row==i+1 & coords$array_col == j+1,]))
      
      neighbors[spot,] <- neigh[1:6]
    }
    
  }
  
} 

# saveRDS(neighbors,'../data/new_samples/neighbors.RDS')
# neighbors <- readRDS('../data/new_samples/neighbors.RDS')

neighbors$new_group <- spatial@meta.data[rownames(neighbors),'new_group']
neighbors$niches_new_group <- spatial@meta.data[rownames(neighbors),'niches_new_group']
neighbors$full_niche <- paste0(neighbors$new_group,'_',neighbors$niches_new_group)

coloc_matrix <- matrix(0,nrow = length(target_niches),ncol = length(target_niches),
                       dimnames = list(target_niches,target_niches))


head(neighbors)
niche_neig <- neighbors
for (i in 1:nrow(neighbors)){
  #i <- 1  
  local <- neighbors[as.character(neighbors[i,1:6]),]
  niche_neig[i,1:6] <- local[1:6,"full_niche"]  
}
niche_neig <- niche_neig[apply(niche_neig,1,function(x) any(grepl('PT_',x))),]

coloc_matrix <-  coloc_matrix + as.matrix(table(niche_neig[,c("full_niche","V1")])[target_niches,target_niches])
coloc_matrix <-  coloc_matrix + as.matrix(table(niche_neig[,c("full_niche","V2")])[target_niches,target_niches])
coloc_matrix <-  coloc_matrix + as.matrix(table(niche_neig[,c("full_niche","V3")])[target_niches,target_niches])
coloc_matrix <-  coloc_matrix + as.matrix(table(niche_neig[,c("full_niche","V4")])[target_niches,target_niches])
coloc_matrix <-  coloc_matrix + as.matrix(table(niche_neig[,c("full_niche","V5")])[target_niches,target_niches])
coloc_matrix <-  coloc_matrix + as.matrix(table(niche_neig[,c("full_niche","V6")])[target_niches,target_niches])

coloc_matrix

diag(coloc_matrix) <- diag(coloc_matrix) / 2

#saveRDS(coloc_matrix,'../data/new_samples/coloc_matrix.RDS')
#coloc_matrix <- readRDS('../data/new_samples/coloc_matrix.RDS')

table(neighbors$full_niche)[target_niches]

df <- reshape2::melt(coloc_matrix)

df$Var1 <- factor(df$Var1,levels=paste0('PT_',c(3,34,22,10,11,4,30,5,38,27,28,48,14,41)))
df$Var2 <- factor(df$Var2,levels=paste0('PT_',c(3,34,22,10,11,4,30,5,38,27,28,48,14,41)))

df <- df[as.numeric(df$Var2) <= as.numeric(df$Var1),]
df$Var2 <- factor(df$Var2,levels=paste0('PT_',rev(c(3,34,22,10,11,4,30,5,38,27,28,48,14,41))))


ggplot(df,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()


df <- reshape2::melt(coloc_matrix / as.numeric(table(neighbors$full_niche)[target_niches]))
df$Var1 <- factor(df$Var1,levels=paste0('PT_',c(3,34,22,10,11,4,30,5,38,27,28,48,14,41)))
df$Var2 <- factor(df$Var2,levels=paste0('PT_',rev(c(3,34,22,10,11,4,30,5,38,27,28,48,14,41))))
ggplot(df,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  coord_fixed()+
  theme_minimal()

####################################################
### Permutation ran on cluster, code added below ###
####################################################

mean_sd_coloc <- readRDS('../data/new_samples/mean_sd_coloc.RDS')

zscore <- (coloc_matrix - mean_sd_coloc[['mean']] ) / mean_sd_coloc[['sd']]
pval <- 2*pnorm(abs(zscore),lower.tail = F)

df <- reshape2::melt(zscore)
colnames(df) <- c('Center','Neighbor','zscore')
df <- cbind(df,reshape2::melt(pval)[,3])
colnames(df) <- c('Center','Neighbor','zscore','pvalue')

df$Center <- factor(df$Center,levels=paste0('PT_',c(3,34,22,10,11,4,30,5,38,27,28,48,14,41)))
df$Neighbor <- factor(df$Neighbor,levels=paste0('PT_',c(3,34,22,10,11,4,30,5,38,27,28,48,14,41)))

df <- df[as.numeric(df$Neighbor) >= as.numeric(df$Center),]
df$Neighbor <- factor(df$Neighbor,levels=paste0('PT_',rev(c(3,34,22,10,11,4,30,5,38,27,28,48,14,41))))

df[df$pvalue==0,'pvalue'] <- min(df[df$pvalue>0,'pvalue'])

ggplot(df,aes(x=Center,y=Neighbor,color=zscore,size=-log10(pvalue)))+
  geom_point()+
  coord_fixed()+
  theme_classic()
ggsave('../figures/new_samples_abridged/dotplot_colocalization_neighbors.pdf',width = 5,height = 4.5)

df <- df[df$Center != df$Neighbor,]
ggplot(df,aes(x=Center,y=Neighbor,color=zscore,size=-log10(pvalue)))+
  geom_point()+
  coord_fixed()+
  theme_classic()
ggsave('../figures/new_samples_abridged/dotplot_colocalization_neighbors_nodiag.pdf',width = 5,height = 4.5)




####################################################
###               Permutation code               ###
####################################################

library(Seurat)
library(ggplot2)

args<-commandArgs(TRUE)
run_num <- as.numeric(args[1])
print(run_num)

spatial <- readRDS('spatial_niches_new_group_FTU.RDS')

target_niches <- paste0('PT_',c(3,34,30,5,22,10,11,4,38,27,28,48,14,41))

neighbors <- readRDS('neighbors.RDS')

neighbors$new_group <- spatial@meta.data[rownames(neighbors),'new_group']
neighbors$niches_new_group <- spatial@meta.data[rownames(neighbors),'niches_new_group']
neighbors$full_niche <- paste0(neighbors$new_group,'_',neighbors$niches_new_group)


randomize_colocalization <- function(association_matrix,target_niches){
  niche_association <- neighbors
  niche_association$full_niche <- sample(neighbors$full_niche,size = nrow(neighbors))
  
  for (i in 1:nrow(association_matrix)){
    #i <- 1  
    local <- association_matrix[as.character(association_matrix[i,1:6]),]
    niche_association[i,1:6] <- local[1:6,"full_niche"]  
  }
  niche_association <- niche_association[apply(niche_association,1,function(x) any(grepl('PT_',x))),]
  
  niche_colocalization_matrix <- matrix(0,nrow = length(target_niches),ncol = length(target_niches),
                                        dimnames = list(target_niches,target_niches))
  niche_colocalization_matrix <-  niche_colocalization_matrix + as.matrix(table(niche_association[,c("full_niche","V1")])[target_niches,target_niches])
  niche_colocalization_matrix <-  niche_colocalization_matrix + as.matrix(table(niche_association[,c("full_niche","V2")])[target_niches,target_niches])
  niche_colocalization_matrix <-  niche_colocalization_matrix + as.matrix(table(niche_association[,c("full_niche","V3")])[target_niches,target_niches])
  niche_colocalization_matrix <-  niche_colocalization_matrix + as.matrix(table(niche_association[,c("full_niche","V4")])[target_niches,target_niches])
  niche_colocalization_matrix <-  niche_colocalization_matrix + as.matrix(table(niche_association[,c("full_niche","V5")])[target_niches,target_niches])
  niche_colocalization_matrix <-  niche_colocalization_matrix + as.matrix(table(niche_association[,c("full_niche","V6")])[target_niches,target_niches])
  
  diag(niche_colocalization_matrix) <- diag(niche_colocalization_matrix) / 2
  return(niche_colocalization_matrix)  
}


random_coloc <- matrix(0,nrow=length(target_niches),ncol=length(target_niches),
                       dimnames = list(target_niches,target_niches))

set.seed(1000+run_num)

random_coloc <- randomize_colocalization(neighbors,target_niches)
saveRDS(random_coloc,paste0('random_coloc_run_',run_num,'.RDS'))
