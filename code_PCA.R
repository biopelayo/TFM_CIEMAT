###########################################################################################
#          Msc. PELAYO G. DE LENA RODR?GUEZ
#            code for FactoMineR PCA
###########################################################################################


# start with a clean slate
rm(list=ls(all=TRUE)) 



#
# load libraries
library(FactoMineR)
library(ggplot2)
library(scales)
library(grid)
library(plyr)
library(gridExtra)
#



# load data express and data labels for cluster


express <- read.table("/home/bioinfo/Escritorio/PROYECTO_TFM_PRAD/tcga_RSEM_gene_fpkm_PRAD_mean_of_repeated_samples_lncRNAs_primary_332_de_cell_combat_batchid_338genes.txt",dec = ",",sep = "\t", stringsAsFactors = TRUE, header = TRUE )
#data1<-read.table("C:/Users/Julia/Desktop/classlabels.txt")
View(express)

# first column of data1 is rownames of express
lnclass$lncRNA_cluster_pelayo_combat_333_k4 ->express$lncRNA
  

View(express)

# columna cluster como factor con 4 niveles

prad$ <-as.factor(prad$cluster)
View(express)

# if all the columns needs to changed to numeric, 
# use lapply to loop over the columns and convert to numeric by 
# first converting it to character class as the columns were factor.

express[] <- lapply(express, function(x) as.numeric(as.character(x)))

# Both the columns in the OP's post are factor because of the string "n/a". 
# This could be easily avoided while reading the file using na.strings = "n/a" 
# in the read.table/read.csv or if we are using data.frame, we can have character columns 
# with stringsAsFactors=FALSE (the default is stringsAsFactors=TRUE)


###########################################################################################
#    P C A----ANALYSIS                                                                                     #
###########################################################################################

#

#
# compute PCA
res.pca <- PCA(PRAD_ml, quali.sup=1, graph = TRUE) # quali.sup es clusterlabel
#
# extract some parts for plotting
PC1 <- res.pca$ind$coord[,1]
PC2 <- res.pca$ind$coord[,2]
labs <- rownames(res.pca$ind$coord)
PCs <- data.frame(cbind(PC1,PC2))
rownames(PCs) <- labs
#
# Just showing the individual samples...
ggplot(PCs, aes(PC1,PC2, label=rownames(PCs))) + 
  geom_text() 
#
#
# Now get supplementary categorical variables
cPC1 <- res.pca$quali.sup$coor[,1]
cPC2 <- res.pca$quali.sup$coor[,2]
clabs <- rownames(res.pca$quali.sup$coor)
cPCs <- data.frame(cbind(cPC1,cPC2))
rownames(cPCs) <- clabs
colnames(cPCs) <- colnames(PCs)
#


# Put samples and categorical variables (ie. grouping
# of samples) all together
p <- ggplot() + theme_bw(base_size = 20) 


# add on data 
p <- p + geom_text(data=PCs, aes(x=PC1,y=PC2,label=rownames(PCs)), size=4) 
p <- p + geom_text(data=cPCs, aes(x=cPC1,y=cPC2,label=rownames(cPCs)),size=15)
p # show plot with both layers
#
# clear the plot
dev.off()
#

# Now do extract variables
#
vPC1 <- res.pca$var$coord[,1]
vPC2 <- res.pca$var$coord[,2]
vlabs <- rownames(res.pca$var$coord)
vPCs <- data.frame(cbind(vPC1,vPC2))
rownames(vPCs) <- vlabs
colnames(vPCs) <- colnames(PCs)
#

# and plot them
#
pv <- ggplot()  + theme_bw(base_size = 20) 

# put a faint circle there, as is customary
angle <- seq(-pi, pi, length = 50) 
df <- data.frame(x = sin(angle), y = cos(angle)) 
pv <- pv + geom_path(aes(x, y), data = df, colour="grey70") 
#
# add on arrows and variable labels
pv <- pv + geom_text(data=vPCs, aes(x=vPC1,y=vPC2,label=vlabs), size=4) + xlab("PC1") + ylab("PC2")
pv <- pv + geom_segment(data=vPCs, aes(x = 0, y = 0, xend = vPC1*0.9, yend = vPC2*0.9), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30")
pv # show plot 
#
# clear the plot
dev.off()
#
# Now put them side by side
#
library(gridExtra)
grid.arrange(p,pv,nrow=1)
# 
# Now they can be saved or exported...
#
# tidy up by deleting the plots
#
dev.off()

#Hierarchical Clustering on Principle Components (HCPC)
hcpc<-HCPC(res.pca, nb.clust=4, consol=FALSE, iter.max=100, min=3,  max=16, metric="manhattan", method="complete", order=TRUE, graph.scale="inertia", nb.par=5, graph=TRUE, proba=0.05,  cluster.CA="rows")
aaa<-summary(res.pcaind)
View(aaa)

## Copy results to files
write.table(express, "c:/Users/Julia/Desktop/PRED/analysis_PCA/dataset.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(vPCs, "c:/Users/Julia/Desktop/PRED/analysis_PCA/genesPCAr.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(PCs, "c:/Users/Julia/Desktop/PRED/analysis_PCA/samplesPCA.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(data1, "c:/Users/Julia/Desktop/PRED/analysis_PCA/dataset_labels.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.table(cPCs, "c:/Users/Julia/Desktop/PRED/analysis_PCA/cluster_PCA.txt", sep="\t", row.names = TRUE, col.names = TRUE)
write.csv(res.pcaind, "c:/Users/Julia/Desktop/PRED/analysis_PCA/summary_PCA.txt", sep="\t", row.names = TRUE, col.names = TRUE)




# To see whether the categories of the supplementary variable are significantly different from each other, we can draw confidence ellipses around them.
#To do so, type:
concat = cbind.data.frame(PRAD_ml$cluster,res.pca$ind$coord)
ellipse.coord = coord.ellipse(concat)
plot.PCA(res.pca,habillage=1,cex=0.5)

