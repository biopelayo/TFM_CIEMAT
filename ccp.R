prad <- as.matrix(tcga_RSEM_gene_fpkm_PRAD_mean_of_repeated_samples_lncRNAs_primary_332_de_cell_combat_batchid_338genes)




cat("\nI. Running Consensus Clustering...\n\n")
title<-c("PRAD_TFM_CCP")
k<-6 ##CHEK TO EXECUTE, maximum # of clusters
res <- ConsensusClusterPlus(prad, ##CHECK TO EXECUTE
                            maxK=k,
                            reps=1000,
                            pItem=0.8,
                            pFeature=1,
                            title=title,
                            innerLinkage="average",
                            finalLinkage="average",
                            clusterAlg="pam",
                            distance="manhattan",
                            writeTable = TRUE,
                            plot="pdf",
                            verbose=TRUE)


## datos simulados CCP, icl, PAC

cat("\nI. Running Consensus Clustering...\n\n")
title<-c("PRAD_simulated_TFM_CCP")
k<-6 ##CHEK TO EXECUTE, maximum # of clusters
ressimulated <- ConsensusClusterPlus(datasimulated, ##CHECK TO EXECUTE
                            maxK=k,
                            reps=1000,
                            pItem=0.8,
                            pFeature=1,
                            title=title,
                            innerLinkage="average",
                            finalLinkage="average",
                            clusterAlg="pam",
                            distance="manhattan",
                            writeTable = TRUE,
                            plot="pdf",
                            verbose=TRUE)





###################################################
### code chunk number 6: ConsensusClusterPlus.Rnw:84-96
###################################################
#consensusMatrix - the consensus matrix.  
#For .example, the top five rows and columns of results for k=2:
results[[2]][["consensusMatrix"]][1:5,1:5]

#consensusTree - hclust object 
results[[2]][["consensusTree"]]

#consensusClass - the sample classifications
results[[2]][["consensusClass"]][1:5]

#ml - consensus matrix result
#clrs - colors for cluster  


###################################################
### code chunk number 7: ConsensusClusterPlus.Rnw:104-105
###################################################
icl = calcICL(ressimulated,title=title,plot="pdf")


###################################################
### code chunk number 8: ConsensusClusterPlus.Rnw:109-110
###################################################
icl[["clusterConsensus"]]


###################################################
### code chunk number 9: ConsensusClusterPlus.Rnw:113-114
###################################################
icl[["itemConsensus"]][1:5,]






## PAC code
  
kvec <- 2:max(k)
x1 <- 0.1; x <- 0.9
PAC <- rep(NA, length(kvec))
names(PAC) <- paste("K=", kvec, sep = " ")
for(i in kvec){M = ressimulated[[i]]$consensusMatrix
Fn = ecdf(M[lower.tri(M)])
PAC[i-1] = Fn(x) - Fn(x1)
}#end for i
# The optimal K
optk = kvec[which.min(PAC)]
View(PAC)




