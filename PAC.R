###########################################################################################
#          Msc. PELAYO G. DE LENA RODRìGUEZ
#            code for Proportion of Ambiguous Clustering (PAC)
###########################################################################################
  
kvec <- 2:max(k) # 
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
