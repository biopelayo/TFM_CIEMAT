tabla_total <- read.delim("C:/Users/Julia/Desktop/Pela/tabla_total.txt", row.names=1)
View(tabla_total)
dim(tabla_total)
df <- t(as.data.frame(tabla_total)
)
rownames(df)<- df$sampleID
View(df)
attach(df)
class(df)
df <- as.data.frame(df)
class(df)
attach(df)
table(c2, Subtype)
tablasubt <- table(c2, Subtype)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
summary(a$observed)
table(c2, iCluster)
tablasubt <- table(c2, iCluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
summary(a$observed)
plot(a)
plot(tablasubt)
table(c2, Age)
tablasubt <- table(c2, Age)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
table(c2, mRNA_cluster)
tablasubt <- table(c2, mRNA_cluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
table(c2, methylation_cluster)
tablasubt <- table(c2, methylation_cluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
summary(c2)
table(c2, SPOP_mut)
tablasubt <- table(c2, SPOP_mut)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
summary(a$observed)
plot(a)
plot(tablasubt)
table(c2, PTEN_CNA)
tablasubt <- table(c2, PTEN_CNA)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
table(c2, SPOPL_CNA)
tablasubt <- table(c2, SPOPL_CNA)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
table(c2, df$Race)
tablasubt <- table(c2, df$Race)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
summary(a$observed)
plot(a)
plot(tablasubt)
table(c2, df$Clinical_Gleason_category)
tablasubt <- table(c2, df$Clinical_Gleason_category)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
table(c2, df$Mutations)
tablasubt <- table(c2, df$Mutations)
chisq.test(tablasubt)
df$Mutations <- as.character(df$Mutations)
table(c2, df$Mutations)
df$Mutations <- as.numeric(df$Mutations)
table(c2, df$Mutations)
tablasubt <- table(c2, df$Mutations)
chisq.test(tablasubt)
plot(tablasubt)
table(c2, df$miRNA_cluster)
tablasubt <- table(c2, df$miRNA_cluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
table(c2, df$RPPA_cluster)
tablasubt <- table(c2, df$RPPA_cluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
summary(a$observed)
plot(a)
plot(tablasubt)
table(c2, df$SCNA_cluster)
tablasubt <- table(c2, df$SCNA_cluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
summary(a$observed)
plot(a)
plot(tablasubt)
table(c2, df$RB1_mut)
tablasubt <- table(c2, df$RB1_mut)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
table(c2, df$TP53_CNA)
tablasubt <- table(c2, df$TP53_CNA)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
table(c2, df$TP53_mut)
tablasubt <- table(c2, df$TP53_mut)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
table(c2, df$PSA_preop)
tablasubt <- table(c2, df$PSA_preop)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
table(c3, df$TP53_mut)
tablasubt <- table(c3, df$TP53_mut)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
table(c1, df$mRNA_cluster)
tablasubt <- table(c1, df$mRNA_cluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
a$residuals
dim(a)
summary(a$observed)
table(c1, df$methylation_cluster)
tablasubt <- table(c1, df$methylation_cluster)
chisq.test(tablasubt)
plot(tablasubt)
summary(tablasubt)
chisq.test(tablasubt, simulate.p.value = TRUE)
summary(tablasubt)
a<-chisq.test(tablasubt, simulate.p.value = TRUE)
summary(a)
a$p.value
a$statistic
a$parameter
a$method
a$data.name
a$expected
a$observed
