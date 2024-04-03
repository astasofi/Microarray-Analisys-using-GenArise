#Mestrado Bioinformatica e Biologia Computacional
#Métodos Estatísticos em Bioinformática 20/21
#Ana Sofia Almeida - 49292

# Packages 
library(genArise)
library(limma)

my_spot_1 <- read.spot("chip1.txt",id=1, cy3=3, cy5=2, bg.cy3=5, bg.cy5=4, header = TRUE, sep = "\t", is.ifc = FALSE)

c1 <- read.table("chip1.txt", header = T)

c2 <- read.table("chip2.txt", header = T)

c3 <- read.table("chip3.txt", header = T)


#a) Plot the intensities observed in the two channels, through the dispersion
#diagram and the MA-plot.

# Dispersion diagram
cys.plot(my_spot_1)
plot(c1$Ven, c1$Art, xlab = "Cy3", ylab = "Cy5")

# Histograms Raw vs Log
#Histograma Art
hist(c1$Art,  xlab = "Cy5")

#Histograma Ven
hist(c1$Ven,  xlab = "Cy3")

#Histograma log(Art)
hist(log(c1$Art), xlab = "log(Cy5)")

#Histograma log(Ven)
hist(log(c1$Ven), xlab = "log(Cy3)")

# MA-plot
ma.plot(my_spot_1)


# b) Subtract the background in both channels, and compare the results using
#the convenient graphical representation.

# Background correction
my_spot_bc <- bg.correct(my_spot_1)

# MA-plot representation
ma_bc <- ma.plot(my_spot_bc)


# c) Normalize the array and compare the results again.

# Normalization
my_spot_norm <- global.norm(mySpot = my_spot_bc)

# MA-plot after the normalization
ma_norm <- ma.plot(my_spot_norm)


# (d) Calculate the Z-score for each gene. Construct the appropriate graph and
#use it to establish a cut-off point that allows you to create lists of genes
#deferentially expressed in Art and Ven.

my_spot_z <- Zscore(my_spot_norm, type="ma") #using M-A values 
Zscore.plot(my_spot_z)

c1_zscores <- data.frame(Ids = my_spot_z@dataSets$Id, Zscore = my_spot_z@dataSets$Zscore)
nrow(c1_zscores)

#mean of z-score
mean(my_spot_z@dataSets$Zscore)
#standard deviation of z-score
sd(my_spot_z@dataSets$Zscore)

c1_zscores <- c1_zscores[c1_zscores$Zscore < -2 | c1_zscores$Zscore > 2, ]
nrow(c1_zscores)

#down regulated - the ones with z-score lower than -2
down_reg <- data.frame(c1_zscores[c1_zscores$Zscore < -2, ])
nrow(down_reg)

#up regulated - the ones with z-scores higher than 2
up_reg <- data.frame(c1_zscores[c1_zscores$Zscore > 2, ])
nrow(up_reg)



# (e) Normalize both arrays taking into account 
#the normalization process considered above.

my_spot_2 <- read.spot("chip2.txt",id = 1, cy3 = 3, cy5 = 2, bg.cy3 = 5, bg.cy5 = 4, header = TRUE, sep = "\t", is.ifc = FALSE)

my_spot_3 <- read.spot("chip3.txt",id = 1, cy3 = 3, cy5 = 2, bg.cy3 = 5, bg.cy5 = 4, header = TRUE, sep = "\t", is.ifc = FALSE)

# background correction for patient 2
my_spot_bc2 <- bg.correct(my_spot_2)
ma.plot(my_spot_bc2)

# intra-array normalization for patient 2
my_spot_norm2 <- global.norm(mySpot = my_spot_bc2)
ma.plot(my_spot_norm2) #MA-plot

# background correction for parient 3
my_spot_bc3 <- bg.correct(my_spot_3)
ma.plot(my_spot_bc3)

# intra-array normalization for patient 3
my_spot_norm3 <- global.norm(mySpot = my_spot_bc3)
ma.plot(my_spot_norm3 )


#(f) Proceed with the normalization between arrays 
#by doing a convenient transformation of the data 
#(suggestion: centerin)

#Logarithm of the ratio (Cy5/Cy3) for the 3 arrays
m1 <- log(my_spot_norm@spotData$Cy5, 2) - log(my_spot_norm@spotData$Cy3, 2)

m2 <- log(my_spot_norm2@spotData$Cy5, 2) - log(my_spot_norm2@spotData$Cy3, 2) 

m3 <- log(my_spot_norm3@spotData$Cy5, 2) - log(my_spot_norm3@spotData$Cy3, 2) 


#Inter-array normalization - centralization + scaling
m1_cent <- (m1-mean(m1))/sd(m1)

m2_cent <- (m2-mean(m2))/sd(m2)

m3_cent <- (m3-mean(m3))/sd(m3)

m_cent_all <- data.frame(Id=my_spot_1@spotData$Id, chip1 = m1_cent, chip2 = m2_cent, chip3 = m3_cent)

#Boxplots
par(mfrow=c(1,2))

#boxplot with the M values prior to the inter-array normalizaation
ma_norm <- data.frame(chip1 = m1, chip2 = m2, chip3 = m3)
boxplot(ma_norm, names = colnames(ma_norm), col = rainbow(5), main = "(a)")

#boxplot with the M values after the inter-array normalizaation
ma_cent <- data.frame(chip1 = m1_cent, chip2 = m2_cent, chip3 = m3_cent)
boxplot(ma_cent, names = colnames(ma_cent), col = rainbow(5), main = "(b)")


#(g) Based on the Z-scores identify the genes 
#with differential expression in the three arrays.
#Comment this procedure.

# DE genes in patient 1
zscores_dif_1 <- m_cent_all[m_cent_all$chip1 < -2 | m_cent_all$chip1 > 2, c(1,2)]
nrow(zscores_dif_1)

# Down regulated genes in patient 1
down_reg_1 <- data.frame(m_cent_all[m_cent_all$chip1 < -2, c(1,2)])
nrow(down_reg_1)

# Up regulated genes in patient 1
up_reg_1 <- data.frame(m_cent_all[m_cent_all$chip1 > 2, c(1,2)])
nrow(up_reg_1)

# DE genes in patient 2
zscores_dif_2 <- m_cent_all[m_cent_all$chip2 < -2 | m_cent_all$chip2 > 2, c(1,3)]
nrow(zscores_dif_2)

# Down regulated genes in patient 2
down_reg_2 <- data.frame(m_cent_all[m_cent_all$chip2 < -2, c(1,3)])
nrow(down_reg_2)

# Up regulated genes in patient 2
up_reg_2 <- data.frame(m_cent_all[m_cent_all$chip2 > 2, c(1,3)])
nrow(up_reg_2)

# DE genes in patient 3
zscores_dif_3 <- m_cent_all[m_cent_all$chip3 < -2 | m_cent_all$chip3 > 2, c(1,4)]
nrow(zscores_dif_3)

# Down regulated genes in patient 3
down_reg_3 <- data.frame(m_cent_all[m_cent_all$chip3 < -2, c(1,4)])
nrow(down_reg_3)

# Up regulated genes in patient 3
up_reg_3 <- data.frame(m_cent_all[m_cent_all$chip3 > 2, c(1,4)])
nrow(up_reg_3)

# DE genes in common in patient 1 and 2 
id_common12 <- merge(x = zscores_dif_1, y = zscores_dif_2, by = "Id")
nrow(id_common12)

# DE genes in common in patient 1 and 3 
id_common13 <- merge(x = zscores_dif_1, y = zscores_dif_3, by = "Id")
nrow(id_common13)

# DE genes in common in patient 2 and 3 
id_common23 <- merge(x = zscores_dif_2, y = zscores_dif_3, by = "Id")
nrow(id_common23)

# DE genes in common in all patients
id_common_all <- merge(x = id_common23, y = zscores_dif_1, by = "Id")
nrow(id_common_all)


#(h) Apply the Bayesian method of L¨onnstedt and Speed (package limma). 
#Establish a few lines to identify the genes with differential expression 
#and justify the cutoff considered.


rownames(m_cent_all) <- m_cent_all$Id
m_cent_all <- m_cent_all[2:4]

par(mfrow=c(1,2))

# comparison between patient 1 and 2
fit12 <- lmFit(m_cent_all[,c(1,2)])

fit12_05 <- eBayes(fit12, proportion = 0.05, stdev.coef.lim = c(-2,2))
table12_05 <- topTable(fit12_05, number=2994, adjust = "BH")
max(table12_05$B)
mean(table12_05$B)
diff12_05 <- table12_05[table12_05$B > 0,]
nrow(diff12_05)


# Volcanoplot when proportion 0.05
volcanoplot(fit12_05, highlight = nrow(diff12_05), style = "B-statistic", main = "(a)")
abline(0,0)


fit12_01 <- eBayes(fit12, proportion = 0.1, stdev.coef.lim = c(-2,2))
table12_01 <- topTable(fit12_01, number=2994, adjust = "BH")
max(table12_01$B)
mean(table12_01$B)
diff12_01 <- table12_01[table12_01$B > 0,]
nrow(diff12_01)


# Volcanoplot when proportion 0.1
volcanoplot(fit12_01, highlight = nrow(diff12_01), style = "B-statistic", main = "(b)")
abline(0,0)

# comparison between patient 1 and 3
par(mfrow=c(1,2))
fit13 <- lmFit(m_cent_all[,c(1,3)])

fit13_05 <- eBayes(fit13, proportion = 0.05, stdev.coef.lim = c(-2,2))
max(table13_05$B)
mean(table13_05$B)
table13_05 <- topTable(fit13_05, number=2994, adjust = "BH")
diff13_05 <- table13_05[table13_05$B > 0,]
nrow(diff13_05)


# Volcanoplot when proportion 0.05
volcanoplot(fit13_05, highlight = nrow(diff13_05), style = "B-statistic", main = "(a)")
abline(0,0)


fit13_01 <- eBayes(fit13, proportion = 0.1, stdev.coef.lim = c(-2,2))
table13_01 <- topTable(fit13_01, number=2994, adjust = "BH")
max(table13_01$B)
mean(table13_01$B)
diff13_01 <- table13_01[table13_01$B > 0,]
nrow(diff13_01)


# Volcanoplot when proportion 0.1
volcanoplot(fit13_01, highlight = nrow(diff13_01), style = "B-statistic", main = "(b)")
abline(0,0)


# comparison between patient 2 and 3
par(mfrow=c(1,2))
fit23 <- lmFit(m_cent_all[,c(2,3)])

fit23_05 <- eBayes(fit23, proportion = 0.05, stdev.coef.lim = c(-2,2))
table23_05 <- topTable(fit23_05, number=2994, adjust = "BH")
max(table23_05$B)
mean(table23_05$B)
diff23_05 <- table23_05[table23_05$B > 0,]
nrow(diff23_05)


# Volcanoplot when proportion 0.05
volcanoplot(fit23_05, highlight = nrow(diff23_05), style = "B-statistic", main = "(a)")
abline(0,0)


fit23_01 <- eBayes(fit23, proportion = 0.1, stdev.coef.lim = c(-2,2))
max(table23_01$B)
mean(table23_01$B)
table23_01 <- topTable(fit23_01, number=2994, adjust = "BH")
diff23_01 <- table23_01[table23_01$B > 0,]
nrow(diff23_01)


# Volcanoplot when proportion 0.1
volcanoplot(fit23_01, highlight = nrow(diff23_01), style = "B-statistic", main = "(b)")
abline(0,0)


# comparison between patient 1, 2 and 3
par(mfrow=c(1,2))
fit_all <- lmFit(m_cent_all)

fit_all_05 <- eBayes(fit_all, proportion = 0.05, stdev.coef.lim = c(-2,2))
max(table_all_05$B)
mean(table_all_05$B)
table_all_05 <- topTable(fit_all_05, number=2994, adjust = "BH")
diff_all_05 <- table_all_05[table_all_05$B > 0,]
nrow(diff_all_05)

# Volcanoplot when proportion 0.05
volcanoplot(fit_all_05, highlight = nrow(diff_all_05), style = "B-statistic", main = "(a)")
abline(0,0)


fit_all_01 <- eBayes(fit_all, proportion = 0.1, stdev.coef.lim = c(-2,2))
table_all_01 <- topTable(fit_all_01, number=2994, adjust = "BH")
max(table_all_01$B)
mean(table_all_01$B)
diff_all_01 <- table_all_01[table_all_01$B > 0,]
nrow(diff_all_01)


# Volcanoplot when proportion 0.1
volcanoplot(fit23_01, highlight = nrow(diff23_01), style = "B-statistic", main = "(b)")
abline(0,0)


# (i) Compare the results obtained by the two methods 

# Check if there are common genes between Z-score and limma methods
intersept_genes_1 <- id_common12[id_common12$Ids %in% diff12_01, ]
nrow(intersept_genes_1)

intersept_genes_2 <- id_common13[id_common13$Ids %in% diff13_01, ]
nrow(intersept_genes_2)

intersept_genes_3 <- id_common23[id_common23$Ids %in% diff23_01, ]
nrow(intersept_genes_3)

citation()
