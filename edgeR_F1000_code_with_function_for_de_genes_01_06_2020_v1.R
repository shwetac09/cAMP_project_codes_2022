source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")


library(edgeR)  
library(limma)

setwd("D:/Shweta/LAB/analysis_with_dcya_4mM_camp/logfc12pval001/")
#dir.create("trial")
#----raw data------
raw_data <- read.table("D:/Shweta/LAB/2018_05/rna_seq_dcya_strains/coding_gene_readcounts.cov")

nrow(raw_data)#4102

#separating gene id and data columns

readcounts <- raw_data[ ,-1]
nrow(readcounts)#4102|after removing ins = 4090



gene_id <- raw_data[ ,1]

rownames(readcounts) <- make.names(gene_id, unique = TRUE)

colnames(readcounts) <- paste(c(rep("cAMP_0.01mM", 2), rep("cAMP_0.05mM", 2), rep("cAMP_0.1mM", 2), rep("cAMP_0.3mM",2),
                                rep("cAMP_0.6mM",2), rep("cAMP_0.8mM", 2), rep("cAMP_1.0mM", 2), rep("cAMP_2.0mM",2), 
                                rep("cAMP_4.0mM", 2), rep("dcyaA", 2), rep("N600Y",2),
                                rep("N600Y_4mM",2), rep("WT",2)), c(rep(1:2, 13)), sep = "")
View(readcounts)


readcounts <- readcounts[ ,-c(21:24)] #removing n600y data

#saving coverage files for submission to GEO
# readcounts$genename <- rownames(readcounts)
# View(readcounts)
# readcounts <- readcounts[ ,c(23,1:22)]
# write.table(readcounts,
#             file = "D:/Shweta/LAB/final_results_v1/CPM/cyaA_readcounts.table",
#             sep = "\t", row.names = F)


group <- c(rep("cAMP_0.01mM", 2), rep("cAMP_0.05mM", 2), rep("cAMP_0.1mM", 2), rep("cAMP_0.3mM",2),
           rep("cAMP_0.6mM",2), rep("cAMP_0.8mM", 2), rep("cAMP_1.0mM", 2), rep("cAMP_2.0mM",2), 
           rep("cAMP_4.0mM", 2), rep("dcyaA", 2),rep("WT",2))


group <- factor(group)

table(group)

View(readcounts)
readcounts[duplicated(rownames(readcounts)), ]
#nrow(readcounts2)
# readcounts[which(readcounts$WT1 > 5e4), ]
# readcounts <- readcounts[-which(readcounts$WT1 > 5e4), ] #removing genes with high reads lpp and cspA/E
nrow(readcounts)
#readcounts <- readcounts2
#View(readcounts2)


#---- minor check for replicate variation----

for(i in seq(1, ncol(readcounts),2)){
  
  x <- readcounts[ ,i]
  y <- readcounts[ ,(i+1)]
  plot(x,y, xlab = colnames(readcounts)[i], ylab = colnames(readcounts)[(i+1)])#,
       #xlim = c(0, 1e05), ylim = c(0,1e05))
  abline(lm(y~x))
  legend("bottomright",c("result"), c(signif(cor.test(x,y)$estimate)), cex = 0.7)
  #summary(lm(y ~ x))
  }



#----creating DGElist ---------


y <- DGEList(readcounts, group = group, genes = rownames(readcounts))
y$samples
#filtering low counts: keeping genes which have more than 0.5CPM in atleast 2 libraries 
10/1.3

keep <- rowSums(cpm(y) > 7.7) >=2 # 7.7|

table(keep) #3490 genes pass the criterion| 3467 pass this criterion | 3478 after removing ins

y <- y[keep, ,keep.lib.sizes = FALSE]
nrow(y)
head(y)

#-------normalising for the RNA composition by TMM---------
View(y$counts)

y <- calcNormFactors(y)

options(digits = 3)

y$samples
View(y$counts)
View(y)


#----plotting mds plots-----

pch <- c(0,1,2,4,5,6,15,18,17,10,7,8,14)
colors <- rep(c(1:13),2)
plotMDS(y, col = colors[group], pch = pch[group], cex = 1.2)
legend("bottomleft", legend = levels(group), pch = pch, col = colors, ncol = 2, cex = .75)

#plotting md plot 

plotMD(y, column = 22)
abline(h=0, col = "red", lty =2, lwd =2)

#----specification of design for linear model. -----

design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

#---estimating dispersion -----

y <- estimateDisp(y, design, robust = TRUE)
plotBCV(y)

#genewise dispersios 
fit <- glmQLFit(y, design, robust = TRUE)

head(fit$coefficients)
plotQLDisp(fit)

summary(fit$df.prior)
#________________________________

#----to get the logCPM values -----
setwd("D:/Shweta/LAB/final_results_v1/CPM/")
View(y)

CPM <- cpm(y, log = TRUE)
View(CPM)
rownames(CPM) <- y$genes$genes
colnames(CPM) <- paste(y$samples$group, 1:2, sep = "_")
head(CPM)
filepath <- "D:/Shweta/LAB/final_results_v1/CPM/"
write.table(CPM, file = paste(filepath, "cpm_values_all_genes.table", sep = ""),
            sep = "\t")
write.table(CPM, file = paste(filepath, "logcpm_values_all_genes.table", sep = ""),
            sep = "\t")
            


#cpm_tab <- sapply(seq(1, ncol(CPM),2), function(x) rowMeans(CPM[ ,x:(x+1)]))
#head(cpm_tab)

#names_cpm_tab <- colnames(CPM)[seq(1, ncol(CPM),2)]
#names_cpm_tab

#colnames(cpm_tab) <- names_cpm_tab
filepath <- "D:/Shweta/LAB/final_results_v1/CPM/"
write.table(cpm_tab, file = paste(filepath, "cpm_values_all_genes.table", sep = ""),
            sep = "\t")


#---post normalisation replicate test again-----
for(i in seq(1,ncol(CPM),2)){
  x <- CPM[ ,i]
  y1 <- CPM[ ,(i+1)]
  plot(x,y1, xlab =  colnames(CPM)[i], ylab = colnames(CPM)[i+1])
  abline(lm(y1~x))
  print(c(colnames(CPM)[i], cor(x,y1)))
}
#_______________________________________________

filepath <- "/Shweta/LAB/final_results_v1/"
write.table(cpm_tab, file = paste(filepath, "cpm_values_all_genes.table", sep = ""),
            sep = "\t")

head(log2(test))
View(test)

logCPM <- cpm(y, prior.count = 2, log = TRUE)
View(logCPM)
rownames(logCPM) <- y$genes$genes
colnames(logCPM) <- paste(y$samples$group, 1:2, sep = "_")

View(logCPM)
write.table(logCPM, file = paste(filepath, "logCPM_all_genes.table"), sep = "\t")

#_________________________________________________________

#----Findind DE genes : wt vs dcyaA exploration --------

wt_dcya <- makeContrasts(WT - dcyaA, levels = design)

res <- glmQLFTest(fit, contrast = wt_dcya)
View(res)

#woTREAT <- topTags(res, sort.by = "logFC", n = 4000)$table

is.de <- decideTestsDGE(res)

summary(is.de)


plotMD(res, status = is.de, values = c(1, -1), col = c("red", "blue"), legend = "topright")

tr <- glmTreat(fit, contrast = wt_dcya, lfc = log2(1))
nrow(tr$table)
wt1 <- tr

#wTREAT <- topTags(tr, sort.by = "logFC", n = 4000)$table

#tabwoTR <- woTREAT[woTREAT$logFC > 1, ]
#nrow(tabwoTR)
#tabwTR <- wTREAT[wTREAT$logFC > 1, ]

# plot(tabwoTR$PValue, tabwTR$PValue)
# plot(tabwoTR$FDR, tabwTR$FDR)
# plot(tabwTR$PValue, tabwTR$FDR)


is.de <- decideTestsDGE(tr)
b <- topTags(tr, n = 5000, adjust.method = "fdr")
View(b$table)
summary(is.de)

plotMD(tr, status = is.de, values=c(1,-1), col = c("red", "blue"), legend = "topright")

View(tr$table)
#to plot volcano plot
x <- tr$table$logFC
z <- -log10(tr$table$PValue)

plot(x,z, xlab = "logFC", ylab = "-log10(Pvalue)", main = "WT vs dCyaA")
View(tr$table)


#DE genes
pval_cut <- 0.01
logfc_cut <- 1  

#+vely regulated genes
wt_dcya_pos <- tr[tr$table$PValue < pval_cut & tr$table$logFC > logfc_cut, ]$table 
nrow(wt_dcya_pos)#365
wt1 <- wt_dcya_pos
View(wt_dcya_pos)
points(wt_dcya_pos$logFC, -log10(wt_dcya_pos$PValue), col = "red")#plotting on volcano plot

# - ve regulated genes 
wt_dcya_neg <- tr[tr$table$PValue < pval_cut & tr$table$logFC < -logfc_cut, ]$table

nrow(wt_dcya_neg)

#----DE gene calculator function-----


DE_gene_calculator <- function(match, glmtr_cutoff){
  res <- glmQLFTest(fit, contrast = match)
  #topTags(res, sort.by = "logFC")
  is.de <- decideTestsDGE(res)
  #plotMD(res, status = is.de, values = c(-1,1), col = c("red", "blue"), legend = "topright")
  tr <- glmTreat(fit, contrast = match, lfc = log2(glmtr_cutoff))
  is.de <- decideTestsDGE(tr)
  #plotMD(tr, status = is.de, values = c(1,-1), col = c("red", "blue"), legend = "topright")
  return(tr$table)
}

#--------WRT dcyaA---------

wt_dcya<- makeContrasts(WT - dcyaA, levels = design)
wt_dcya_tr <- DE_gene_calculator(wt_dcya, glmtr_cutoff = 1)
nrow(wt_dcya_tr)

camp0.01_dcya<- makeContrasts(cAMP_0.01mM - dcyaA, levels = design)
camp0.01_dcya_tr <- DE_gene_calculator(camp0.01_dcya, glmtr_cutoff = 1)

camp0.05_dcya<- makeContrasts(cAMP_0.05mM - dcyaA, levels = design)
camp0.05_dcya_tr <- DE_gene_calculator(camp0.05_dcya, glmtr_cutoff = 1)

camp0.1_dcya<- makeContrasts(cAMP_0.1mM - dcyaA, levels = design)
camp0.1_dcya_tr <- DE_gene_calculator(camp0.1_dcya, glmtr_cutoff = 1)

camp0.3_dcya<- makeContrasts(cAMP_0.3mM - dcyaA, levels = design)
camp0.3_dcya_tr <- DE_gene_calculator(camp0.3_dcya, glmtr_cutoff = 1)

camp0.6_dcya<- makeContrasts(cAMP_0.6mM - dcyaA, levels = design)
camp0.6_dcya_tr <- DE_gene_calculator(camp0.6_dcya, glmtr_cutoff = 1)

camp0.8_dcya<- makeContrasts(cAMP_0.8mM - dcyaA, levels = design)
camp0.8_dcya_tr <- DE_gene_calculator(camp0.8_dcya, glmtr_cutoff = 1)

camp1.0_dcya<- makeContrasts(cAMP_1.0mM - dcyaA, levels = design)
camp1.0_dcya_tr <- DE_gene_calculator(camp1.0_dcya, glmtr_cutoff = 1)

camp2.0_dcya<- makeContrasts(cAMP_2.0mM - dcyaA, levels = design)
camp2.0_dcya_tr <- DE_gene_calculator(camp2.0_dcya, glmtr_cutoff = 1)

camp4.0_dcya<- makeContrasts(cAMP_4.0mM - dcyaA, levels = design)
camp4.0_dcya_tr <- DE_gene_calculator(camp4.0_dcya, glmtr_cutoff = 1)


tstab_wrt_dcya <- list(camp0.01_dcya_tr, camp0.05_dcya_tr, camp0.1_dcya_tr, camp0.3_dcya_tr,
                       camp0.6_dcya_tr, camp0.8_dcya_tr, camp1.0_dcya_tr, camp2.0_dcya_tr,
                       camp4.0_dcya_tr, wt_dcya_tr)

names(tstab_wrt_dcya) <- c("0.01", "0.05", "0.1", "0.3", "0.6", "0.8", "1.0",
                           "2.0", "4.0", "WT" )

View(tstab_wrt_dcya)


DE_wrt_dcya_all <- do.call(cbind.data.frame , tstab_wrt_dcya)
cols_req <- grep("[[:digit:]].logFC|PValue|WT.logFC|WT.P", colnames(DE_wrt_dcya_all))

logfc_dcya_table <- DE_wrt_dcya_all[ ,cols_req]

View(logfc_dcya_table)

filepath <- "D:/Shweta/LAB/final_results_v1/rna_seq_checks/rna_seq_without_removing_genes/"
write.table(logfc_dcya_table, file = paste(filepath, "logfc_pval_all_genes.table", sep = ""), sep = "\t")

pval_cut <- 0.01
logfc_cut <- 1

pos_genes_dt <- logfc_dcya_table[logfc_dcya_table$WT.logFC > 1 & logfc_dcya_table$WT.PValue< 0.01, ]
nrow(pos_genes_dt)#365

plot(logfc_dcya_table$`2.0.logFC`, -log10(logfc_dcya_table$`2.0.PValue`))
pos_genes_logfc_2 <- logfc_dcya_table[logfc_dcya_table$`2.0.PValue` < pval_cut & logfc_dcya_table$`2.0.logFC`> logfc_cut, ]

pos_pos_2diff <- subset(pos_genes_logfc_2, !(rownames(pos_genes_logfc_2) %in% cmn))
View(pos_pos_2diff)

plot(pos_pos_2diff$WT.logFC, pos_pos_2diff$`2.0.logFC`, ylim = c(0,1.5))

nrow(pos_genes_logfc_2)
neg_genes_logfc_2 <- logfc_dcya_table[logfc_dcya_table$`2.0.PValue` < pval_cut & logfc_dcya_table$`2.0.logFC`< - logfc_cut, ]
nrow(neg_genes_logfc_2)
points(neg_genes_logfc_2$`2.0.logFC`, -log10(neg_genes_logfc_2$`2.0.PValue`), col = "red")

points(pos_genes_logfc_2$`2.0.logFC`, -log10(pos_genes_logfc_2$`2.0.PValue`), col = "red")

pos_genes_logfc_wt <- logfc_dcya_table[logfc_dcya_table$WT.PValue < pval_cut & logfc_dcya_table$WT.logFC > logfc_cut, ]
nrow(pos_genes_logfc_wt)
neg_genes_logfc_wt<- logfc_dcya_table[logfc_dcya_table$WT.PValue < pval_cut & logfc_dcya_table$WT.logFC < -logfc_cut, ]
nrow(neg_genes_logfc_wt)

cmn <- intersect(rownames(pos_genes_logfc_2), rownames(pos_genes_logfc_wt))
length(cmn)

de_genes_wt_cya <- logfc_dcya_table[logfc_dcya_table$WT.PValue < 0.05, ]

cmn <- intersect(rownames(pos_genes_logfc_wt), rownames(pos_genes_logfc_2))
length(cmn)

only_wt <- subset(pos_genes_logfc_wt, !(rownames(pos_genes_logfc_wt) %in% cmn))
nrow(only_wt)

only_2 <- subset(pos_genes_logfc_2, !(rownames(pos_genes_logfc_2) %in% cmn))
View(only_2)

cmn_tab <- subset(logfc_dcya_table, rownames(logfc_dcya_table) %in% cmn)
nrow(cmn_tab)

plot(pos_genes_logfc_wt$WT.logFC, pos_genes_logfc_wt$`2.0.logFC`)
points(pos_genes_logfc_2$WT.logFC, pos_genes_logfc_2$`2.0.logFC`, col = "blue")
points(cmn_tab$WT.logFC, cmn_tab$`2.0.logFC`, col = "red")



plot(logfc_dcya_table$`2.0.logFC`, -log10(logfc_dcya_table$`2.0.PValue`))
points(pos_genes_logfc$`2.0.logFC`, -log10(pos_genes_logfc$`2.0.PValue`), col = "red")




View(pos_genes_logfc)
nrow(pos_genes_logfc)

getwd()
write.table(pos_genes_logfc, file = paste(filepath, "logfc_pos_genes.table", sep = ""), sep = "\t")


set1 <- as.numeric(as.character(pos_genes_logfc$X4.0.logFC))     
View(set1)

set2 <- as.numeric(as.character(pos_genes_logfc$X2.0.logFC))     
View(set1)

plot(density(set1), xlim = c(-4, 10))
lines(density(set2), col = "red")
t.test(set1, set2)


plot(pos_genes_logfc$`4.0.logFC`, pos_genes_logfc$WT.logFC)
abline(h = 1.5, v = 1)

diff <- pos_genes_logfc[pos_genes_logfc$`4.0.logFC` < 1 & pos_genes_logfc$WT.logFC > 1.5, ]
View(diff)
nrow(diff)




#making a file with strong cut offs 
logfc_dcya_table <- read.table("logfc_pval_all_genes.table")

pval_cut <- 0.05
logfc_cut <- 1.5

pos_genes_logfc <- logfc_dcya_table[logfc_dcya_table$WT.PValue < pval_cut & logfc_dcya_table$WT.logFC > logfc_cut, ]
nrow(pos_genes_logfc)


write.table(pos_genes_logfc, file = paste(filepath, "logfc_pos_genes1.5.table", sep = ""), sep = "\t")




#_____________________________________________________________________

#-------- WRT WT ----------
dcya_wt<- makeContrasts(dcyaA - WT, levels = design)
dcya_wt_tr <- DE_gene_calculator(dcya_wt, glmtr_cutoff = 1)

camp0.01_wt<- makeContrasts(cAMP_0.01mM - WT, levels = design)
camp0.01_wt_tr <- DE_gene_calculator(dcya_wt, glmtr_cutoff = 1)

camp0.05_wt<- makeContrasts(cAMP_0.05mM - WT, levels = design)
camp0.05_wt_tr <- DE_gene_calculator(camp0.05_wt, glmtr_cutoff = 1)

camp0.1_wt<- makeContrasts(cAMP_0.1mM - WT, levels = design)
camp0.1_wt_tr <- DE_gene_calculator(camp0.1_wt, glmtr_cutoff = 1)

camp0.3_wt<- makeContrasts(cAMP_0.3mM - WT, levels = design)
camp0.3_wt_tr <- DE_gene_calculator(camp0.3_wt, glmtr_cutoff = 1)

camp0.6_wt<- makeContrasts(cAMP_0.6mM - WT, levels = design)
camp0.6_wt_tr <- DE_gene_calculator(camp0.6_wt, glmtr_cutoff = 1)

camp0.8_wt<- makeContrasts(cAMP_0.8mM - WT, levels = design)
camp0.8_wt_tr <- DE_gene_calculator(camp0.8_wt, glmtr_cutoff = 1)

camp1.0_wt<- makeContrasts(cAMP_1.0mM - WT, levels = design)
camp1.0_wt_tr <- DE_gene_calculator(camp1.0_wt, glmtr_cutoff = 1)


camp2.0_wt<- makeContrasts(cAMP_2.0mM - WT, levels = design)
camp2.0_wt_tr <- DE_gene_calculator(camp2.0_wt, glmtr_cutoff = 1)


camp4.0_wt<- makeContrasts(cAMP_4.0mM - WT, levels = design)
camp4.0_wt_tr <- DE_gene_calculator(camp4.0_wt, glmtr_cutoff = 1)

tstab_wrt_wt <- list(camp0.01_wt_tr, camp0.05_wt_tr, camp0.1_wt_tr, camp0.3_wt_tr,
                       camp0.6_wt_tr, camp0.8_wt_tr, camp1.0_wt_tr, camp2.0_wt_tr,
                       camp4.0_wt_tr, dcya_wt_tr)

names(tstab_wrt_wt) <- c("0.01", "0.05", "0.1", "0.3", "0.6", "0.8", "1.0",
                           "2.0", "4.0", "dcya" )

View(tstab_wrt_wt)


DE_wrt_wt_all <- do.call(cbind.data.frame , tstab_wrt_wt)
View(DE_wrt_wt_all)
cols_req <- grep("[[:digit:]].logFC|PValue|dcya.logFC|dcya.P", colnames(DE_wrt_wt_all))

logfc_wt_table <- DE_wrt_wt_all[ ,cols_req]

View(logfc_wt_table)
logfc_wt_table$genename <- rownames(logfc_wt_table)
write.table(logfc_wt_table, 
            file = "D:/Shweta/LAB/final_results_v1/All_tables_lists/all_logfc_pval_wrt_WT.table",
            sep = "\t")


plot(logfc_wt_table$`2.0.logFC`, -log10(logfc_wt_table$`2.0.PValue`))

pval_cut <- 0.05
logfc_cut <- 1

pos_genes_wt_logfc <- logfc_wt_table[logfc_wt_table$`2.0.PValue` < pval_cut & logfc_wt_table$`2.0.logFC` > logfc_cut, ]
nrow(pos_genes_wt_logfc)
neg_genes_wt_logfc <- logfc_wt_table[logfc_wt_table$`2.0.PValue` < pval_cut & logfc_wt_table$`2.0.logFC` < -logfc_cut, ]

nrow(neg_genes_wt_logfc) #70 fit hills 
points(neg_genes_wt_logfc$`2.0.logFC`, -log10(neg_genes_wt_logfc$`2.0.PValue`), col = "red")
View(pos_genes_wt_logfc)

de_genes_2vsWT <- logfc_wt_table[logfc_wt_table$`2.0.PValue` < 0.05, ]
nrow(de_genes_2vsWT)




#----
pval_cut <- 0.05
logfc_cut <- 1.5

pos_genes_wt_logfc <- logfc_wt_table[logfc_wt_table$WT.PValue < pval_cut & logfc_wt_table$WT.logFC < -logfc_cut, ]
nrow(pos_genes_logfc)

logfc_wt_table
#plotting logfc vs pvalue cut offs for all cAMP concentrations

plot(pos_genes_wt_logfc$`0.3.logFC`, -log10(pos_genes_wt_logfc$`0.3.PValue`))
abline(v = -1.5, h = -log10(0.05))
plot(pos_genes_wt_logfc$`0.6.logFC`, -log10(pos_genes_wt_logfc$`0.6.PValue`))
abline(v = -1.5, h = -log10(0.05))
plot(pos_genes_wt_logfc$`0.8.logFC`, -log10(pos_genes_wt_logfc$`0.8.PValue`), col = "blue")
abline(v = -1.5, h = -log10(0.05))
plot(pos_genes_wt_logfc$`1.0.logFC`, -log10(pos_genes_wt_logfc$`1.0.PValue`), col = "yellow")
abline(v = -1.5, h = -log10(0.05))
plot(pos_genes_wt_logfc$`2.0.logFC`, -log10(pos_genes_wt_logfc$`2.0.PValue`), col = "red")
abline(v = -1.5, h = -log10(0.05))
points(pos_genes_wt_logfc$`4.0.logFC`, -log10(pos_genes_wt_logfc$`4.0.PValue`), col = "magenta")
View(pos_genes_wt_logfc)
abline(v = -1.5, h = -log10(0.05))

plot(logfc_wt_table$WT.logFC, -log10(logfc_wt_table$WT.PValue))
plot(logfc_wt_table$`0.6.logFC`, -log10(logfc_wt_table$`0.6.PValue`))
plot(logfc_wt_table$`0.8.logFC`, -log10(logfc_wt_table$`0.8.PValue`))
plot(logfc_wt_table$`1.0.logFC`, -log10(logfc_wt_table$`1.0.PValue`))
plot(logfc_wt_table$`2.0.logFC`, -log10(logfc_wt_table$`2.0.PValue`))
plot(logfc_wt_table$`4.0.logFC`, -log10(logfc_wt_table$`4.0.PValue`))


#----dcya- 4mM comparison, exploratory----
dcya_4 <- makeContrasts(cAMP_4.0mM - dcyaA, levels = design)


res <- glmQLFTest(fit, contrast = dcya_4)

topTags(res, sort.by = "logFC")

is.de <- decideTestsDGE(res)

summary(is.de)

plotMD(res, status = is.de, values = c(1, -1), col = c("red", "blue"), legend = "topright")

tr <- glmTreat(fit, contrast = dcya_4, lfc = log2(1.2))
topTags(tr, sort.by = "logFC")

is.de <- decideTestsDGE(tr)
summary(is.de)
plotMD(tr, status = is.de, values=c(1,-1), col = c("red", "blue"), legend = "topright")

#to plot volcano plot
x <- tr$table$logFC
z <- -log10(tr$table$PValue)

plot(x,z, xlab = "logFC", ylab = "-log10(Pvalue)", main = "dcya vs 4mM cAMP")
#View(tr)

pval_co <- 0.05
logfc_co <- 1.2
  
pos_g_dcya_4camp <-subset(tr$table, tr$table$PValue < pval_co & tr$table$logFC > logfc_co)
nrow(pos_g_dcya_4camp)
#View(pos_g_dcya_4camp)
pos_g_4 <- rownames(pos_g_dcya_4camp)
head(pos_g_4)
neg_g_dcya_4camp <-subset(tr$table, tr$table$PValue < pval_co & tr$table$logFC < -logfc_co)
nrow(neg_g_dcya_4camp)

write.table(de_genes_dcya_4camp,
            file = "D:\Shweta/LAB/analysis_with_dcya_4mM_camp/all_DE_genes_dcya_4camp.table")




##================================================================================



dcya_4 <- makeContrasts(dcyaA - cAMP_4.0mM, levels = design)
dcya_4_tr <- DE_gene_calculator(dcya_4)
View(dcya_4_tr)
nrow(dcya_4_tr)
tmp <- dcya_4_tr[dcya_4_tr$PValue<0.05, ]
nrow(tmp)

camp0.01_4 <- makeContrasts(cAMP_0.01mM - cAMP_4.0mM, levels = design)
camp0.01_4_tr <- DE_gene_calculator(camp0.01_4)
View(camp0.01_4_tr)

camp0.05_4 <- makeContrasts(cAMP_0.05mM - cAMP_4.0mM, levels = design)
camp0.05_4_tr <- DE_gene_calculator(camp0.05_4)

camp0.1_4 <- makeContrasts(cAMP_0.1mM - cAMP_4.0mM, levels = design)
camp0.1_4_tr <- DE_gene_calculator(camp0.1_4)

camp0.3_4 <- makeContrasts(cAMP_0.3mM - cAMP_4.0mM, levels = design)
camp0.3_4_tr <- DE_gene_calculator(camp0.3_4)

camp0.6_4 <- makeContrasts(cAMP_0.6mM - cAMP_4.0mM, levels = design)
camp0.6_4_tr <- DE_gene_calculator(camp0.6_4)

camp0.8_4 <- makeContrasts(cAMP_0.8mM - cAMP_4.0mM, levels = design)
camp0.8_4_tr <- DE_gene_calculator(camp0.8_4)

camp1.0_4 <- makeContrasts(cAMP_1.0mM - cAMP_4.0mM, levels = design)
camp1.0_4_tr <- DE_gene_calculator(camp1.0_4)

camp2.0_4 <- makeContrasts(cAMP_2.0mM - cAMP_4.0mM, levels = design)
camp2.0_4_tr <- DE_gene_calculator(camp2.0_4)


###=========================================================================



##===================================================




#-------heirarchical expression---------- 

all_genes_all_conc <- list(camp0.01_dcya_tr, camp0.05_dcya_tr, camp0.1_dcya_tr, camp0.3_dcya_tr, camp0.6_dcya_tr,
                           camp0.8_dcya_tr, camp1.0_dcya_tr, camp2.0_dcya_tr, camp4.0_dcya_tr, wt_dcya_tr)


View(all_genes_all_conc)
de_gene_list <- vector("list", 10)

for(i_list in 1:length(all_genes_all_conc)){
  de_genes <- subset(all_genes_all_conc[[i_list]], all_genes_all_conc[[i_list]]$PValue < 0.01)
  de_gene_list[[i_list]] <- de_genes
}
View(de_gene_list)



## WITH A MORE RELAXED CUT OFF DO I GET A BETTER GRASP: P <0.02
de_gene_list_0.02 <- vector("list", 10)
for(i_list in 1:length(all_genes_all_conc)){
  de_genes <- subset(all_genes_all_conc[[i_list]], all_genes_all_conc[[i_list]]$PValue < 0.02)
  de_gene_list_0.02[[i_list]] <- de_genes
}
View(de_gene_list_0.02)




##to determine the overlap within these genes 


merged_tab <- data.frame(rownames(de_gene_list[[10]]))
rownames(merged_tab) <- rownames(de_gene_list[[10]])
for (i in 1:9){
  name2 <- data.frame(rownames(de_gene_list[[i]]))
  rownames(name2) <- rownames(de_gene_list[[i]])
  
  merged_tab <- merge(merged_tab, name2, by = "row.names", all.x = TRUE, all.y = TRUE)
  rownames(merged_tab) <- merged_tab$Row.names
  merged_tab <- merged_tab[ ,-1]
  
}

write.table(merged_tab, file = "merged_tab_trial.table", sep = "\t")
nrow(merged_tab)

##============================
##to write the table starting from 0.01

merged2_tab <- data.frame(rownames(de_gene_list[[1]]))
rownames(merged2_tab) <- rownames(de_gene_list[[1]])
for (i in 2:9){
  name2 <- data.frame(rownames(de_gene_list[[i]]))
  rownames(name2) <- rownames(de_gene_list[[i]])
  
  merged2_tab <- merge(merged2_tab, name2, by = "row.names", all.x = TRUE, all.y = TRUE)
  rownames(merged2_tab) <- merged2_tab$Row.names
  merged2_tab <- merged2_tab[ ,-1]
  
}

write.table(merged2_tab, file = "merged_tab_trial2.table", sep = "\t")
nrow(merged2_tab)
View(merged2_tab)
colnames(merged2_tab) <- c("0.01", "0.05", "0.1", "0.3", "0.6", "0.8", "1.0", "2.0", "4.0")


##to count stuff in merged2_tab

#to count freq of genes per row. if its present in all 5 (0.6 to 4mM) the score ll be 5 and so on. 
foo <- as.matrix.data.frame(merged2_tab)
View(foo)

foo[is.na(foo)] <- 0

foo[ foo != 0] <- 1

row_sum <- NULL
for(i in 1:nrow(foo)){
  sums <- sum(as.numeric(foo[i, ]))
  row_sum <- c(row_sum, sums)
}

boo <- data.frame(row_sum)
View(boo)
factor(boo$row_sum)

f <- table(boo$row_sum)


##To calculate number of genes in each column, ie, DE genes per camp conc 
col_sum <- NULL
for(i in 1:ncol(foo)){
  sums <- sum(as.numeric(foo[ ,i]))
  col_sum <- c(col_sum, sums)
}
col_sum

##-------------------
##To find heirarchy across camp conc. the result of this is written in lab note boook and excel sheet
common_6 <- Reduce(intersect, list(merged2_tab$`0.3`,merged2_tab$`0.6`, merged2_tab$`0.8`, merged2_tab$`1.0`,
                                   merged2_tab$`2.0`, merged2_tab$`4.0`))

common_5 <- Reduce(intersect, list(merged2_tab$`0.6`, merged2_tab$`0.8`, merged2_tab$`1.0`,
                                   merged2_tab$`2.0`, merged2_tab$`4.0`))
length(common_5)
View(common_5)

common_4 <- Reduce(intersect, list(merged2_tab$`0.8`, merged2_tab$`1.0`,
                                   merged2_tab$`2.0`, merged2_tab$`4.0`))
common_4

length(common_4)

length(common_4[is.na(common_4)])

common_3 <- Reduce(intersect, list(merged2_tab$`1.0`,
                                   merged2_tab$`2.0`, merged2_tab$`4.0`))
common_3

length(common_3)
length(common_3[is.na(common_3)])

common_2<- Reduce(intersect, list(merged2_tab$`2.0`, merged2_tab$`4.0`))

length(common_2)

common_2[is.na(common_2)])

length(intersect(common_2, merged2_tab$`4.0`))

length(setdiff(merged2_tab$`4.0`,intersect(common_5, common_4)))

####=======================================================

##to see common genes per conc common with the wt vs dcya set 
View(de_gene_list)

cmn_de_genes <- vector("list", 9)

for(i in 1:9){
  cmn_de_genes[[i]] <- intersect(rownames(de_gene_list[[i]]), rownames(de_gene_list[[10]]))
  
}











