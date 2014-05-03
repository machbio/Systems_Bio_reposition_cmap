source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
biocLite("simpleaffy")
biocLite("affyPLM")
biocLite("affy")
biocLite("gcrma")
biocLite("hgu133plus2.db")
biocLite("hgu133a.db")
biocLite("limma")

library(GEOquery)
library(affy)
library(gcrma)
library(simpleaffy)

# GSE20437 Healthy Individuals
getGEOSuppFiles("GSE20437")
untar("GSE20437/GSE20437_RAW.tar", exdir="healthy_data")
cels <- list.files("healthy_data/", pattern = "[gz]")
sapply(paste("healthy_data", cels, sep="/"), gunzip)

# GSE31519 Cancer Individuals
getGEOSuppFiles("GSE31519")
untar("GSE31519/GSE31519_RAW.tar", exdir="diseased_data")
cels <- list.files("diseased_data/", pattern = "[gz]")
sapply(paste("diseased_data", cels, sep="/"), gunzip)


celfiles <- read.affy(covdesc="phenodata.txt", path="data")
celfiles.gcrma <- gcrma(celfiles)

library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
boxplot(celfiles, col=cols)
library(affyPLM)
boxplot(celfiles.gcrma, col=cols)
hist(celfiles, col=cols)
hist(celfiles.gcrma, col=cols)

celfiles.qc <- fitPLM(celfiles)
image(celfiles.qc, which=1, add.legend=TRUE)
image(celfiles.qc, which=4, add.legend=TRUE)
RLE(celfiles.qc, main="RLE")
NUSE(celfiles.qc, main = "NUSE")

#Removing Outliers
#nusevalues <- NUSE(celfiles.qc, type="stats")
#celbols <- nusevalues[1,] > 1.05
#celbols[celbols==TRUE]

eset <- exprs(celfiles.gcrma)
eset=format(eset, digits=5)
probes=row.names(eset)

distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters)

celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)
celfiles.filtered$filter.log

samples <- celfiles.gcrma$Target
samples <- as.factor(samples)
design <- model.matrix(~0 + samples)
colnames(design) <- c("diseased","healthy")
#design

library(limma)
fit <- lmFit(exprs(celfiles.filtered$eset), design)
contrast.matrix <- makeContrasts(diseased_healthy = diseased - healthy, levels=design)
#contrast.matrix
huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)
topTable(huvec_ebFit, number=10, coef=1)

nrow(topTable(huvec_ebFit, coef=1, number=10000, lfc=2))
probeset.list <- topTable(huvec_ebFit, coef=1, number=10000, lfc=2)
#length(probeset.list)

library(hgu133a.db)
library(annotate)

#head(probeset.list)
#attributes(probeset.list)

probematrix = as.matrix(probeset.list)
probenames = row.names(probematrix)
gene.symbols <- getSYMBOL(probenames, "hgu133a")

####
#gene.symbols <- getSYMBOL(probeset.list$ID, "hgu133a")

results <- cbind(probeset.list, gene.symbols)
head(results)

write.table(results, "results.txt", sep="\t", quote=FALSE)


#########

contrast.matrix <- makeContrasts(diseased_control = diseased - control, levels=design)

tnbc_fits <- contrasts.fit(fit, contrast.matrix)
tnbc_ebFit <- eBayes(tnbc_fits)
topTable(tnbc_ebFit, number=10, coef=1)


#nrow(topTable(tnbc_ebFit, coef=1, number=10000, lfc=5))
#nrow(topTable(tnbc_ebFit, coef=1, number=10000, lfc=4))
#nrow(topTable(tnbc_ebFit, coef=1, number=10000, lfc=3))
#nrow(topTable(tnbc_ebFit, coef=1, number=10000, lfc=2))

