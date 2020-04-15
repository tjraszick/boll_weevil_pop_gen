## load the vcfR package and adegenet packages
library("vcfR")
library("adegenet")
library("genepop")
library("poppr")
library("ggplot2")

## import the vcf file of interest
vcf_file <- "./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose=FALSE)

## convert vcfR to genlight object
BW_dDocent_Dovetail_phylo_genlight <- vcfR2genlight(vcf, n.cores = 12)

popmap <- read.table("./BW_dDocent_Dovetail_popmap.txt", sep="\t")
pop_vector <- as.vector(popmap[,2])

## append to the genlight object
pop(BW_dDocent_Dovetail_phylo_genlight) <- as.factor(c(pop_vector))

## Load the file of interest and check the basic statistics about the file
BW_dDocent_Dovetail_phylo_genlight
nLoc(BW_dDocent_Dovetail_phylo_genlight)
pop(BW_dDocent_Dovetail_phylo_genlight)
table(pop(BW_dDocent_Dovetail_phylo_genlight))
indNames(BW_dDocent_Dovetail_phylo_genlight)

## run a PCA on the genlight object
BW_dDocent_Dovetail_phylo_light_PCA <- glPca(BW_dDocent_Dovetail_phylo_genlight)
## Select the number of axes: 6

## export PCA scores and population info
BW_dDocent_Dovetail_phylo_light_PCA_scores <- BW_dDocent_Dovetail_phylo_light_PCA$scores[,c("PC1","PC2")]
Collection <- pop_vector
BW_dDocent_Dovetail_phylo_light_PCA_scores <- cbind(BW_dDocent_Dovetail_phylo_light_PCA_scores,Collection)
BW_dDocent_Dovetail_phylo_light_PCA_scores

## the scatterplot
metadata <- read.table("./BW_dDocent_Dovetail_metadata.txt", sep=",", header = TRUE)
metadata

ggplot(as.data.frame(BW_dDocent_Dovetail_phylo_light_PCA_scores), aes(x=as.numeric(BW_dDocent_Dovetail_phylo_light_PCA_scores[,1]), y=as.numeric(BW_dDocent_Dovetail_phylo_light_PCA_scores[,2]))) + geom_point(aes(color = as.factor(metadata$Population), shape = as.factor(metadata$Year)), size = 3) + scale_shape_manual(values = c(16, 5, 3))

## find groups
grp <- find.clusters(BW_dDocent_Dovetail_phylo_genlight, max.n.clust=40, n.iter = 1e9, n.start = 1e4)
## Keep 300 PCs, 6 clusters

## evaluate the clusters
head(grp$Kstat, 20)
grp$stat
grp$size
table(pop(BW_dDocent_Dovetail_phylo_genlight), grp$grp)
pops <- levels(popmap$V2)
table.value(table(pop(BW_dDocent_Dovetail_phylo_genlight), grp$grp), col.lab=paste("Group", 1:6), row.lab=paste(pops))

## run the DAPC
BW_dDocent_Dovetail_phylo_DAPC <- dapc(BW_dDocent_Dovetail_phylo_genlight, grp$grp)
## Keep 6 PCs, 2 DFs

myCol <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
group_labels <- c("AWC", "SMX", "CMX", "ARG", "RGV-Tex (2016)", "RGV")

scatter(BW_dDocent_Dovetail_phylo_DAPC, scree.da=TRUE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4, cex=3,clab=0, leg=TRUE, posi.leg="topleft", posi.da="bottomleft", txt.leg=paste(group_labels))

scatter(BW_dDocent_Dovetail_phylo_DAPC,1,1, col=myCol, bg="white", scree.da=FALSE, leg=TRUE, posi.leg="topleft", txt.leg=paste(group_labels), solid=.4)

table.value(table(pop(BW_dDocent_Dovetail_phylo_genlight), grp$grp), col.lab=paste(group_labels), row.lab=paste(pops))

## run an analysis of molecular variance (AMOVA) create input csv file on your own
stratalist <- read.csv("./BW_dDocent_Dovetail_phylo_metadata.csv")
strata_df <- data.frame(stratalist)
head(strata_df)
strata(BW_dDocent_Dovetail_phylo_genlight) <- strata_df
BW_dDocent_Dovetail_phylo_genlight
head(strata(BW_dDocent_Dovetail_phylo_genlight))
popNames(BW_dDocent_Dovetail_phylo_genlight)
nameStrata(BW_dDocent_Dovetail_phylo_genlight) <- ~Individual/Collection/Year/State/Cluster/Host
BW_dDocent_Dovetail_phylo_genlight

BW_dDocent_Dovetail_phylo_AMOVA <- poppr.amova(BW_dDocent_Dovetail_phylo_genlight, hier = ~Year/Collection, clonecorrect = FALSE, within = TRUE,
                                               dist = NULL, squared = TRUE, freq = TRUE,
                                               correction = "quasieuclid", sep = "_", filter = FALSE,
                                               threshold = 0, algorithm = "farthest_neighbor", threads = 1L,
                                               missing = "loci", cutoff = 0.05, quiet = FALSE,
                                               method = c("ade4", "pegas"), nperm = 0)
BW_dDocent_Dovetail_phylo_AMOVA
amova.test <- randtest(BW_dDocent_Dovetail_phylo_AMOVA)
amova.test
plot(amova.test)

## basic info
BW_dDocent_Dovetail_phylo_genepop_basic_info <- basic_info("./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop.txt", outputFile = "./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop_basic_info.txt")

## gene diversity and Fis
BW_dDocent_Dovetail_phylo_genepop_gendivFis <- genedivFis("./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop.txt", sizes = FALSE, outputFile = "./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop_gene_diversity.txt", dataType = "Diploid")

## run a global Fst analysis using genepop (create input file using PGDSpider)
BW_dDocent_Dovetail_phylo_genepop_global_Fst <- Fst("./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop.txt", sizes = FALSE, pairs = FALSE, outputFile = "./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop_global_Fst.txt", dataType = "Diploid")

## run a pairwise Fst analysis using genepop (create input file using PGDSpider)
BW_dDocent_Dovetail_phylo_genepop_Fst <- Fst("./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop.txt", sizes = FALSE, pairs = TRUE, outputFile = "./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop_pairwise_Fst.txt", dataType = "Diploid")

## run a pairwise differentiation test  using genepop (create input file using PGDSpider)
BW_dDocent_Dovetail_phylo_genepop_test_diff <- test_diff("./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop.txt", genic = FALSE, pairs = TRUE, outputFile = "./BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_biallelic_genepop_test_diff.txt", dememorization = 1000, batches = 10, iterations = 500)