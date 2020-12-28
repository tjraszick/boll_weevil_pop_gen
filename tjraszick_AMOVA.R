## load the vcfR package and adegenet packages
library("vcfR")
library("adegenet")
library("genepop")
library("poppr")
library("ggplot2")

## import the vcf file of interest
vcf_file <- "BW_dDocent_Dovetail_phylo_lowDP_0.5_phred30_final.recode.vcf"
vcf <- read.vcfR(vcf_file, verbose=FALSE)

## convert vcfR to genlight object
BWpopgen_genlight <- vcfR2genlight(vcf, n.cores = 12)

popmap <- read.table("BW_dDocent_Dovetail_popmap.txt", sep="\t")
pop_vector <- as.vector(popmap[,2])

## append to the genlight object
pop(BWpopgen_genlight) <- as.factor(c(pop_vector))

## Load the file of interest and check the basic statistics about the file
BWpopgen_genlight
nLoc(BWpopgen_genlight)
pop(BWpopgen_genlight)
table(pop(BWpopgen_genlight))
indNames(BWpopgen_genlight)

## run an analysis of molecular variance (AMOVA) create input csv file on your own
stratalist <- read.csv("./BW_dDocent_Dovetail_phylo_metadata.csv")
strata_df <- data.frame(stratalist)
head(strata_df)
strata(BWpopgen_genlight) <- strata_df
BWpopgen_genlight
head(strata(BWpopgen_genlight))
popNames(BWpopgen_genlight)
nameStrata(BWpopgen_genlight) <- ~Individual/Collection/Year/State/Cluster/Host
BWpopgen_genlight

BW_dDocent_Dovetail_phylo_AMOVA <- poppr.amova(BWpopgen_genlight, hier = ~Year/Collection, clonecorrect = FALSE, within = TRUE,
                                               dist = NULL, squared = TRUE, freq = TRUE,
                                               correction = "quasieuclid", sep = "_", filter = FALSE,
                                               threshold = 0, algorithm = "farthest_neighbor", threads = 1L,
                                               missing = "loci", cutoff = 0.05, quiet = FALSE,
                                               method = c("ade4", "pegas"), nperm = 10000)
BW_dDocent_Dovetail_phylo_AMOVA
amova.test <- randtest(BW_dDocent_Dovetail_phylo_AMOVA, nrepet=9999)
amova.test
plot(amova.test)
