## load the necessary packages
library(dplyr)
library(car)
library(multcomp)
library(multcompView)
library(ggpubr)

relatedness <- read.table("withinPopulationRelatedness", sep="\t", header = TRUE)
relatedness

group_by(relatedness, Population) %>%
  summarise(
    count = n(),
    mean = mean(Relatedness, na.rm = TRUE),
    sd = sd(Relatedness, na.rm = TRUE)
  )

## run a quick ANOVA
res.aov <- aov(relatedness$Relatedness ~ relatedness$Population, data = relatedness)
summary(res.aov)

## run a quick Tukey's HSD test
Tukey.matrix <- (TukeyHSD(res.aov))
Tukey.matrix

## Levene's test for homogeneity of variance
leveneTest(Relatedness ~ Population, data = relatedness)

## Welch's pairwise t-test with no assumption of equal variances
pairwise.t.test(relatedness$Relatedness, relatedness$Population,
                p.adjust.method = "bonferroni", pool.sd = FALSE)

## check the distribution of the residuals
plot(res.aov, 2)
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

## kruskal-Wallis rank sum test
kruskal.test(relatedness$Relatedness ~ relatedness$Population, data = relatedness)

## Wilcoxon rank test
wilcox.test <- pairwise.wilcox.test(jitter(relatedness$Relatedness), relatedness$Population, p.adjust.method = "bonferroni")
wilcox.semi.matrix <- as.data.frame(wilcox.test$p.value)
wilcox.semi.matrix
Export(wilcox.semi.matrix, file = "wilcox_semi_matrix.txt", format = "txt")

## use Excel to transpose semimatrix to full double matrix, then run the following code
wilcox.matrix <- (read.table("wilcox_matrix.txt", header = TRUE, row.names = 1, sep = "\t"))
wilcox.matrix

wilcox.letters <- multcompLetters(wilcox.matrix,
                compare="<",
                threshold=0.05,
                Letters=LETTERS,
                reversed = FALSE)

## boxplot
ggboxplot(relatedness, x = "Population", y = "Relatedness", 
          order = c("SMX-Caj (2014)", "CMX-Chi (2014)", "CMX-Dur (2014)", "RGV-Tam (2014)", "RGV-Tex (2014)", "AWC-Lem (2016)", "AWC-Sah (2016)", "AWC-H83 (2016)", "AWC-Cal (2016)", "AWC-Bi1 (2016)", "AWC-Bi2 (2016)", "RGV-Tam (2016)", "RGV-Tex (2016)", "SMX-Caj (2017)", "CMX-Coa (2017)", "ARG-Cha (2017)", "ARG-Sae (2017)", "ARG-Sal (2017)", "ARG-San (2017)", "ARG-For (2017)"),
          ylab = "Relatedness", xlab = "Population") + rotate_x_text(270) + stat_compare_means(method = "kruskal.test", label.x.npc = 0.35)

## add the letters generated in an image editor
wilcox.letters
