# The purpose of this experiment is to determine if PA outperforms OPV
# AND assess its effects on genetic gain and genetic diversity

# load various libraries we need
library(Rmisc)
library(ggplot2)
library(ggpubr)

# set working directory
setwd("C:/Users/rs14/Documents/research/breeding_optimization/PyBrOpt/experiments/2020-Jan-05/")

# read file
pa_opv_df <- read.csv("pa_opv_pop.tsv", sep='\t', header=T)

# grab gebv header names
gebv_cols <- names(pa_opv_df)[7:56]

# convert ints to factors
pa_opv_df[,"cycle"] <- as.factor(pa_opv_df[,"cycle"])

# calculate row means
pa_opv_df[,"gebv_mean"] <- rowMeans(pa_opv_df[,gebv_cols])

# make plot of gain in GEBV over time
ggline(pa_opv_df, x="cycle", y="gebv_mean", add="mean_se", color="method") +
    ggtitle("Effect of selection on increase in GEBV over time") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Breeding cycle") +
    ylab("Mean population GEBV") +
    stat_compare_means(aes(group = method), method="t.test", label="p.signif",
        label.y=c(7,15,22,28,33,38,43,47,51,54))
ggsave(width=8, height=6, filename="mean_gebv.png")
dev.off()

# make plot of genetic diversity over time
ggline(pa_opv_df, x="cycle", y="opv_score", add="mean_se", color="method") +
    ggtitle("Effect of selection on genetic diversity over time") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Breeding cycle") +
    ylab("Maximum haploid GEBV") +
    stat_compare_means(aes(group = method), method="t.test", label="p.signif",
        label.y=c(53,47,44,43,42,41,40,39,38,37))
ggsave(width=8, height=6, filename="genetic_diversity.png")
dev.off()
