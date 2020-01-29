# The purpose of this experiment is to determine if PAv2 outperforms OPV
# AND assess its effects on genetic gain and genetic diversity

# load various libraries we need
library(Rmisc)
library(ggplot2)
library(ggpubr)

# set working directory
setwd("C:/Users/rs14/Documents/research/breeding_optimization/PyBrOpt/experiments/2020-Jan-18_bug_testing/")

# read file
cgs_opv_pa_v2_df <- read.csv("cgs_opv_pa_v2.tsv", sep='\t', header=T)

# grab gebv header names
gebv_cols <- names(cgs_opv_pa_v2_df)[7:ncol(cgs_opv_pa_v2_df)]

# convert ints to factors
cgs_opv_pa_v2_df[,"bcycle"] <- as.factor(cgs_opv_pa_v2_df[,"bcycle"])

# calculate row means
cgs_opv_pa_v2_df[,"gebv_mean"] <- rowMeans(cgs_opv_pa_v2_df[,gebv_cols])

# make plot of gain in GEBV over time
ggline(cgs_opv_pa_v2_df, x="bcycle", y="gebv_mean", add="mean_se", color="method") +
    ggtitle("Effect of selection on increase in GEBV over time") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Breeding cycle") +
    ylab("Mean population GEBV") # +
    # stat_compare_means(aes(group = method), method="t.test", label="p.signif",
        # label.y=c(7,15,22,28,33,38,43,47,51,54))
ggsave(width=8, height=6, filename="mean_gebv.png")
dev.off()

# make plot of genetic diversity over time
ggline(cgs_opv_pa_v2_df, x="bcycle", y="opv_score", add="mean_se", color="method") +
    ggtitle("Effect of selection on valuable genetic diversity over time") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Breeding cycle") +
    ylab("Maximum GEBV (upper selection limit)") # +
    # stat_compare_means(aes(group = method), method="t.test", label="p.signif",
        # label.y=c(53,47,44,43,42,41,40,39,38,37))
ggsave(width=8, height=6, filename="genetic_diversity.png")
dev.off()
