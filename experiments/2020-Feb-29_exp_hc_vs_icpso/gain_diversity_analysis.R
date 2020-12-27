library(Rmisc)
library(ggplot2)
library(ggpubr)
library(gridExtra)


# set working directory
setwd("C:/Users/rs14/Documents/research/PyBrOpt/experiments/2020-Feb-29_exp_hc_vs_icpso")

# read file
pop_df <- read.csv("pop_summary_uniq.tsv", sep = "\t", header = T)
sel_df <- read.csv("sel_summary_uniq.tsv", sep = "\t", header = T)

# convert ints to factors
pop_df[,"bcycle"] <- as.factor(pop_df[,"cycle"])
sel_df[,"bcycle"] <- as.factor(sel_df[,"cycle"])
pop_df[,"trial"] <- as.factor(pop_df[,"trial"])
sel_df[,"trial"] <- as.factor(sel_df[,"trial"])

cycle0_df <- pop_df[pop_df$cycle == 0,]
cycle1_df <- pop_df[pop_df$cycle == 1,]
cycle2_df <- pop_df[pop_df$cycle == 2,]
cycle3_df <- pop_df[pop_df$cycle == 3,]
cycle4_df <- pop_df[pop_df$cycle == 4,]
cycle5_df <- pop_df[pop_df$cycle == 5,]

delta1 <- cycle1_df$mean_gebv - cycle0_df$mean_gebv
delta2 <- cycle2_df$mean_gebv - cycle1_df$mean_gebv
delta3 <- cycle3_df$mean_gebv - cycle2_df$mean_gebv
delta4 <- cycle4_df$mean_gebv - cycle3_df$mean_gebv
delta5 <- cycle5_df$mean_gebv - cycle4_df$mean_gebv

pop_gain_df <- data.frame(
    method = c(cycle1_df$method, cycle2_df$method, cycle3_df$method, cycle4_df$method, cycle5_df$method),
    algorithm = c(cycle1_df$algorithm, cycle2_df$algorithm, cycle3_df$algorithm, cycle4_df$algorithm, cycle5_df$algorithm),
    cycle = c(cycle1_df$cycle, cycle2_df$cycle, cycle3_df$cycle, cycle4_df$cycle, cycle5_df$cycle),
    delta = c(delta1, delta2, delta3, delta4, delta5)
)
# R converts factors into numbers for some dumb reason
pop_gain_df[pop_gain_df$method == 1,]$method <- "opv"
pop_gain_df[pop_gain_df$algorithm == 1,]$algorithm <- "hc_sa_set"
pop_gain_df[pop_gain_df$algorithm == 2,]$algorithm <- "icpso"
pop_gain_df[,"method"] <- as.factor(pop_gain_df[,"method"])
pop_gain_df[,"algorithm"] <- as.factor(pop_gain_df[,"algorithm"])

# make plot of GEBV over time
gebv_plot <- ggline(pop_df, x="cycle", y="mean_gebv", add="mean_se", color="algorithm") +
    ggtitle("Effect of Algorithm\non GEBV") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_discrete(
        name="Algorithm",
        labels=c("HC", "ICPSO")
    ) +
    xlab("Breeding Cycle") +
    ylab("Mean Population GEBV") +
    stat_compare_means(
        aes(group = algorithm),
        data = pop_df[pop_df$cycle > 0,],
        method="t.test",
        label="p.signif",
        label.y=c(49,53,57,61,65)
    )

# make plot of genetic gain over time
gain_plot <- ggline(pop_gain_df, x="cycle", y="delta", add="mean_se", color="algorithm") +
    ggtitle("Effect of Algorithm\non Genetic Gain") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_discrete(
        name="Algorithm",
        labels=c("HC", "ICPSO")
    ) +
    xlab("Breeding Cycle") +
    ylab("Mean Gain in GEBV") +
    stat_compare_means(
        aes(group = algorithm),
        method="t.test",
        label="p.signif",
        label.y=c(68,9,9,9,9)
    )

# make plot of genetic diversity over time
div_plot <- ggline(pop_df, x="cycle", y="score", add="mean_se", color="algorithm") +
    ggtitle("Effect of Algorithm on\nValuable Genetic Diversity") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_discrete(
        name="Algorithm",
        labels=c("HC", "ICPSO")
    ) +
    xlab("Breeding Cycle") +
    ylab("Maximum GEBV (Upper Selection Limit)") +
    stat_compare_means(
        aes(group = algorithm),
        data = pop_df[pop_df$cycle > 0,],
        method="t.test",
        label="p.signif",
        label.y=c(5800,5600,5400,5300,5200)
    )

# put plots together
whole_plot <- grid.arrange(gebv_plot, gain_plot, div_plot, nrow = 1)

ggsave(width=12, height=4, plot=whole_plot, filename="algorithm_effects_2.png")
dev.off()
