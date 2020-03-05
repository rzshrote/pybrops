library(Rmisc)
library(ggplot2)
library(ggpubr)

setwd("C:/Users/rs14/Documents/research/breeding_optimization/PyBrOpt/experiments/2020-Feb-29_exp_hc_vs_icpso/")

n_trials <- 120
n_phases <- 2
n_sel <- 10
n_loci <- 10000
n_cycle <- 5

compare_allele_props <- function(icpso_vec, hc_vec, n1, n2, sign) {
    pvals <- rep(NA, length(icpso_vec))
    for(i in 1:length(icpso_vec)) {
        test_result <- prop.test(
            c(icpso_vec[i], hc_vec[i]),
            c(n1[i], n2[i]),
            #alternative = "two.sided"
            alternative = if(sign[i] < 0) "less" else "greater"
        )
        pvals[i] <- test_result$p.value
    }
    return(pvals)
}

# this function does not work because R is not logical.
fdr_correction <- function(pvals, cycles, loci) {
    fdr <- c()
    for(i in 1:cycles) {
        adj <- p.adjust(c(pvals[(loci*(i-1))+1:loci*i]), method="fdr")
        fdr <- append(fdr, adj)
    }
    return(fdr)
}

allele_df <- read.csv("sel_allele_counts.tsv", sep="\t", header=T)

hc_sa_set_df <- allele_df[allele_df$algorithm == "hc_sa_set",]
icpso_df <- allele_df[allele_df$algorithm == "icpso",]

n <- rep(n_sel*n_phases*n_trials, n_loci*n_cycle)

pvals <- compare_allele_props(
    icpso_df$marker_count,
    hc_sa_set_df$marker_count,
    n,
    n,
    icpso_df$marker_effect
)

fdr1 <- p.adjust(pvals[1:10000], method="fdr")
fdr2 <- p.adjust(pvals[10001:20000], method="fdr")
fdr3 <- p.adjust(pvals[20001:30000], method="fdr")
fdr4 <- p.adjust(pvals[30001:40000], method="fdr")
fdr5 <- p.adjust(pvals[40001:50000], method="fdr")
fdr_vals <- c(fdr1, fdr2, fdr3, fdr4, fdr5)

# make summarizing df
summary_df <- data.frame(icpso_df)
summary_df <- summary_df[,-c(1,7)]
summary_df[,"pval"] <- pvals
summary_df[,"fdr"] <- fdr_vals


cycle5_df <- summary_df[summary_df$cycle == 5,]
# plot(density(cycle5_df[cycle5_df$fdr < 0.05,]$marker_effect))
# lines(density(cycle5_df$marker_effect))
# abline(v=0)


# do not use this
# density_plot <- ggplot(cycle5_df, aes(x=marker_effect)) +
#     labs(
#         title = "Marker Effects of Loci Under Differential Selection at Cycle 5",
#         subtitle = "Proportion differences between ICPSO and HC algorithm strategies"
#     ) +
#     theme(
#         plot.title = element_text(hjust = 0.5),
#         plot.subtitle = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA)
#     ) +
#     xlab("Marker Effect") +
#     ylab("Probability Density") +
#     geom_vline(xintercept = 0) +
#     geom_hline(yintercept = 0) +
#     geom_density(
#         aes(color = "All markers"),
#         fill="blue",
#         alpha=0.4,
#         show.legend = TRUE
#     ) +
#     geom_density(
#         data=cycle5_df[cycle5_df$fdr < 0.05,],
#         aes(color = "Selected markers"),
#         fill="red",
#         alpha=0.4,
#         show.legend = TRUE
#     ) +
#     scale_color_manual(
#         name = "Distribution",
#         values = c("black", "black")
#     ) +
#     scale_fill_manual(
#         values = c("red","blue")
#     )

# make df with significant and all.
all_vec <- cycle5_df$marker_effect
sig_vec <- cycle5_df[cycle5_df$fdr < 0.05,]$marker_effect

density_df <- data.frame(
    Distribution = c(rep("All markers", length(all_vec)), rep("Selected markers", length(sig_vec))),
    effect = c(all_vec, sig_vec)
)

density_plot <- ggplot() +
    labs(
        title = "Marker Effects of Loci Under Differential Selection at Breeding Cycle 5",
        subtitle = "From allele frequency differences between ICPSO and HC algorithm strategies"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(0.85,0.9)
    ) +
    xlab("Marker Effect") +
    ylab("Probability Density") +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_density(
        data = density_df,
        aes(x = effect, color = Distribution, fill = Distribution),
        alpha=0.5
    )
ggsave(width=6.5, height=6, plot=density_plot, filename="marker_density.png")


library(qqman)


cycle5_df[,"pos_int"] <- as.integer(cycle5_df[,"gmap_pos"] * 1000000)

sig_pval = -log10(max(cycle5_df[cycle5_df$fdr < 0.05,]$pval))


png("manhattan.png", units="in", width=6, height=6, res=300)
manhattan(
    cycle5_df,
    main = "Loci Under Differential Selection at Breeding Cycle 5",
    chr="chr",
    bp="pos_int",
    p="pval",
    snp="marker_id",
    suggestiveline = FALSE,
    genomewideline = sig_pval,
    col = c("#00b050", "#00b0f0", "#0070c0", "#002060", "#7030a0")
)
dev.off()

col = c("#000000", "#ff0000", "#ffc000", "#ffff00", "#92d050", "#00b050", "#00b0f0", "#0070c0", "#002060", "#7030a0")

    labs(
        title = "Loci Under Differential Selection at Breeding Cycle 5",
        subtitle = "From allele frequency differences between ICPSO and HC algorithm strategies"
    ) +
