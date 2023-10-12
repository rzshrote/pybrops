#!/usr/bin/env Rscript

mocgsrs_df <- read.csv("dh_mocgsrs_program.csv", header = TRUE)
mocgsrs_df$t_cur = as.factor(mocgsrs_df$t_cur)

library(ggplot2)
library(ggpubr)
library(gridExtra)

mocgsrs_plot = ggplot(
    data = mocgsrs_df,
    aes(x = main_true_mean_syn1, y = main_true_mean_syn2)
) + geom_line(
) + geom_point(
    aes(color = t_cur)
) + ggtitle(
    "Multi-Objective Conventional Genomic Selection Over Time",
    subtitle = "Population Mean Trait Values"
) + xlab(
    "Synthetic Trait 1 (h² = 0.4)"
) + ylab(
    "Synthetic Trait 2 (h² = 0.6)"
) + theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
) + scale_colour_discrete(
    "Generation"
)

ggsave(
    "breeding_value_plot.png",
    mocgsrs_plot,
    width = 6,
    height = 6,
    units = "in",
    dpi = 300
)

mocgsrs_plot <- ggline(
    data = mocgsrs_df,
    x = "t_cur",
    y = c("main_mehe"),
    add = "mean_sd",
    color = "blue",
    merge = TRUE,
    title = "Multi-Objective Conventional Genomic Selection Over Time",
    subtitle = "Mean Expected Heterozygosity",
    xlab = "Generation",
    ylab = "Mean expected heterozygosity"
) + theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
)

ggsave(
    "mean_expected_heterozygosity_plot.png",
    mocgsrs_plot,
    width = 6,
    height = 6,
    units = "in",
    dpi = 300
)

# mocgsrs_plot = ggplot(
#     data = mocgsrs_df,
#     aes(x = main_true_var_A_syn1, y = main_true_var_A_syn2)
# ) + geom_line(
# ) + geom_point(
#     aes(color = t_cur)
# ) + ggtitle(
#     "Multi-Objective Conventional Genomic Selection Over Time",
#     subtitle = "Population Additive Genetic Variance"
# ) + xlab(
#     "Synthetic Trait 1 (h² = 0.4)"
# ) + ylab(
#     "Synthetic Trait 2 (h² = 0.6)"
# ) + theme(
#     plot.title = element_text(hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5)
# ) + scale_colour_discrete(
#     "Generation"
# )
#
# ggsave(
#     "genic_genetic_variance_plot.png",
#     mocgsrs_plot,
#     width = 6,
#     height = 6,
#     units = "in",
#     dpi = 300
# )
