#!/usr/bin/env Rscript
library(ggplot2)

frontier_df = read.csv("multiobjective_genomic_selection_program_frontier.csv", header = TRUE)
# convert time to factor
frontier_df$t_cur = as.factor(frontier_df$t_cur)
# negate the trait values since they are minimizing in nature
frontier_df$syn1 = -frontier_df$syn1
frontier_df$syn2 = -frontier_df$syn2

frontier_plot = ggplot(
    data = frontier_df,
    aes(x = syn1, y = syn2, color = t_cur)
) + geom_point(
) + ggtitle(
    "Pareto Frontier for Multi-Objective Conventional Genomic Selection Over Time"
) + xlab(
    "Synthetic Trait 1 (h² = 0.4)"
) + ylab(
    "Synthetic Trait 2 (h² = 0.6)"
) + theme(
    plot.title = element_text(hjust = 0.5)
) + scale_colour_discrete(
    "Generation"
)

ggsave(
    "multiobjective_genomic_selection_program_frontier.png",
    frontier_plot,
    width = 7,
    height = 7,
    units = "in",
    dpi = 300
)
