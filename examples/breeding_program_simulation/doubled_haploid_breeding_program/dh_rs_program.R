#!/usr/bin/env Rscript

prs_df <- read.csv("dh_rs_program.csv", header = TRUE)

library(ggplot2)
library(ggpubr)
library(gridExtra)

prs_plot <- ggline(
    data = prs_df,
    x = "t_cur",
    y = c("cand_true_mean", "cand_true_usl", "cand_true_lsl"),
    add = "mean_sd",
    color = "blue",
    merge = TRUE,
    title = "Test of Phenotypic Recurrent Selection",
    xlab = "Breeding Cycle",
    ylab = "GEBV"
)

ggsave("breeding_value_plot.png", prs_plot, width = 5, height = 5, units = "in", dpi = 300)

prs_plot <- ggline(
    data = prs_df,
    x = "t_cur",
    y = c("cand_mehe"),
    add = "mean_sd",
    color = "blue",
    merge = TRUE,
    title = "Test of Phenotypic Recurrent Selection",
    xlab = "Breeding Cycle",
    ylab = "Mean expected heterozygosity"
)

ggsave("mean_expected_heterozygosity_plot.png", prs_plot, width = 5, height = 5, units = "in", dpi = 300)

prs_plot <- ggline(
    data = prs_df,
    x = "t_cur",
    y = c("cand_true_var_A", "cand_true_var_a"),
    add = "mean_sd",
    color = "blue",
    merge = TRUE,
    title = "Test of Phenotypic Recurrent Selection",
    xlab = "Breeding Cycle",
    ylab = "Genic and Genetic Variance"
)

ggsave("genic_genetic_variance_plot.png", prs_plot, width = 5, height = 5, units = "in", dpi = 300)
