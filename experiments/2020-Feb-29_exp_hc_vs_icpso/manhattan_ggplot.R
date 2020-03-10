library(ggplot2)

manhattan_ggplot <- function(df, chr = "chr", pos = "pos", pval = "pval",
   color = NULL, upbuf = 0.05, ylim = NULL, mai = NULL, omi = NULL,
   nrow = NULL, ncol = NULL, shape = 1, title = NULL, subtitle = NULL,
   background = ggplot2::element_blank(),
   chr_lab_pos = "bottom", chr_lab_placement = "outside", chr_legend = FALSE,
   chr_background = ggplot2::element_blank(),
   xlab = "Chromosome Position", ylab = expression(-log[10](p)),
   sig_line = NULL, xmargin = ggplot2::unit(0, "lines"), ymargin = NULL,
   xtext = ggplot2::element_blank()
) {
    # create internal data.frame
    data <- data.frame(
        chr = df[,chr],                 # get the chromosome row as a vector
        pos = df[,pos],                 # get positions as a vector
        pval = df[,pval],               # get p value vector
        nlog_pval = -log10(df[,pval])   # get the negative log10 of p values
    )

    n_uniq_chr <- length(unique(data$chr))  # get number of unique chromsomes
    nlog_pval_max <- max(data$nlog_pval)    # get max negative log p value

    # get graph limits
    if(is.null(ylim)) {
        ylim <- c(0, (nlog_pval_max * (1 + upbuf)))
    }

    # fill colors to have the number of chromosomes
    color <- rep_len(
        if(is.null(color)) "#000000" else color, # color values to repeat
        n_uniq_chr # number of times to repeat until
    )


    # plot chromosome positions on x-axis, -log(p) on y-axis
    outman <- ggplot2::ggplot(
            data = data,
            mapping = aes(pos, nlog_pval, color = factor(chr))
        ) +
        ggplot2::geom_point(     # makes the points display
            shape = shape        # point display type
        ) +
        ggplot2::geom_hline(
            yintercept = sig_line,
            color = "red"
        ) +
        ggplot2::facet_wrap(
            ~ chr,
            nrow = nrow,
            ncol = ncol,
            strip.position = chr_lab_pos
        ) +
        ggplot2::labs(
            title = title,
            subtitle = subtitle
        ) +
        ggplot2::theme(
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            panel.background = background,
            strip.placement = chr_lab_placement,
            strip.background = chr_background,
            panel.spacing.x = xmargin,
            panel.spacing.y = ymargin,
            axis.text.x=element_blank()
        ) +
        ggplot2::scale_color_manual(
            name = "Chromosome",
            values = color,
            guide = chr_legend
        ) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::scale_y_continuous(
            limits = ylim,
            expand = c(0,0)
        )


    return(outman)
}

cycle5_df <- read.csv("R/test.csv")
sig_line <- -log10(max(cycle5_df[cycle5_df$fdr < 0.05,"pval"]))

manplot <- manhattan_ggplot(cycle5_df, pos="gmap_pos", nrow = 1,
                 sig_line = sig_line,
                 title = element_text("Loci Under Differential Selection at Breeding Cycle 5"),
                 xlab = "Chromosome",
                 color = c("#00b050", "#00b0f0", "#0070c0", "#002060", "#7030a0"),
                 shape = 20,
                 ylab = expression(-log[10](p))
)
manplot
ggsave(width = 6.5, height = 6, plot = manplot, filename = "manhattan.png", dpi=400)
