# The purpose of this experiment is to examine the relationship between
# allele availability and allele frequency (spoiler: they're very weakly
# correlated). This has applications for identifying an optimal point on a
# Pareto front.

#set working directory
setwd("C:/Users/rs14/Documents/research/breeding_optimization/PyBrOpt/experiments/2020-Jan-07/")

# load data frame
pts_df <- read.csv("paa_pad_1578411829.tsv", sep='\t', header=T)

# construct linear models
paa_pad_lm <- lm(pad ~ paa, data=pts_df)
paa_pad_prime_lm <- lm(pad_prime ~ paa, data=pts_df)
pad_pad_prime_lm <- lm(pad_prime ~ pad, data=pts_df)

# write model summaries to file
sink("paa_pad_1578411829_lm.txt")
summary(paa_pad_lm)
summary(paa_pad_prime_lm)
summary(pad_pad_prime_lm)
sink()

# extract coefficients
paa_pad_coeff <- round(coef(paa_pad_lm), 4)
paa_pad_prime_coeff <- round(coef(paa_pad_prime_lm), 4)
pad_pad_prime_coeff <- round(coef(pad_pad_prime_lm), 4)

# construct equations
paa_pad_eq <- paste0(
    "PAD = ",
    paa_pad_coeff[1],
    ifelse(sign(paa_pad_coeff[2]) == 1, " + ", " - "),
    abs(paa_pad_coeff[2]),
    " PAA"
)
paa_pad_prime_eq <- paste0(
    "PAD` = ",
    paa_pad_prime_coeff[1],
    ifelse(sign(paa_pad_prime_coeff[2]) == 1, " + ", " - "),
    abs(paa_pad_prime_coeff[2]),
    " PAA"
)
pad_pad_prime_eq <- paste0(
    "PAD` = ",
    pad_pad_prime_coeff[1],
    ifelse(sign(pad_pad_prime_coeff[2]) == 1, " + ", " - "),
    abs(pad_pad_prime_coeff[2]),
    " PAD"
)


# make plots
png("paa_pad_1578411829_lm1.png", width=700, height=700)
with(pts_df, plot(paa, pad, xlab="PAA", ylab="PAD", main="PAA vs PAD", pch=3))
abline(paa_pad_lm, col="red", lwd=2)
mtext(paa_pad_eq)
dev.off()
png("paa_pad_1578411829_lm2.png", width=700, height=700)
with(pts_df, plot(paa, pad_prime, xlab="PAA", ylab="PAD`", main="PAA vs PAD`", pch=3))
abline(paa_pad_prime_lm, col="red", lwd=2)
mtext(paa_pad_prime_eq)
dev.off()
png("paa_pad_1578411829_lm3.png", width=700, height=700)
with(pts_df, plot(pad, pad_prime, xlab="PAD", ylab="PAD`", main="PAD vs PAD`", pch=3))
abline(pad_pad_prime_lm, col="red", lwd=2)
mtext(pad_pad_prime_eq)
dev.off()
