#set working directory
setwd("C:/Users/rs14/Documents/research/breeding_optimization/PyBrOpt/experiments/2020-Jan-09_exp01/")

# load data frame
pts_df <- read.csv("paf_topology_struct_unequalafreq.tsv", sep='\t', header=T)

model1 <- lm(paa ~ paf, data=pts_df)
with(pts_df, plot(paa, paf, xlab="PAA", ylab="PAF", main="PAA vs PAF", pch=3))
abline(model1, col="red", lwd=2)

model2 <- lm(paa ~ pafd, data=pts_df)
with(pts_df, plot(paa, pafd, xlab="PAA", ylab="PAFD", main="PAA vs PAFD", pch=3))
abline(model2, col="red", lwd=2)
