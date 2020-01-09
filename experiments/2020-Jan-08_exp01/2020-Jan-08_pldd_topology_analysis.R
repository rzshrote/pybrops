#set working directory
setwd("C:/Users/rs14/Documents/research/breeding_optimization/PyBrOpt/experiments/2020-Jan-08_exp01/")

# load data frame
pts_df <- read.csv("pldd_topology.tsv", sep='\t', header=T)

model1 <- lm(paa ~ pldd, data=pts_df)

with(pts_df, plot(paa, pldd, xlab="PAA", ylab="PLDD", main="PAA vs PLDD", pch=3))
abline(model1, col="red", lwd=2)
