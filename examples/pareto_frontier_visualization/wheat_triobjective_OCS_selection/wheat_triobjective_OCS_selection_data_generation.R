#!/usr/bin/env Rscript

# load libraries
library(BGLR)
library(rrBLUP)

# load wheat dataset from BGLR
data(wheat)

# get marker matrix in {0,1,2} encoding
Z = 2 * wheat.X

# fit rrBLUP models for yield data in location 1, 2, 3, 4
fm1 = mixed.solve(y = wheat.Y[,1], Z = Z)
fm2 = mixed.solve(y = wheat.Y[,2], Z = Z)
fm3 = mixed.solve(y = wheat.Y[,3], Z = Z)
fm4 = mixed.solve(y = wheat.Y[,4], Z = Z)

# extract fixed effects (intercepts) for yield models in locations 1, 2, 3, 4
beta_df = data.frame(
    "fixed_effect" = c("intercept"),
    "yield1" = fm1$beta,
    "yield2" = fm2$beta,
    "yield3" = fm3$beta,
    "yield4" = fm4$beta
)

# extract additive marker effects for yield models in locations 1, 2, 3, 4
u_a_df = data.frame(
    "marker_name" = colnames(Z),
    "yield1" = fm1$u,
    "yield2" = fm2$u,
    "yield3" = fm3$u,
    "yield4" = fm4$u
)

# save the marker matrix to csv file in {0,1,2} format
write.csv(Z, "wheat_markers.csv")

# save the fixed effects to csv file
write.csv(beta_df, "wheat_intercepts.csv", row.names = FALSE)

# save the marker effects to csv file
write.csv(u_a_df, "wheat_marker_effects.csv", row.names = FALSE)
