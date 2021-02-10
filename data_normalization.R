# Data normalisation
library(microbiomeSeq)

pseq <- normalise_data(pseq, norm.method = "relative", norm.meta = T)
pseq <- normalise_data(pseq, norm.method = "scale", type = "log")
pseq <- normalise_data(pseq, norm.method = "scale", type = "sqrt")


pseq <- normalise_data(pseq,method = "randomsubsample")
pseq <- normalise_data(pseq,method = "edgeRnorm")
pseq <- normalise_data(pseq,method = "proportion")