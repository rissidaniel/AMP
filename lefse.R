# lefse

# to use lefse wrapper the lefse needs to be installed in the pc and in the path
# https://github.com/ying14/yingtools2/wiki/Using-the-%60lefse%60-wrapper-function


library(devtools)
# install.packages("adespatial", dependencies = TRUE)
# install_github("umerijaz/microbiomeSeq", force = TRUE)
library(microbiomeSeq)
library(yingtools2)
library(phyloseq)
library(remotes)

# install.packages(c("mvtnorm", "modeltools", "coin"))
# remotes::install_github("ying14/yingtools2", dependencies = TRUE)

#remotes::install_github("umerijaz/microbiomeSeq", dependencies = TRUE)
#install.packages("lefse", dependencies = TRUE)

data("enterotype")

phy_obj1 <- (enterotype)


library(devtools)
library(ggplot2)


lefse.tbl <- lefse(physeq,class="Country", levels = "Genus")

data(pitlatrine)

physeq <- (pitlatrine)  

physeq <- normalise_data(physeq, norm.method = "relative", norm.meta = T)


lefs <- lefse(physeq,class="Country")


