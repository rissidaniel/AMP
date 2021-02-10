library(phyloseq)
#data(enterotype)
#plot_net(ps, type = "samples", color="Source", maxdist = 0.3)
#plot_net(enterotype, color="SeqTech", maxdist = 0.3, laymeth = "auto")
#plot_net(enterotype, color="SeqTech", maxdist = 0.3, laymeth = "svd")
#plot_net(enterotype, color="SeqTech", maxdist = 0.3, laymeth = "circle")
#plot_net(enterotype, color="SeqTech", shape="Enterotype", maxdist = 0.3, laymeth = "circle")


plot_net(ps, color="Source", distance="bray", )


# https://web.stanford.edu/class/bios221/Pune/Labs/Lab_networks/LabNetworks.html