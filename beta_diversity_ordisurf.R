
library(microbiomeSeq)

# reference: http://userweb.eng.gla.ac.uk/umer.ijaz/projects/microbiomeSeq_Tutorial.html

# Beta-Diversity

bray_dist <- phyloseq::distance(ps, "bray")
nmds_bray <- ordinate(ps, method = "NMDS", distance = bray_dist, grouping_column = "genotype", 
                      pvalue.cutoff = 0.05)
plot_ordination(ps, nmds_bray, color = "genotype")


# Ordisurf
# ord.res for nmds_bray

p <- plot_ordisurf(ord.res, meta_table, variable = "pH")
print(p)


# fuzy


p <- generateFSO(ps, grouping_column = "Souce", method = 2, indices = 2, 
                 filename = NULL)
print(p)


# Canonical Correspondence Analysis

plot_cca(physeq = physeq, grouping_column = "Country", pvalueCutoff = 0.01, 
         env.variables = NULL, num.env.variables = NULL, exclude.variables = "Country", 
         draw_species = F)




# other nice ordination plot with the physico-chemical parameters
# http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#data_import

