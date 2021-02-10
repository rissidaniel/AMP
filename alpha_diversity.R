library(microbiome)
library(knitr)
library(phyloseq)
library(kableExtra)
library(magrittr)
library(reshape2)
library(gridExtra)
library(ggplot2)


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("AnnotationDbi", "DESeq2", "GO.db", "impute", "phyloseq", "preprocessCore"))
#library(devtools)
#install_github("umerijaz/microbiomeSeq") 

library("microbiomeSeq")

# print all alpha diversity indicies
  values_diversity <-microbiome::alpha(ps, index = "all")
  write.csv(values_diversity, "diversity.csv")
  
  # Plot some alpha diverstiy
  
  # reference: http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/MicrobiomeWorkshop__MicrobiomeWorkshop/
  
  # replace "sex" for country
  
  # remove missing values (geom_errorbar)
  
  ps %<>%
    taxa_sums() %>%
    is_greater_than(0) %>%
    prune_taxa(ps)
  sample_data(ps)$Source = factor(sample_data(ps)$Source)
    # this plot the alpha diversity
  p = plot_richness(ps, x="Source", color="Source", measures=c("Observed", "InvSimpson", "Shannon", "Chao1")) + geom_jitter()
  p + geom_boxplot(data = p$data, aes(x = Source, y = value, color = NULL), alpha = 0.1)
  
  
   library("microbiomeSeq")
  # plot alpha diversity with statistical significancy using microbiomeSEQ (not installed)
  p<-plot_anova_diversity(ps, method = c("richness","simpson", "shannon"),grouping_column =  "Source",pValueCutoff=0.05)
  plot_anova_diversity(ps, method = c("richness","simpsom","shannon"), grouping_column = "Source")
  print(p)

