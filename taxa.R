
library(microbiomeSeq)
p <- plot_taxa(pseq, grouping_column = "SampleType", method = "hellinger", number.taxa = 21, 
               filename = NULL)
print(p)