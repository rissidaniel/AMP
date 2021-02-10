# in R
write.table(t(seqtab.nochim), "dada2-analysis/seqtab-nochim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

# in R
uniquesToFasta(seqtab.nochim, fout='dada2-analysis/rep-seqs.fna', ids=colnames(seqtab.nochim))

# in q2
qiime tools import \
--input-path dada2-analysis/rep-seqs.fna \
--type 'FeatureData[Sequence]' \
--output-path rep-seqs.qza

# in q2
echo -n "#OTU Table" | cat - dada2-analysis/seqtab-nochim.txt > dada2-analysis/biom-table.txt

# in q2
biom convert -i dada2-analysis/biom-table.txt -o dada2-analysis/table.biom --table-type="OTU table" --to-hdf5

# in q2
qiime tools import \
--input-path dada2-analysis/table.biom \
--type 'FeatureTable[Frequency]' \
--source-format BIOMV210Format \
--output-path table.qza