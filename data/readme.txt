The data frame includes the gene coordinates, PRO-seq read counts, and summerized values of different features.
seqname: chromosome name 
start: the start position of the site on a specific gene
end: the start position of the site on a specific gene 
strand: the + or - strand that the gene is on 
ensembl_gene_id: the id for each gene; it is gene-specific
score: the read counts of PRO-seq data
All the columns after "score" are columns for different features. For example, "gc[,1]" means the normalized column of the feature "G+C content".

