The data frame includes the gene coordinates, PRO-seq read counts, and summerized values of different features.
seqname: chromosome name 
start: the start position of the site on a specific gene
end: the end position of the site on a specific gene 
strand: indicating the gene is either on the + or - strand of DNA 
ensembl_gene_id: the id that marks each gene; one gene, one specific id
score: the read counts of PRO-seq experiment
All the columns after "score" are columns for different features. For example, "gc[,1]" means the normalized column of the feature "G+C content".

