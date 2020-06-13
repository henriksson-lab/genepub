wget https://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz
wget https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz

gunzip gene2ensembl.gz
gunzip Mus_musculus.gene_info.gz

grep ENSMUSG gene2ensembl > gene2ensembl.mouse
