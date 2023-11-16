# Download the Gene Ontology
wget --timestamping --directory-prefix data/input/ http://purl.obolibrary.org/obo/go/go-basic.obo

# Download the gene2go file producted by Entrez Gene
wget --timestamping --directory-prefix data/input/ ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

# Download the Human gene_info.gz file of Entrez Genes
wget --timestamping --directory-prefix data/input/ ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz