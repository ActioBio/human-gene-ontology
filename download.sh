# Download the Gene Ontology
wget --timestamping --directory-prefix data/input/ http://purl.obolibrary.org/obo/go/go-basic.obo

# Download the gene2go file producted by Entrez Gene
wget --timestamping --directory-prefix data/input/ ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

# Download the gene_info.gz file of Entrez Genes
wget --timestamping --directory-prefix data/input/ ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
