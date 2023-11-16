# Human Gene Ontology annotations

The parent project provides easy-to-use [Gene Ontology](http://geneontology.org/) annotations.

This repository hosts a modified version of the gene ontology annotation script, which has been adapted to focus exclusively on human gene annotations by filtering for tax_id 9606—Homo sapiens.

### Execution

```
conda env create -f environment.yml 

conda activate human-gene-ontology

bash download.sh

python process.py
```

### Input

- go-basic.obo
  - The Gene Ontology file in OBO format. It contains ontological information in a structured form, describing gene products in terms of their associated biological processes, cellular components, and molecular functions in a species-independent manner.
- gene2go.gz
  - The file from NCBI that links genes from the Entrez Gene database to Gene Ontology (GO) terms, establishing a connection between gene identifiers and their functional annotations.
- Homo_sapiens.gene_info.gz
  - The file from NCBI is a compressed archive containing detailed information on genes.

### Output

- node_{domain}.csv.gz
  - Lists GO terms and names for specific domains (biological process, cellular component and molecular function).
- edge_gene_to_{domain}.csv.gz
  - Shows gene-GO term relationships for each domain, indicating if annotations are direct or inferred and if the evidence is experimental.

Note: These CSV files are formatted for easy import into Neo4j.