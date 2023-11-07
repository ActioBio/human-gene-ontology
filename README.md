# Human Gene Ontology annotations

The parent project provides easy-to-use [Gene Ontology](http://geneontology.org/) annotations.

This repository hosts a modified version of the gene ontology annotation script, which has been adapted to focus exclusively on human gene annotations by filtering for tax_id 9606—Homo sapiens.

### Execution

```
conda env create -f environment.yml 

conda activate gene-ontology

bash download.sh

python process.py
```

### Input

- go-basic.obo
  - The Gene Ontology file in OBO format. It contains ontological information in a structured form, describing gene products in terms of their associated biological processes, cellular components, and molecular functions in a species-independent manner.
- gene2go.gz
  - The file from NCBI that links genes from the Entrez Gene database to Gene Ontology (GO) terms, establishing a connection between gene identifiers and their functional annotations.
- gene_info.gz
  - The file from NCBI is a compressed archive containing detailed information on genes.

### Output
- GO_annotations-9606-inferred-allev.tsv
  - inferred + direct annotations with all evidence
- GO_annotations-9606-direct-allev.tsv
  - direct annotations with all evidence
- GO_annotations-9606-inferred-expev.tsv
  - inferred + direct annotations with experimental evidence
- GO_annotations-9606-direct-expev.tsv
  - direct annotations with experimental evidence
