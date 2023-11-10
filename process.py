import os
import pandas as pd
import networkx as nx
import obonet
import csv
import gzip

HUMAN_TAX_ID = 9606

INPUT_FILE = {
    "GENE2GO": "data/input/gene2go.gz",
    "GENE_INFO": "data/input/gene_info.gz",
    "GO_OBO": "data/input/go-basic.obo"
}

INPUT_COLUMN = {
    "GENE_INFO": ['tax_id', 'GeneID', 'Symbol'],
    "GENE2GO": ['tax_id', 'GeneID', 'GO_ID', 'Evidence', 'Qualifier']
}

INPUT_DTYPE = {
    'tax_id': 'Int64', 'GeneID': 'Int64', 'Symbol': str, 'Qualifier': str
}

OUTPUT_DIR = "data/output"

REMOVE_SUBSETS = {'goantislim_grouping', 'gocheck_do_not_annotate', 'gocheck_do_not_manually_annotate'}
EXPERIMENTAL_CODES = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'}
PROPAGATE_ALONG = {'is_a', 'part_of'}
GO_DOMAINS = {'biological_process', 'molecular_function', 'cellular_component'}

def join_list_with_pipe(lst):
    return '|'.join(map(str, lst))

def load_filtered_dataframe_iteratively(file_key):
    path = INPUT_FILE[file_key]
    selected_data = []

    with gzip.open(path, 'rt') as f:
        header_line = ''
        while not header_line.startswith('#'):
            header_line = f.readline().strip()
        adjusted_header = header_line.replace('#', '').split('\t')
        col_indices = [adjusted_header.index(col) for col in INPUT_COLUMN[file_key]]
        tax_id_index = adjusted_header.index('tax_id') 

        for line in f:
            split_line = line.strip().split('\t')
            if split_line[tax_id_index] == str(HUMAN_TAX_ID):
                selected_data.append([split_line[i] for i in col_indices])

    return pd.DataFrame(selected_data, columns=INPUT_COLUMN[file_key])

def read_go_to_graph():
    with open(INPUT_FILE["GO_OBO"]) as read_file:
        graph = obonet.read_obo(read_file)
        graph.remove_nodes_from({
            node for node, data in graph.nodes(data=True)
            if REMOVE_SUBSETS & set(data.get('subset', []))
        })
        graph.remove_edges_from([
            (u, v, key) for u, v, key in graph.edges(data=False, keys=True)
            if key not in PROPAGATE_ALONG
        ])
        assert nx.is_directed_acyclic_graph(graph)
        return graph

def go_graph_to_dataframe(graph):
    rows = [(node, data['name'], data['namespace']) for node, data in graph.nodes(data=True)]
    go_df = pd.DataFrame(rows, columns=['go_id', 'go_name', 'go_domain'])
    return go_df.sort_values('go_id')

def setup_annotation_keys(graph):
    keys = ('direct_annotations', 'direct_not_annotations', 'inferred_annotations')
    for node in graph:
        graph.nodes[node].update({key: set() for key in keys})

def process_annotations(graph, goa_df):
    not_qualifier = goa_df['Qualifier'].apply(lambda qualifier: qualifier.startswith('NOT'))
    # Filtering direct annotations (where Qualifier is not 'NOT')
    direct_annotations = goa_df.loc[~not_qualifier, ['GO_ID', 'GeneID']]
    direct_not_annotations = goa_df.loc[not_qualifier, ['GO_ID', 'GeneID']]

    for go_id, gene in direct_annotations.itertuples(index=False):
        if go_id in graph:
            graph.nodes[go_id]['direct_annotations'].add(gene)

    for go_id, gene in direct_not_annotations.itertuples(index=False):
        if go_id in graph:
            graph.nodes[go_id]['direct_not_annotations'].add(gene)

def propagate_direct_annotations(graph):
    for node in nx.topological_sort(graph):
        data = graph.nodes[node]
        data['inferred_annotations'].difference_update(data['direct_not_annotations'])
        data['inferred_annotations'].update(data['direct_annotations'])
        # Propagate to successor nodes
        for successor_node in graph.successors(node):
            graph.nodes[successor_node]['inferred_annotations'].update(data['inferred_annotations'])

def annotate_and_propagate(graph, goa_df):
    graph = graph.copy()
    setup_annotation_keys(graph)
    process_annotations(graph, goa_df)
    propagate_direct_annotations(graph)
    return graph

def create_annotations_df(graph_annot, annotation_key, exclude_keys=None):
    exclude_keys = exclude_keys or []
    annotations = [
        (node, gene) for node, genes in graph_annot.nodes(data=annotation_key)
        for gene in genes if not any(gene in graph_annot.nodes[node][key] for key in exclude_keys)
    ]
    return pd.DataFrame(annotations, columns=['go_id', 'GeneID'])

def merge_gene_info(annotations_df, gene_df):
    return annotations_df.merge(gene_df, on='GeneID')

def aggregate_gene_ids_and_symbols(annotations_df):
    return annotations_df.groupby('go_id').agg({
        'GeneID': lambda x: join_list_with_pipe(sorted(set(x))),
        'Symbol': lambda x: join_list_with_pipe(sorted(set(x)))
    }).reset_index()

def extract_annotation_df(graph_annot, gene_df, go_df):
    direct_df = create_annotations_df(graph_annot, 'direct_annotations')
    inferred_df = create_annotations_df(graph_annot, 'inferred_annotations', ['direct_annotations'])

    direct_df = merge_gene_info(direct_df, gene_df)
    inferred_df = merge_gene_info(inferred_df, gene_df)

    direct_agg = aggregate_gene_ids_and_symbols(direct_df)
    inferred_agg = aggregate_gene_ids_and_symbols(inferred_df)

    final_df = direct_agg.merge(inferred_agg, on='go_id', how='outer', suffixes=('_direct', '_inferred'))
    final_df = go_df.merge(final_df, on='go_id', how='right')
    final_df = final_df[['go_id', 'go_name', 'go_domain', 'GeneID_direct', 'GeneID_inferred', 'Symbol_direct', 'Symbol_inferred']]

    return final_df

def create_node_csv_files(annotation_df):
    for domain in GO_DOMAINS:
        domain_df = annotation_df[annotation_df['go_domain'] == domain][['go_id', 'go_name']]
        node_file_name = f'node_{domain}.csv.gz'
        path = os.path.join(OUTPUT_DIR, node_file_name)
        domain_df.to_csv(path, index=False, compression='gzip') 

def create_edge_csv_files(annotation_df, gene2go_df):
    gene2go_evidence_df = gene2go_df[gene2go_df['Evidence'].isin(EXPERIMENTAL_CODES)]
    gene_go_identifier_set = set(gene2go_evidence_df.apply(lambda row: f"{row['GeneID']}_{row['GO_ID']}", axis=1))

    files = {domain: gzip.open(os.path.join(OUTPUT_DIR, f'edge_gene_to_{domain}.csv.gz'), 'wt', newline='') for domain in GO_DOMAINS} 
    writers = {domain: csv.writer(file) for domain, file in files.items()}

    for writer in writers.values():
        writer.writerow(['gene_id', 'go_id', 'type', 'experimental'])

    for _, row in annotation_df.iterrows():
        go_id = row['go_id']
        direct_gene_ids = extract_gene_ids(row['GeneID_direct'])
        inferred_gene_ids = extract_gene_ids(row['GeneID_inferred'])
        domain = row['go_domain']

        if direct_gene_ids:
            is_experimental = [str(gene_id) + '_' + go_id in gene_go_identifier_set for gene_id in direct_gene_ids]
            for i, gene_id in enumerate(direct_gene_ids):
                writers[domain].writerow([gene_id, go_id, 'direct', is_experimental[i]])

        if inferred_gene_ids:
            for gene_id in inferred_gene_ids:
                writers[domain].writerow([gene_id, go_id, 'inferred', False])

    for file in files.values():
        file.close()

def extract_gene_ids(gene_ids_str):
    return [int(gene_id) for gene_id in gene_ids_str.split('|')] if not pd.isna(gene_ids_str) and gene_ids_str else []

def main():
    gene_df = load_filtered_dataframe_iteratively("GENE_INFO")
    gene2go_df = load_filtered_dataframe_iteratively("GENE2GO")
    go_graph = read_go_to_graph()
    go_df = go_graph_to_dataframe(go_graph)
    graph_annot = annotate_and_propagate(go_graph, gene2go_df)
    annotation_df = extract_annotation_df(graph_annot, gene_df, go_df)
    create_node_csv_files(annotation_df)
    create_edge_csv_files(annotation_df, gene2go_df)

if __name__ == "__main__":
    main()