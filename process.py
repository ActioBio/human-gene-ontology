import os
import pandas as pd
import networkx as nx
import obonet

OUTPUT_DIR = "data/output"
HUMAN_TAX_ID = 9606

FILE_PATHS = {
    "GENE2GO": "data/input/gene2go.gz",
    "GENE_INFO": "data/input/gene_info.gz",
    "GO_OBO": "data/input/go-basic.obo"
}
COLUMNS = {
    "GENE_INFO": ['tax_id', 'GeneID', 'Symbol', 'type_of_gene'],
    "GENE2GO": ['tax_id', 'GeneID', 'GO_ID', 'Evidence', 'Qualifier']
}
DTYPE = {
    'tax_id': 'Int64', 'GeneID': 'Int64', 'Symbol': str, 'type_of_gene': str, 'Qualifier': str
}
REMOVE_SUBSETS = {'goantislim_grouping', 'gocheck_do_not_annotate', 'gocheck_do_not_manually_annotate'}
EXPERIMENTAL_CODES = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'}
PROPAGATE_ALONG = {'is_a', 'part_of'}

def join_list_with_pipe(lst):
    return '|'.join(map(str, lst))

def load_filtered_dataframe(file_key):
    path = FILE_PATHS[file_key]
    return pd.read_csv(
        path,
        sep='\t',
        compression='gzip',
        comment='#',
        names=COLUMNS[file_key],
        usecols=COLUMNS[file_key],
        na_values=['-'],
        dtype=DTYPE
    ).query('tax_id == @HUMAN_TAX_ID')

def is_not_qualifier(qualifier):
    return not pd.isnull(qualifier) and qualifier.upper().startswith('NOT')

def read_go_to_graph():
    with open(FILE_PATHS["GO_OBO"]) as read_file:
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
    direct_annotations = goa_df.loc[~goa_df['Qualifier'].str.contains('NOT'), ['GO_ID', 'GeneID']]
    direct_not_annotations = goa_df.loc[goa_df['Qualifier'].str.contains('NOT'), ['GO_ID', 'GeneID']]

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

def main():
    gene_df = load_filtered_dataframe("GENE_INFO")
    gene2go_df = load_filtered_dataframe("GENE2GO")
    go_graph = read_go_to_graph()
    go_df = go_graph_to_dataframe(go_graph)

    for ev_type in ('allev', 'expev'):
        goa_subset_df = gene2go_df[gene2go_df['Evidence'].isin(EXPERIMENTAL_CODES)] if ev_type == 'expev' else gene2go_df
        graph_annot = annotate_and_propagate(go_graph, goa_subset_df)
        annotation_df = extract_annotation_df(graph_annot, gene_df, go_df)
        file_name = f'GO_annotations-{HUMAN_TAX_ID}-{ev_type}.tsv.gz'
        path = os.path.join(OUTPUT_DIR, file_name)
        annotation_df.to_csv(path, sep='\t', index=False, compression='gzip')

if __name__ == "__main__":
    main()