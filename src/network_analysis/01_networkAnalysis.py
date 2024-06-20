"""
@author: modified version of Gehad Youssef's code, by Anna Chechenina

Running time: depends on input, m * n = x combinations (start proteins: m, end proteins: n), 
    ca. 4 min for x = 1,800 on MacbookPro
    i.e. roughly x/500 = y minutes 
"""

import argparse
import networkx as nx
import pickle
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import numpy as np
import mygene

# Make sure the working folder is the src folder!

def parse_args():
    parser = argparse.ArgumentParser(usage='01_networkAnalysis.py [args]  ...',
                                     description='Perform protein network analysis')
    
    # INPUT ADDRESSES
    parser.add_argument('-s', '--start_address', default='../data/input_comp.txt', nargs='?',
                        help='a text file containing a start protein list in the Ensembl proein ID format, ENSP')
    parser.add_argument('-e', '--end_address', default='../data/input_comp.txt', nargs='?',
                        help='a text file containing an end protein list in the Ensembl proein ID format, ENSP')
    parser.add_argument('--edges', default='../databases/9606.protein.links.v11.0.noscores.threshold400.txt', nargs='?',
                        help='a list of edges from the STRING database')
    parser.add_argument('--prot_gene_list', default='../databases/list_protein_to_gene.txt', nargs='?',
                        help='a text file with a prote.in (ENSP) - gene (gene symbol) correspondence; default - based on Jensen Lab dictionary which is based on alias file from STRING database')
    parser.add_argument('--address_compartment_membrane_genes', default='../databases/COMPARTMENT_membrane_gene_names.txt', nargs='?',
                        help='a list of genes coding for membrane proteins')
    parser.add_argument('--address_compartment_nucleus_genes', default='../databases/COMPARTMENT_membrane_gene_names.txt', nargs='?',
                        help='a list of genes coding for nucleus proteins')
    
    # OUTPUT ADDRESSES
    parser.add_argument('--address_network_genes', default='../data/network_genes', nargs='?',
                        help='a network with protein IDs as nodes, .pkl will be added')
    parser.add_argument('--address_network_proteins', default='../data/network_proteins', nargs='?',
                        help='a network with gene symbols as nodes, .pkl will be added')
    parser.add_argument('--address_network_nodes_genes', default='../data/network_nodes_genes.txt', nargs='?',
                        help='a list of nodes as genes')
    parser.add_argument('--address_network_nodes_proteins', default='../data/network_nodes_proteins.txt', nargs='?',
                        help='a list of nodes as proteins')
    parser.add_argument('--address_network_edges_genes', default='../data/network_edges_genes.txt', nargs='?',
                        help='a list of edges as genes')
    parser.add_argument('--address_network_edges_proteins', default='../data/network_edges_proteins.txt', nargs='?',
                        help='a list of edges as proteins')
    parser.add_argument('--address_network_membrane_genes', default='../data/network_membrane_gene_names.txt', nargs='?',
                        help='a list of network genes in membrane')
    parser.add_argument('--address_network_nucleus_genes', default='../data/network_nucleus_gene_names.txt', nargs='?',
                        help='a list of network genes in nucleus')
    parser.add_argument('--address_centrality_results', default='../data/Centrality_RWR_results.txt', nargs='?',
                        help='a table of network parameters')
    parser.add_argument('--address_network_plot_genes', default='../visualisation/network_visualisation_genes.tiff', nargs='?',
                        help='a list of network genes in nucleus')
    parser.add_argument('--address_network_plot_proteins', default='../visualisation/network_visualisation_proteins.tiff', nargs='?',
                        help='a list of network genes in nucleus')

    return parser.parse_args()


# Function to load pickle object 
def load_obj(file_addr):
    with open(file_addr+ '.pkl', 'rb') as f:
        return pickle.load(f)


# Function to save pickle object
def save_obj(obj, file_addr ):
    with open(file_addr + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# Function to reformat iterable 
# s -> (s0,s1), (s1,s2), (s2, s3), ...
def pairwise(iterable):
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


def make_shortest_path_sources_to_targets(pair, string_G):
    """
    Function to find shortest paths on STRING between all start/end protein combinations
    Args:
        pair (tuple of str) - pair of start/end proteins
    Returns:
        edgelist of all shortest paths between start/end protein pair if existent
    ~ takes 4 min for an initial list of 43 proteins, i.e. 43x43 = ~1,850 combinations
    """

    # find shortest paths between start and end protein
    all_path_list = []
    try:
        # exclude self-loops
        if pair[0] != pair[1]:
            for p in nx.all_shortest_paths(string_G, source=pair[0], target=pair[1]):
                pairs_from_p = pairwise(p)
                all_path_list += pairs_from_p
                print(pair, 'HIT!')
    except:
        print(pair, "No Path")

    # return new edge list
    return all_path_list


def mapping_genes_id_to(prot_list,id_from, id_to, species):
    mg = mygene.MyGeneInfo()
    return mg.querymany(prot_list, scopes=id_from, fields= id_to, species=species, as_dataframe=True, verbose = True)


def calc_network_centrality_RWR(network, start_list, end_list, result_save_dir):
    """
    Function to calculate eigenvector centrality, degree centrality
    Args:
        network (networkx graph)
        start_list, end_list (list) - start and end lists to  use as source and target nodes
        result_save_dir (str) - path to save results
    Returns:
        csv file containing centrality and RWR-values for each node, txt-files containing above threshold genes
    """
    
    network_nodes = network.nodes()

    # eigenvector centrality
    try:
        eigenvector_centrality_dict = nx.eigenvector_centrality(network) #Returns:Dictionary of nodes with eigenvector centrality as the value.
    except:
        eigenvector_centrality_dict = dict()
        for node in network_nodes:
            eigenvector_centrality_dict[node] = 0.0
            
    # degree centrality
    try:
        degree_centrality_dict = nx.degree_centrality(network) #Returns:Dictionary of nodes with degree centrality as the value.
    except:
        degree_centrality_dict = dict()
        for node in network_nodes:
            degree_centrality_dict[node] = 0.0
            
    # betweeness centrality
    try:
        between_centrality_dict = nx.algorithms.betweenness_centrality(network) #Returns: Dictionary of nodes with betweenness centrality as the value.

    except:
        between_centrality_dict = dict()
        for node in network_nodes:
            between_centrality_dict[node] = 0.0
            
    # betweeness centrality subset: starts at membrane bound proteins ends in nuclear proteins
    try:
        betweensub_centrality_dict = nx.algorithms.betweenness_centrality_subset(network, sources=start_list, targets=end_list) #Returns: Dictionary of nodes with betweenness centrality as the value.

    except:
        betweensub_centrality_dict = dict()
        for node in network_nodes:
            betweensub_centrality_dict[node] = 0.0

    # edge betweenness centrality
    edge_BC_no_Zero = dict()
    edge_BC = nx.edge_betweenness_centrality(network) #?

    for k,v in edge_BC.items():
        edge_BC_no_Zero[k] = v+0.1 # add 0.1 to each edge score to make non-zero
    
    # set non-zero edge betweenness centrality as edge weight
    nx.set_edge_attributes(network, edge_BC_no_Zero, 'weight') #set edge weights as non-zero edge betweenness-centrality score

    # pagerank 
    # returns dictionary of nodes with PageRank as value
    PR_score = nx.pagerank(network)
    # personalised pagerank: jumps back to membrane proteins when jumping
    # in COVID-paper: start_genes_for_PR = start_list #{'COVID19':1}, here:
    dict_membrane = dict.fromkeys(start_list,1)
    PRsub_score = nx.pagerank(network, personalization=dict_membrane) 
    
    # export results as csv
    network_property_df = pd.DataFrame(
        columns=['Eigen', 'Degree', 'Between', 'Between Sub', 'RWR', 'RWR Sub'])
    for node in network_nodes:
        network_property_df.loc[node] = [eigenvector_centrality_dict[node], degree_centrality_dict[node],
                                         between_centrality_dict[node], betweensub_centrality_dict[node], 
                                         PR_score[node], PRsub_score[node]]

    network_property_df.to_csv(result_save_dir)


def make_circular_plot_network(nodes, G_gene, address_network_plot, start_list):
    """
    Function to make a circular plot of the network, save as tiff in visualisation folder
    Args:
        nodes (list) - list of nodes to use
        G_gene (networkx graph)
        address_network_plot (str) - path to save the plot
        start_list (list) - start nodes for highlighting
    Returns:
        none
    """
    
    color_map = []
    for node in nodes:
        if node in start_list:
            color_map.append('lightcoral')
        else:
            color_map.append('lightcyan')   
            
    plt.subplots(figsize=(30,30))
    
    n=len(nodes)
    angle = []
    angle_dict = {}
    for i, node in zip(range(n),nodes):
        theta = 2.0*np.pi*i/n
        angle.append((np.cos(theta),np.sin(theta)))
        angle_dict[node] = theta
    pos = {}
    for node_i, node in enumerate(nodes):
        pos[node] = angle[node_i]
        
    description = nx.draw_networkx_labels(G_gene, pos, font_size=20 )
    for node, t in description.items():
        t.set_rotation(angle_dict[node]*360.0/(2.0*np.pi))
                                          
    nx.draw_circular(G_gene, node_color=color_map, with_labels=False, font_size=18, font_color ='k', edge_color = 'grey')
    plt.savefig(address_network_plot)


#%% MAIN
# NETWORK RECONSTRUCTION

if __name__ == '__main__': 
    start_addr = parse_args().start_address
    end_addr = parse_args().end_address
    string_addr = parse_args().edges
    address_dictionary_protein_to_gene = parse_args().prot_gene_list
    address_compartment_membrane_genes = parse_args().address_compartment_membrane_genes
    address_compartment_nucleus_genes = parse_args().address_compartment_nucleus_genes
    address_network_genes = parse_args().address_network_genes
    address_network_proteins = parse_args().address_network_proteins
    address_network_nodes_genes = parse_args().address_network_nodes_genes
    address_network_nodes_proteins = parse_args().address_network_nodes_proteins
    address_network_edges_genes = parse_args().address_network_edges_genes
    address_network_edges_proteins = parse_args().address_network_edges_proteins
    address_network_membrane_genes = parse_args().address_network_membrane_genes
    address_network_nucleus_genes = parse_args().address_network_nucleus_genes
    address_centrality_results = parse_args().address_centrality_results
    address_network_plot_genes = parse_args().address_network_plot_genes
    address_network_plot_proteins = parse_args().address_network_plot_proteins

    print("Start program. Program will tell you when finished.")

    # open start / end protein list for network reconstruction
    start_df = pd.read_csv(start_addr, sep='\t', header = None)
    start_list_prot = start_df.iloc[:,0].to_list()
    end_df = pd.read_csv(end_addr, sep='\t', header = None)
    end_list_prot = end_df.iloc[:,0].to_list()

    print("Create network:")

    with open(string_addr) as network_f:
        string_network_edges = [x.strip().split(',') for x in network_f.readlines()]
    print('open STRING list done')

    # STRING edgelist to network
    string_G = nx.Graph(string_network_edges)
    print('make STRING network done')

    ppi_deg_pair_list = list(itertools.product(start_list_prot, end_list_prot))
    print('start/end protein combinations list done')

    shortest_paths_result =[]
    for i in range(len(ppi_deg_pair_list)):
        shortest_paths_result += make_shortest_path_sources_to_targets(ppi_deg_pair_list[i], string_G)

    # make new graph with hidden layers
    G_prot = nx.Graph()
    G_prot.add_edges_from(shortest_paths_result)

    # save new network with protein IDs as nodes
    save_obj(G_prot, address_network_proteins)
    #%% MAIN
    # TRANSLATION PROTEIN TO GENE

    print('\nRenaming network nodes (Protein ID to gene symbol).')

    # translate  protein IDs using local dictionary from JENSEN lab
    # open dictionary (ensembl protein IDs to gene symbols, JENSEN lab)
    network_dict = pd.read_csv(address_dictionary_protein_to_gene,  delimiter=',', header=None)
    network_dict = dict(sorted(network_dict.values.tolist()))
    #rename nodes in network
    G_gene = nx.relabel.relabel_nodes(G_prot, mapping=network_dict, copy=True)

    if len(G_gene.nodes) == len(G_prot.nodes):
        print('Renamed all nodes successfully.')
    else:
        print('Not all nodes renamed, or more than one protein corresponds to the same gene.')

    # also tranlsate the original input lists of proteins for later use
    start_list_gene = []
    for protein in start_list_prot:
        if protein in network_dict:
            start_list_gene += [network_dict[protein]]
        else:
            start_list_gene += [protein]
    end_list_gene = []
    for protein in end_list_prot:
        if protein in network_dict:
            end_list_gene += [network_dict[protein]]
        else:
            end_list_gene += [protein]

    # save list of nodes and network with gene names
    network_nodes_genes = pd.DataFrame(G_gene.nodes)
    network_nodes_genes.to_csv(address_network_nodes_genes, header=None, index=False)
    network_nodes_proteins = pd.DataFrame(G_prot.nodes)
    network_nodes_proteins.to_csv(address_network_nodes_proteins, header=None, index=False)
    network_edges_genes = pd.DataFrame(G_gene.edges)
    network_edges_genes.to_csv(address_network_edges_genes, header=None, index=False)
    network_edges_proteins = pd.DataFrame(G_prot.edges)
    network_edges_proteins.to_csv(address_network_edges_proteins, header=None, index=False)
    save_obj(G_gene, address_network_genes)

    #%% MAIN
    # PROTEIN LOCALISATION

    print('Finding proteins that are localised in nucleus/membrane.')
    # open pre-edited lists of genes in membrane and in nucleus (based on COMPARTMENT database on 30th July 2021)
    df_comp_mem = pd.read_csv(address_compartment_membrane_genes, delimiter='\t', header=None)
    df_comp_nuc = pd.read_csv(address_compartment_nucleus_genes, delimiter='\t', header=None)
    list_comp_mem = df_comp_mem.iloc[:,0].to_list()
    list_comp_nuc = df_comp_nuc.iloc[:,0].to_list()

    # convert network nodes to list
    nodes = list(G_gene.nodes)

    # make list of genes (=nodes) that are in membrane/nucleus
    network_membrane_gene_names = []
    network_nucleus_gene_names = []

    for node in nodes:
        if node in list_comp_mem:
            network_membrane_gene_names += [node]
        if node in list_comp_nuc:
            network_nucleus_gene_names += [node]

    # save lists as txt
    network_membrane_gene_names_list = pd.DataFrame(list(network_membrane_gene_names))
    network_nucleus_gene_names_list = pd.DataFrame(list(network_nucleus_gene_names))
    network_membrane_gene_names_list.to_csv(address_network_membrane_genes, header=None, index=False)
    network_nucleus_gene_names_list.to_csv(address_network_nucleus_genes, header=None, index=False)

    membrane_list = list(network_membrane_gene_names)
    nucleus_list = list(network_nucleus_gene_names)

    #%% MAIN
    # ANALYSE NETWORK

    # give lists, network, output address to function
    print("Analyzing pathway network using multiple centrality methods.")
    calc_network_centrality_RWR(G_gene, membrane_list, nucleus_list, address_centrality_results)

    #%% MAIN
    # VISUALISE NETWORK

    print("Generating circular plot of network.")
    make_circular_plot_network(list(G_prot.nodes), G_prot, address_network_plot_proteins, start_list_prot)
    make_circular_plot_network(list(G_gene.nodes), G_gene, address_network_plot_genes, start_list_gene)
    print("Program finished.")
