# code available at https://github.com/MarcoSolarino/pages_graph_comments

import os
import json
import random

import networkx as nx
import re
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import collections


path = 'sequences/sequences/'


def get_sequences(limit):
    """
    Read the json files of the sequences from the folder.
    :param limit: int number of sequences that will be read
    :return: ndarray of sequences
    """
    # get all sequences names
    filenames = np.array([filename for filename in os.listdir(path)])
    if limit != 0:
        filenames = np.random.choice(filenames, size=limit, replace=False)
    # build sequences dictionaries
    sequences = np.empty(0)
    for f in filenames:
        with open(path + f) as json_file:
            sequences = np.append(sequences, json.load(json_file))
    return sequences


def generate_graph(limit=100):
    """
    Generates graph from sequences given by get_sequences().
    :param limit: integer the numer of sequences taken randomly from all sequences
    :return: networkx graph build from sequences
    """
    sequences = get_sequences(limit)
    graph = nx.Graph()  # empty graph

    ids = np.array([re.search('A[0-9]{6}', (s["query"])).group() for s in sequences])  # get all ids
    graph.add_nodes_from(ids)  # add all nodes at once

    edges = find_edges(ids, sequences)
    graph.add_edges_from(edges)  # add all edges at once
    return graph


def find_edges(ids, sequences):
    """
    Find all references to other sequences in comments section and build edges.
    :param ids: ndarray list of all sequence indexes
    :param sequences: ndarray of all sequences
    :return: ndarray of edges
    """
    edges = np.empty((0, 2))
    for s in sequences:
        qid = re.search('A[0-9]{6}', (s["query"])).group()
        if 'comment' in s['results'][0]:
            comments = (s["results"][0])["comment"]
            for c in comments:
                references = re.findall('A[0-9]{6}', c)
                for n in references:
                    # check if the reference is a node in the graph and if it is not a self reference
                    if n in ids and qid != str(n):
                        edges = np.append(edges, [[qid, str(n)]], axis=0)
    return edges


def find_maximal_clique(graph, n):
    """
    Find a single maximal clique (if it exists) which contains node n. This algorithm works in a similar way of the
    Bron-Kerbosh algorithm, but it is not recursive and stops after one maximal clique is found
    (if there is one that includes node n).
    :param graph: NetworkX graph
    :param n: node of the graph
    :return: np-array of the maximal clique or an empty np-array
    """
    p = np.array(list(graph.neighbors(n)))
    r = np.empty(0)
    r = np.append(r, n)
    while len(p) != 0:
        # choosing one of the nodes having more neighbors in p. This guarantees to find a non-trivial maximal clique
        # containing n (if it exists) but does not guarantee that such maximal clique has the maximum size possible
        nn = find_pivot(graph, p, [])
        nn_neighbors = np.array(list(graph.neighbors(nn)))
        r = np.append(r, nn)
        p = np.intersect1d(p, nn_neighbors)
    if len(r) > 2:
        return r
    else:
        return np.empty(0)


def bron_kerbosch(graph, p, r=np.empty(0), x=np.empty(0), verbose=True):
    """
    Executes Bron-Kerbosch algorithm with no pivoting and prints all the maximal cliques (with repetition).
    :param graph: networkx graph of sequences
    :param p: list of nodes
    :param r: list of nodes candidates to be a maximal clique
    :param x: list of nodes already explored
    :param verbose: bool if True the method prints all te maximal cliques
    :return:
    """
    if len(p) == 0 and len(x) == 0:
        if len(r) > 2 and verbose:
            print(f'Maximal clique: {r}')
        return

    for v in p:
        new_r = np.union1d(r, v)
        new_p = np.intersect1d(p, np.array(list(graph.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(graph.neighbors(v))))
        bron_kerbosch(graph, new_p, new_r, new_x, verbose=verbose)
        p = np.delete(p, np.where(p == v))
        x = np.union1d(x, v)
    return


def find_pivot(graph, p, x):
    """
    Find the pivot using Tomita algorithm
    :param graph: networkx graph of sequences
    :param p: list of nodes
    :param x: list of nodes
    :return: networkx node: pivot
    """
    pux = np.union1d(p, x)
    num_neighbors = 0
    pivot = None
    for node in pux:
        neighbors_in_p = np.intersect1d(p, list(graph.neighbors(node)))  # get the number of neighbors of node in p
        if len(neighbors_in_p) >= num_neighbors:
            pivot = node
            num_neighbors = len(neighbors_in_p)
    return pivot


def bron_kerbosch_with_pivot(graph, p, r=np.empty(0), x=np.empty(0), verbose=True):
    """
    Executes Bron-Kerbosch algorithm using Tomita pivot.
    :param graph: networkx graph of sequences
    :param p: list of nodes
    :param r: list of nodes candidates to be a maximal clique
    :param x: list of nodes already explored
    :param verbose: bool if True the method prints all te maximal cliques
    :return:
    """
    if len(p) == 0 and len(x) == 0:
        if len(r) > 2 and verbose:
            print(f'Maximal clique: {r}')
        return

    pivot = find_pivot(graph, p, x)
    for v in np.setdiff1d(p, list(graph.neighbors(pivot))):
        new_p = np.intersect1d(p, np.array(list(graph.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(graph.neighbors(v))))
        bron_kerbosch(graph, new_p, np.union1d(r, v), new_x)
        p = np.delete(p, np.where(p == v))
        x = np.union1d(x, v)
    return


def build_neighbors_dict(graph):
    """
    build a dictionary (hash) with key = degree of a node, values = list of nodes with key-degree
    :param graph: NetworkX graph
    :return: OrderedDict, the dictionary ordered by key values
    """
    # build empty dictionary
    hash_neighbors = {}
    # build a np-array of couples (degree, node)
    num_neighbors = np.array([(len(list(graph.neighbors(node))), node) for node in (graph.nodes())])
    for n in num_neighbors:
        # the key is the degree of the node, the value is the node itself
        key, value = int(n[0]), n[1]
        if key in hash_neighbors:
            hash_neighbors[key] = np.append(hash_neighbors[key], value)
        else:
            hash_neighbors[key] = np.array([value])
    # ordering the dictionary by its key values
    hash_neighbors = collections.OrderedDict(sorted(hash_neighbors.items()))
    return hash_neighbors


def degeneracy_order(graph):
    """
    This is a greedy algorithm for computing the degeneracy ordering of a graph. At each iterations it choose the node
    with the smallest degree
    :param graph: NetworkX graph
    :return: np-array, the degeneracy order of the graph
    """
    node_degrees = build_neighbors_dict(graph)
    deg_order = np.empty(0)
    degree = 0
    graph_copy = graph.copy()
    # loop until deg_order is complete
    while len(deg_order) < len(list(graph.nodes())):
        # select the smallest key, which is the smallest key that holds a non-empty list of nodes
        smallest_degree = int((list(node_degrees.keys()))[0])
        while smallest_degree not in node_degrees or len(node_degrees[smallest_degree]) == 0:
            smallest_degree += 1
        # update degree which is the maximum of the degrees of the nodes at the time they are removed from the graph
        if degree < smallest_degree:
            degree = smallest_degree
        # pop the last node in the list
        node = (node_degrees[smallest_degree])[-1]
        node_degrees.update({smallest_degree: np.delete(node_degrees[smallest_degree],
                                                        np.where(node_degrees[smallest_degree] == node))})
        # add node to the degeneracy ordering
        deg_order = np.append(deg_order, node)

        # update dictionary for every neighbor of node: move each neighbor in the previous list
        node_neighbors = list(graph_copy.neighbors(node))
        for nn in node_neighbors:
            # the key will be the degree of the node
            key = len(list(graph_copy.neighbors(nn)))
            node_degrees.update({key: np.delete(node_degrees[key], np.where(node_degrees[key] == nn))})

            # reduce the degree of nn by one (move it to the list in the previous position
            key -= 1
            if key in node_degrees:
                node_degrees.update({key: np.append(node_degrees[key], nn)})
            # if the key does not exist we need to add it to the dictionary and sort it again
            else:
                node_degrees[key] = np.array([nn])
                node_degrees = collections.OrderedDict(sorted(node_degrees.items()))
        # remove the node from the graph
        graph_copy.remove_node(node)
    print(f'\ndegeneracy order for graph is {degree}')
    return deg_order


def bron_kerbosch_degeneracy(graph):
    """
    Executes Bron-Kerbosch algorithm using a degeneracy order of the graph
    :param graph: networkx graph
    :return:
    """
    deg_order = degeneracy_order(graph)
    nodes = np.array(list(graph.nodes))
    for v in deg_order:
        v_neighbors = np.array(list(graph.neighbors(v)))
        first_half_nodes, second_half_nodes = np.hsplit(nodes, np.where(nodes == v)[0])
        second_half_nodes = np.delete(second_half_nodes, np.where(second_half_nodes == v))
        p = np.intersect1d(v_neighbors, first_half_nodes)
        x = np.intersect1d(v_neighbors, second_half_nodes)
        bron_kerbosch_with_pivot(graph, p, v, x)


def get_max_clique(graph, p, r=np.empty(0), x=np.empty(0)):
    """
    Search a max clique between all maximal cliques using Bron-Kerbosch with pivoting.
    :param graph: networkx graph
    :param p: list of nodes
    :param r: list of nodes candidates to be a maximal clique
    :param x: list of nodes already explored
    :return: ndarray of a max clique
    """
    if len(p) == 0 and len(x) == 0:
        return r

    max_clique = np.empty(0)
    pivot = find_pivot(graph, p, x)
    for v in np.setdiff1d(p, list(graph.neighbors(pivot))):
        new_p = np.intersect1d(p, np.array(list(graph.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(graph.neighbors(v))))
        clique = get_max_clique(graph, new_p, np.union1d(r, v), new_x)
        if len(clique) > len(max_clique):
            max_clique = clique
        p = np.delete(p, np.where(p == v))
        x = np.union1d(x, v)
    return max_clique


if __name__ == '__main__':
    # generate graph and save it

    # g = generate_graph(0)
    # g.remove_nodes_from(list(nx.isolates(g)))
    # print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
    # nx.readwrite.write_adjlist(g, 'graph/big_graph.adjlist')

    # ------------------------------------------------------------------------------------------------------------------
    # read previously saved graph, this graph is sparse, selecting random subgraphs from it leads to a graph
    # without maximal cliques

    g = nx.readwrite.read_adjlist('graph/big_graph.adjlist')
    print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
    density = g.number_of_edges() / ((g.number_of_nodes() * (g.number_of_nodes() - 1))/2)
    print(f'The density of the graph is {density}')
    print("---------------------------------------------------------------------------------")

    # ------------------------------------------------------------------------------------------------------------------
    # generate a random graph of 100 nodes and execute bron-kerbosch

    print("Using a random graph with 100 nodes:\n")
    random_graph = nx.fast_gnp_random_graph(100, 0.4)
    nx.draw(random_graph, with_labels=True)
    plt.show()

    nodes = np.array(list(random_graph.nodes()))
    # choosing the node with the highest degree to improve the probability to find a non-trivial maximal clique
    pivot_node = find_pivot(random_graph, nodes, np.empty(0))
    print(f' Maximal Clique including node {pivot_node} : {find_maximal_clique(random_graph, pivot_node)}')

    print("---------------------------------------------------------------------------")
    print('Executing Bron-Kerbosch without pivoting:')
    start = time.time()
    bron_kerbosch(graph=random_graph, p=list(random_graph.nodes()))
    print(f'\nTime: {time.time() - start} s')
    print("---------------------------------------------------------------------------")
    print('Executing Bron-Kerbosch with pivoting::')
    start = time.time()
    bron_kerbosch_with_pivot(graph=random_graph, p=list(random_graph.nodes()))
    print(f'\nTime: {time.time() - start} s')
    print("---------------------------------------------------------------------------")
    print('Executing Bron-Kerbosch using a degeneracy order of the nodes:')
    start = time.time()
    bron_kerbosch_degeneracy(random_graph)
    print(f'\nTime: {time.time() - start} s')

    print("---------------------------------------------------------------------------")
    start = time.time()
    mc = get_max_clique(random_graph, p=list(random_graph.nodes()))
    print(f'Max Clique method 1: {mc}')
    print(f'Time: {time.time() - start} s')
