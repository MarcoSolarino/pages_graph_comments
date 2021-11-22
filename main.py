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


def bron_kerbosch(graph, p, r=np.empty(0), x=np.empty(0), verbose=True):
    """
    Executes Bron-Kerbosch algorithm with no pivoting ant prints all the maximal cliques (with repetition).
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
        if len(list(graph.neighbors(node))) > num_neighbors:
            pivot = node
            neighbors_p = np.intersect1d(p, list(graph.neighbors(node)))  # get neighbors of node in p
            num_neighbors = len(neighbors_p)
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


def degeneracy_order(graph):
    """
    Find a degeneracy order of the graph.
    :param graph: networkx graph
    :return: ndarray of the degeneracy ordering
    """
    d = 1
    while d <= graph.number_of_nodes():
        # print(f'order {d}')
        subgraph = graph.copy()
        deg_order = np.empty(0)

        while subgraph is not None:
            # np array that stores node indexes and their num of neighbors
            num_neighbors = np.array(
                [(node, len(list(subgraph.neighbors(node)))) for node in (subgraph.nodes())])
            # deleting nodes with more num_neighbors > d
            num_neighbors = np.delete(num_neighbors, np.where(num_neighbors[:, 1] > d)[0], axis=0)
            if len(num_neighbors) == 0:
                subgraph = None
                deg_order = np.empty(0)
            else:
                num_neighbors = np.sort(num_neighbors, axis=0)
                id, __ = num_neighbors[-1]
                deg_order = np.append(deg_order, id)

                # remove node from graph
                subgraph.remove_node(id)
                if len(list(subgraph.nodes())) == 0:
                    subgraph = None
        if deg_order.any():
            return deg_order
        d += 1


def deg_order_hash(graph):
    hash_neighbors = {}
    num_neighbors = np.array(
        [(len(list(graph.neighbors(node))), node) for node in (graph.nodes())])
    for n in num_neighbors:
        if n[0] in hash_neighbors:
            hash_neighbors[n[0]] = np.append(hash_neighbors[n[0]], [n[1]])
        else:
            hash_neighbors[n[0]] = np.array([n[1]])
    hash_neighbors = sorted(hash_neighbors.items())
    print('ciao')


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


def get_max_clique2(graph):
    """
    Search a max clique between all maximal cliques using Bron-Kerbosch using a degeneracy order.
    :param graph: networkx graph
    :return: ndarray of a max clique
    """
    max_clique = np.empty(0)
    deg_order = degeneracy_order(graph)
    nodes = np.array(list(graph.nodes))
    for v in deg_order:
        v_neighbors = np.array(list(graph.neighbors(v)))
        first_half_nodes, second_half_nodes = np.hsplit(nodes, np.where(nodes == v)[0])
        second_half_nodes = np.delete(second_half_nodes, np.where(second_half_nodes == v))
        p = np.intersect1d(v_neighbors, first_half_nodes)
        x = np.intersect1d(v_neighbors, second_half_nodes)
        clique = get_max_clique(graph, p, v, x)
        if len(clique) > len(max_clique):
            max_clique = clique
    return max_clique


if __name__ == '__main__':
    # generate graph and write it

    # g = generate_graph(0)
    # g.remove_nodes_from(list(nx.isolates(g)))
    # print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
    # nx.readwrite.write_adjlist(g, 'graph/big_graph.adjlist')

    # ------------------------------------------------------------------------------------------------------------------
    # read graph and execute bron-kerbosch

    # g = nx.readwrite.read_adjlist('graph/big_graph.adjlist')
    # print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
    # nodes_to_remove = [n for n in g.nodes() if g.degree(n) < 100]
    # g.remove_nodes_from(nodes_to_remove)
    # print(f'g nodes = {g.number_of_nodes()}  g edges = {g.number_of_edges()}')
    #
    # # random_nodes = random.choices(list(g.nodes()), k=100)
    # # g_100 = nx.subgraph(g, random_nodes)
    # g.remove_nodes_from(list(nx.isolates(g)))
    # print(f'g nodes = {g.number_of_nodes()}  g edges = {g.number_of_edges()}')
    # nx.draw(g, with_labels=True)
    # plt.show()
    #
    # maximal_c1 = nx.find_cliques(g)
    # maximal_c1 = [mc for mc in maximal_c1 if len(mc) > 2]
    # print(len(list(maximal_c1)))
    # print('bk1:')
    # bron_kerbosch(g, list(g.nodes()))
    # print('\nbk2:')
    # bron_kerbosch_with_pivot(g, list(g.nodes()))
    # print('\nbk3:')
    # bron_kerbosch_degeneracy(g)

    # ------------------------------------------------------------------------------------------------------------------
    # generate a random graph of 100 nodes and execute bron-kerbosch

    random_graph = nx.fast_gnp_random_graph(10, 0.3)
    nx.draw(random_graph, with_labels=True)
    plt.show()
    #
    # print('bk1:')
    # start = time.time()
    # bron_kerbosch(graph=random_graph, p=list(random_graph.nodes()))
    # print(f'\nTime: {time.time() - start} s')
    # print('\nbk2:')
    # start = time.time()
    # bron_kerbosch_with_pivot(graph=random_graph, p=list(random_graph.nodes()))
    # print(f'\nTime: {time.time() - start} s')
    # print('\nbk3')
    # start = time.time()
    # bron_kerbosch_degeneracy(random_graph)
    # print(f'\nTime: {time.time() - start} s')
    #
    # mc = get_max_clique(random_graph, p=list(random_graph.nodes()))
    # print(f'\nMax Clique method 1: {mc}')
    #
    # mc = get_max_clique2(random_graph)
    # print(f'\nMax Clique method 2: {mc}')

    deg_order_hash(random_graph)

