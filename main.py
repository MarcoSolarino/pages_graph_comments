import os
import json
import networkx as nx
import re
import numpy as np
import matplotlib.pyplot as plt
import sys


path = 'sequences/sequences/'


def get_sequences(limit):
    # get all sequences names
    filenames = np.array([filename for filename in os.listdir(path)])
    if limit is not None:
        filenames = filenames[0:limit]
    # build sequences dictionaries
    sequences = np.empty(0)
    for f in filenames:
        with open(path + f) as json_file:
            sequences = np.append(sequences, json.load(json_file))
    return sequences


def generate_graph(limit=100):
    sequences = get_sequences(limit)
    graph = nx.Graph()  # empty directed graph

    ids = np.array([re.search('A[0-9]{6}', (s["query"])).group() for s in sequences])  # get all ids
    graph.add_nodes_from(ids)  # add all nodes at once

    edges = find_edges(ids, sequences)
    graph.add_edges_from(edges)  # add all edges at once
    return graph


def find_edges(ids, sequences):
    edges = np.empty((0, 2))
    for s in sequences:
        qid = re.search('A[0-9]{6}', (s["query"])).group()
        if 'comment' in s['results'][0]:
            comments = (s["results"][0])["comment"]
            for c in comments:
                references = re.findall('A[0-9]{6}', c)
                for n in references:
                    # check if the reference is a node in the graph
                    if n in ids:
                        edges = np.append(edges, [[qid, str(n)]], axis=0)
    return edges


def bron_kerbosch(graph, p, r=np.empty(0), x=np.empty(0)):
    if len(p) == 0 and len(x) == 0:
        if len(r) > 2:
            print(f'Maximal clique: {r}')
        return

    for v in p:
        new_r = np.union1d(r, v)
        new_p = np.intersect1d(p, np.array(list(graph.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(graph.neighbors(v))))
        bron_kerbosch(graph, new_p, new_r, new_x)
        p = np.delete(p, np.where(p == v))
        x = np.union1d(x, v)
    return


def find_pivot(graph, p, x):
    pux = np.union1d(p, x)
    num_neighbors = 0
    pivot = None
    for node in pux:
        if len(list(graph.neighbors(node))) > num_neighbors:
            pivot = node
            neighbors_p = np.intersect1d(p, list(graph.neighbors(node)))  # get neighbors of node in p
            num_neighbors = len(neighbors_p)
    return pivot


def bron_kerbosch_with_pivot(graph, p, r=np.empty(0), x=np.empty(0)):
    if len(p) == 0 and len(x) == 0:
        if len(r) > 2:
            print(f'Maximal clique: {r}')
        return

    pivot = find_pivot(graph, p, x)
    for v in np.delete(p, np.where(p == pivot)):
        new_p = np.intersect1d(p, np.array(list(graph.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(graph.neighbors(v))))
        bron_kerbosch(graph, new_p, np.union1d(r, v), new_x)
        p = np.delete(p, np.where(p == v))
        x = np.union1d(x, v)
    return


if __name__ == '__main__':
    # sys.setrecursionlimit(3000)
    # g = generate_graph(None)
    # print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
    # maximal_c1 = nx.find_cliques(g)
    # maximal_c1 = [mc for mc in maximal_c1 if len(mc) > 2]
    # bron_kerbosch(p=list(g.nodes()))
    # print(len(list(maximal_c1)))
    random_graph = nx.fast_gnp_random_graph(8, 0.65)
    nx.draw(random_graph, with_labels=True)
    plt.show()
    print('bk1:')
    bron_kerbosch(graph=random_graph, p=list(random_graph.nodes()))
    print('\nbk2:')
    bron_kerbosch_with_pivot(graph=random_graph, p=list(random_graph.nodes()))
