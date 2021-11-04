import os
import json
import networkx as nx
import re
import numpy as np
import sys


path = 'sequences/sequences/'


def get_sequences(limit=100):
    # get all sequences names
    filenames = [filename for filename in os.listdir(path)]
    if limit is not None:
        filenames = filenames[0:limit]
    # build sequences dictionaries
    sequences = []
    for f in filenames:
        with open(path + f) as json_file:
            sequences.append(json.load(json_file))
    return sequences


def generate_graph(limit=100):
    sequences = get_sequences(limit)
    graph = nx.Graph()  # empty directed graph

    ids = [(s["query"]).split(":")[1] for s in sequences]  # get all ids
    graph.add_nodes_from(ids)  # add all nodes at once

    edges = find_edges(sequences)
    # counter = Counter(edges)
    graph.add_edges_from(edges)  # add all edges at once
    return graph


def find_edges(sequences):
    edges = []
    for s in sequences:
        id = (s["query"]).split(":")[1]
        if 'comment' in s['results'][0]:
            comments = (s["results"][0])["comment"]
            for c in comments:
                numbers_in_comment = re.findall('[0-9]+', c)
                numbers_in_comment = list(filter(filter_id, numbers_in_comment))
                for n in numbers_in_comment:
                    edges.append((id, 'A' + str(n)))
    return edges


def filter_id(n):
    if len(n) == 6:
        return True
    else:
        return False


def bron_kerbosch(p, r=np.empty(0), x=np.empty(0)):
    if len(p) == 0 and len(x) == 0:
        if len(r) > 2:
            print(f'Maximal clique: {r}')
            return
    # p is list(g.nodes()) at first iter.
    for v in p:
        new_p = np.intersect1d(p, np.array(list(g.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(g.neighbors(v))))
        bron_kerbosch(new_p, np.append(r, v), new_x)
        p = np.delete(p, np.where(p == v))
        x = np.append(x, v)


def find_pivot(p, x):
    pux = np.union1d(p, x)
    num_neighbors = 0
    pivot = None
    for node in pux:
        if len(list(g.neighbors(node))) >= num_neighbors:
            pivot = node
            num_neighbors = len(list(g.neighbors(node)))
    return pivot


def bron_kerbosch_with_pivot(p, r=np.empty(0), x=np.empty(0)):
    if len(p) == 0 and len(x) == 0:
        if len(r) > 2:
            print(f'Maximal clique: {r}')
            return
    # p is list(g.nodes()) at first iter.
    pivot = find_pivot(p, x)
    for v in np.delete(p, np.where(p == pivot)):
        new_p = np.intersect1d(p, np.array(list(g.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(g.neighbors(v))))
        bron_kerbosch(new_p, np.append(r, v), new_x)
        p = np.delete(p, np.where(p == v))
        x = np.append(x, v)


if __name__ == '__main__':
    sys.setrecursionlimit(2000)
    g = generate_graph()
    # for n in g.nodes():
    #     neighbors = g.neighbors(n)
    #     for neigh in neighbors:
    #         print(neigh)
    print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
    # subgraph_nodes = np.random.choice(list(g.nodes()), 1000)
    # subgraph = g.subgraph(subgraph_nodes)
    # print(f'subgraph edges: {list(subgraph.edges())}')
    # bron_kerbosch(list(g.nodes()))
    bron_kerbosch_with_pivot(list(g.nodes()))

