import os
import json
import networkx as nx
import re
import numpy as np
import sys
# from collections import Counter


path = 'sequences/sequences/'


def get_sequences():
    # get all sequences names
    filenames = [filename for filename in os.listdir(path)]

    # build sequences dictionaries
    sequences = []
    for f in filenames:
        with open(path + f) as json_file:
            sequences.append(json.load(json_file))
    return sequences


def generate_graph():
    sequences = get_sequences()
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
        print(f'Maximal clique: {r}')
    # p is list(g.nodes()) at first iter.
    for v in p:
        new_p = np.intersect1d(p, np.array(list(g.neighbors(v))))
        new_x = np.intersect1d(x, np.array(list(g.neighbors(v))))
        bron_kerbosch(new_p, np.append(r, v), new_x)
        p = np.delete(p, np.where(p == v))
        x = np.append(x, v)


if __name__ == '__main__':
    sys.setrecursionlimit(30000)
    g = generate_graph()
    # for n in g.nodes():
    #     neighbors = g.neighbors(n)
    #     for neigh in neighbors:
    #         print(neigh)
    print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
    subgraph_nodes = np.random.choice(list(g.nodes()), 100)
    subgraph = g.subgraph(subgraph_nodes)
    bron_kerbosch(list(subgraph.nodes()))

