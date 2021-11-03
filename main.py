import os
import json
import networkx as nx
import re
from collections import Counter


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
    graph = nx.Graph()

    ids = [(s["query"]).split(":")[1] for s in sequences]  # get all ids
    graph.add_nodes_from(ids)  # add all nodes at once

    edges = find_edges(sequences)
    counter = Counter(edges)
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


if __name__ == '__main__':
    g = generate_graph()
    print(f'graph nodes = {g.number_of_nodes()}  graph edges = {g.number_of_edges()}')
