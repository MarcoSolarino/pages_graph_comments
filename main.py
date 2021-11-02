import os
import json
import networkx as nx


path = 'sequences/sequences/'

# get all sequences names
filenames = [filename for filename in os.listdir(path)]
num_seq = len(filenames)

# build sequences dictionaries
sequences = []
for f in filenames:
    with open(path + f) as json_file:
        sequences.append(json.load(json_file))

print(f'len dictionary: {len(sequences)}')
print(f'id: {((sequences[0])["query"]).split(":")[1]}')
print(f'comment: {((sequences[0])["results"][0])["comment"]}')

if __name__ == '__main__':
    pass

