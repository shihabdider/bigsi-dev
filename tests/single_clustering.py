from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt
import numpy as np
import json


def load_json(path):
    with open(path) as f:
        data = json.load(f)
        return np.asarray(data)


X = load_json('./data/c_elegans_c_briggsae_dist.json')
# print(X)
# X = [[i] for i in list(np.random.randn(190))]

Z = linkage(X, 'single')
n = len(X)
distances = np.array(Z)[:, 2]
print(distances)
print('sum of branch lengths: ', distances.sum())

fig = plt.figure(figsize=(25, 10))
dn = dendrogram(Z)

plt.show()
# plt.savefig('human_slnk_tree.png')


