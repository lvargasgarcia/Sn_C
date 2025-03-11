from itertools import dropwhile
from permutaciones import express_into_adyacent_transpositions
from matrix_utils import decompress_YOR, decompress, print_matrix
from functools import reduce


def get_subpartition(yt):
    partition = [0 for _ in range(yt.n-1)]
    for i in range(len(yt.shape)):
        for j in range(yt.shape[i]):
            if yt.tableau[i][j] == 0 and j < yt.shape[i]:
                partition[i] += 1
    return list(dropwhile(lambda x: x == 0, partition[::-1]))[::-1]

def ordered_1stlevel_tableaux(partition):
    dict = {}
    for i in range(len(partition)):
        for j in range(sum(partition)):
            if j == (partition[i] - 1):
                no_cells_below = True
                k = i + 1
                while no_cells_below and k<len(partition):
                    if partition[k] >= j + 1:
                        no_cells_below = False
                    k += 1
                if no_cells_below:
                    yt = YoungTableau(partition)
                    yt.tableau[i][j] = yt.n
                    dict[(i)] = yt
    resp = []
    ordered_keys = sorted(dict.keys(), key= lambda t: t)
    for k in ordered_keys:
        resp.append(get_subpartition(dict[k]))
    return resp

def youngTableaux_1stlevel(partition):
    resp = []
    for i in range(len(partition)):
        for j in range(sum(partition)):
            if j == (partition[i] - 1):
                no_cells_below = True
                k = i + 1
                while no_cells_below and k<len(partition):
                    if partition[k] >= j + 1:
                        no_cells_below = False
                    k += 1
                if no_cells_below:
                    yt = YoungTableau(partition)
                    yt.tableau[i][j] = yt.n
                    resp.append((yt, (i,j)))
    return resp

def ordered_2ndlevel_tableaux(partition):
    resp = []
    for t in youngTableaux_1stlevel(partition):
        ((yt, (a,b))) = t
        for obj in youngTableaux_1stlevel(get_subpartition(yt)):
            ((yt, (c,d))) = obj
            resp.append(((a,b), (c,d)))
    return sorted(resp, key=lambda t: int(f"{t[0][0]}{t[1][0]}"))

class YoungTableau:

    def __init__(self, partition):
        self.shape = partition
        self.n = sum(partition)
        self.tableau = [[0 for _ in range(self.n)] for _ in range(len(partition))]
        self.height = len(partition)

    def __eq__(self, other):
        return self.tableau == other.tableau and self.shape == other.shape