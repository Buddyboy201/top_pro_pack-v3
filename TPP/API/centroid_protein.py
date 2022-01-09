
import TPP.API.atom as atom
import TPP.API.residue as residue
import numpy as np
import scipy.spatial
import networkx as nx
import math
import time
import json
from pathlib import Path

#python version=3.7.7+


def get_dist(coord1, coord2):
    return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (coord1[2] - coord2[2]) ** 2)

AAs = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
               "SER", "THR", "TRP", "TYR", "VAL"]

ref = {
            "GLY": 0,
            "PRO": 1,
            "ASP": 2,
            "GLU": 3,
            "LYS": 4,
            "ARG": 5,
            "HIS": 6,
            "SER": 7,
            "THR": 8,
            "ASN": 9,
            "GLN": 10,
            "ALA": 11,
            "MET": 12,
            "TYR": 13,
            "TRP": 14,
            "VAL": 15,
            "ILE": 16,
            "LEU": 17,
            "PHE": 18,
            "CYS": 19
        }

class CentroidProtein:
    def __init__(self, name, file_path, exclude_backbone=False, distance_cutoff=6, filter_bfactor=60):


        self.name = name
        self.exclude_backbone = exclude_backbone
        self.distance_cutoff = distance_cutoff
        self.filter_bfactor = filter_bfactor
        self.file_path = file_path
        self.residues = {}
        self.centroid_cliques = []
        self.read_pdb()

    def read_pdb(self):
        atom_count = 0
        res_count = -1
        prev_res = -1
        with open(self.file_path) as pdb_file:

            for line in pdb_file:
                original_res_name = line[16:20].strip(" ")
                if line[0:4] == "ATOM" and (original_res_name in AAs or original_res_name[1:] in AAs):
                    res_name = line[17:20].strip(" ")
                    res_id = int(line[22:26].strip(" "))
                    if prev_res != res_id:
                        prev_res = res_id
                        res_count += 1
                    old_res_id = res_id
                    res_id = res_count
                    atom_id = atom_count
                    atom_count += 1
                    coordx = float(line[30:38].strip(" "))
                    coordy = float(line[38:46].strip(" "))
                    coordz = float(line[46:54].strip(" "))
                    bfactor = float(line[61:66].strip(" "))
                    atom_name = line[12:16].strip(" ")
                    symbol = atom_name[0]
                    coords = (coordx, coordy, coordz)
                    chain = str(line[21].strip(" "))
                    atm = atom.Atom(symbol, atom_name, atom_id, coords, bfactor)
                    if self.residues.get(res_id) is None:
                        self.residues[res_id] = residue.Residue(res_name, res_id, [atm], chain,
                                                                old_resid=old_res_id)
                    else:
                        self.residues[res_id].add_atom(atm)

    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def generate_centroid_cliques(self):
        centroids = []
        centroid_res = {}
        for res in self.residues:
            centroid = self.residues[res].get_centroid(exclude_backbone=self.exclude_backbone)
            if centroid is not None:
                centroids.append(centroid)
                centroid_res[centroid] = self.residues[res]
        centroids = np.array(centroids)
        tri = scipy.spatial.qhull.Delaunay(centroids)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        self.centroid_cliques = list(nx.find_cliques(graph))
        for protein in range(len(self.centroid_cliques)):
            for res in range(len(self.centroid_cliques[protein])):
                self.centroid_cliques[protein][res] = centroid_res[tuple(list(centroids[self.centroid_cliques[protein][res]]))]
        self.centroid_cliques = np.array(self.centroid_cliques)
        return self.centroid_cliques

    def get_clique_frequency(self): # TODO: convert freq_arr returned parameter to np array for better integration with upped level classes
        if self.centroid_cliques is None: # TODO: instead of implicitly taking care of centroid clique gen, have user explicitly gen cliques before calling method by raising exception
            self.generate_centroid_cliques()
        freq_arr = [0, 0, 0, 0, 0, 0, 0]
        for i in self.centroid_cliques:
            freq_arr[len(i)] += 1
        return freq_arr

    def get_centroid_clique_distances(self):
        distances = []
        for i in range(len(self.centroid_cliques)):
            clique = self.centroid_cliques[i]
            coords = []
            for j in clique:
                coords.append(j.get_centroid())
            for x in range(len(coords) - 1):
                for y in range(x + 1, len(coords)):
                    d = get_dist(coords[x], coords[y])
                    distances.append(round(d, 3))
        return distances # TODO: convert distances to np array

    def get_heatmap_data_centroid(self):
        arr = [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]

        for clique in self.centroid_cliques:
            for i in range(len(clique)):
                for j in range(i + 1, len(clique)):
                    if ref[clique[i].get_name()] == ref[clique[j].get_name()]:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                    else:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                        arr[ref[clique[j].get_name()]][ref[clique[i].get_name()]] += 1
        return np.array(arr)


def get_protein_pairs_matrix(cliques):
    arr = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]

    for clique in cliques:
        for i in range(len(clique)):
            for j in range(i + 1, len(clique)):
                if ref[clique[i].get_name()] == ref[clique[j].get_name()]:
                    arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                else:
                    arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                    arr[ref[clique[j].get_name()]][ref[clique[i].get_name()]] += 1
    return np.array(arr)

def get_protein_pairs_matrix_str(cliques):
    arr = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]

    for clique in cliques:
        for i in range(len(clique)):
            for j in range(i + 1, len(clique)):
                if ref[clique[i]] == ref[clique[j]]:
                    arr[ref[clique[i]]][ref[clique[j]]] += 1
                else:
                    arr[ref[clique[i]]][ref[clique[j]]] += 1
                    arr[ref[clique[j]]][ref[clique[i]]] += 1
    return np.array(arr)
