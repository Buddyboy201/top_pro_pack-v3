
import TPP.API.atom as atom
import TPP.API.residue as residue
import os
import sys
import numpy as np
import scipy.spatial
import networkx as nx
import csv
import math
import time
import matplotlib.pyplot as plt
import statistics
#import streamlit
import pandas as pd
import json
import requests
from pathlib import Path

#python version=3.7.7+


def binary_search(arr, val):
    start = 0
    end = len(arr) - 1
    while True:
        mid = int((start + end) / 2)
        # print(start, end)
        if start >= end and arr[mid] != val: return -1
        if arr[mid] == val:
            return mid
        elif val > arr[mid]:
            start = mid + 1
        else:
            end = mid - 1


def get_one_var_stats(data):
    data.sort()
    x_bar = statistics.mean(data)
    med = statistics.median(data)
    std_dev = statistics.stdev(data, xbar=x_bar)
    mode = statistics.mode(data)
    data_range = data[len(data) - 1] - data[0]
    return x_bar, med, mode, data_range, std_dev


def get_dist(coord1, coord2):
    return math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (coord1[2] - coord2[2]) ** 2)


class CentroidProtein:
    def __init__(self, name, file_path, exclude_backbone=False):
        AAs = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
               "SER", "THR", "TRP", "TYR", "VAL"]

        if Path(file_path).suffix == ".json":
            with open(Path(file_path), "rt") as file:
                data = json.load(file)
                data = data["protein"]
                self.deserialize_res_json(data)
                self.centroid_cliques = []
                for clique in data["centroid_cliques"]:
                    for i in range(len(clique)):
                        clique[i] = self.residues[clique[i]]
                        self.centroid_cliques.append(clique)
                self.centroid_clique_distances = list(data["centroid_clique_distances"])
                self.centroid_clique_frequency = list(data["centroid_clique_freq"])
        else:
            self.name = name
            self.exclude_backbone = exclude_backbone
            self.file_path = file_path
            self.residues = {}
            self.centroid_cliques = []
            # self.centroid_graph = None
            self.centroid_clique_distances = []
            self.centroid_clique_frequency = []
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
                        # atom_id = int(line[6:11].strip(" "))
                        old_res_id = res_id
                        res_id = res_count
                        atom_id = atom_count
                        atom_count += 1
                        coordx = float(line[30:38].strip(" "))
                        coordy = float(line[38:46].strip(" "))
                        coordz = float(line[46:54].strip(" "))
                        bfactor = float(line[61:66].strip(" "))
                        # symbol = line[76:78].strip(" ")
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

    def deserialize_res_json(self, data):
        self.name = data["name"]
        self.exclude_backbone = data["exclude_backbone"]
        self.file_path = data["file_path"]
        self.residues = {}
        for res in data["residues"]:
            res_name = res["name"]
            res_resid = res["resid"]
            res_atoms = []
            res_centroid = res["centroid"]
            for atm in res["atoms"]:
                atm_symbol = atm["symbol"]
                atm_name = atm["name"]
                atm_id = atm["atomid"]
                atm_coords = tuple(atm["coords"])
                atm_mc_sc = atm["mc_sc"]
                atm_mass = atm["atomic_mass"]
                A = atom.Atom("", "", "", "", load_json=True)
                A.symbol = atm_symbol
                A.name = atm_name
                A.atomid = atm_id
                A.coords = atm_coords
                A.mc_sc = atm_mc_sc
                A.atomic_mass = atm_mass
                res_atoms.append(A)
            R = residue.Residue("", "", "", load_json=True)
            R.name = res_name
            R.resid = res_resid
            R.atoms = res_atoms
            R.centroid = res_centroid
            self.residues[res_resid] = R

    def convert_clique_to_json(self, clique):
        return tuple((res.get_json_dict() for res in clique))

    def convert_clique_to_json_id(self, clique):
        return tuple((res.resid for res in clique))

    def get_json_dict(self):
        start_time = time.time()
        data = {
            "protein":
            {
                "name": self.name,
                "exclude_backbone": self.exclude_backbone,
                "file_path": self.file_path,
                "residues": tuple((self.residues[res_id].get_json_dict() for res_id in self.residues)),
                #"centroid_cliques": tuple((self.convert_clique_to_json(clique) for clique in self.centroid_cliques)),
                "centroid_cliques": tuple((self.convert_clique_to_json_id(clique) for clique in self.centroid_cliques)),
                "centroid_clique_distances": tuple(self.centroid_clique_distances),
                "centroid_clique_freq": tuple(self.centroid_clique_frequency)
            }
        }
        print("json generating time: {}".format(time.time()-start_time))
        return data

    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def generate_centroid_cliques(self, distance_cutoff=6):
        self.updated = True
        centroids = []
        centroid_res = {}
        for i in self.residues:
            centroid = self.residues[i].get_centroid(exclude_backbone=self.exclude_backbone)
            if centroid is not None:
                centroids.append(centroid)
                centroid_res[centroid] = self.residues[i]
        centroids = np.array(centroids)
        tri = scipy.spatial.qhull.Delaunay(centroids)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        #self.centroid_graph = graph
        self.centroid_cliques = list(nx.find_cliques(graph))
        for i in range(len(self.centroid_cliques)):
            for j in range(len(self.centroid_cliques[i])):
                self.centroid_cliques[i][j] = centroid_res[tuple(list(centroids[self.centroid_cliques[i][j]]))]
        self.centroid_cliques = np.array(self.centroid_cliques)
        return self.centroid_cliques

    def get_clique_frequency(self):
        self.updated = True
        if len(self.centroid_clique_frequency) > 0:
            return self.centroid_clique_frequency
        if self.centroid_cliques is None:
            self.generate_centroid_cliques()
        freq_arr = [0, 0, 0, 0, 0, 0, 0]
        for i in self.centroid_cliques:
            freq_arr[len(i)] += 1
        self.centroid_clique_frequency = freq_arr
        return freq_arr

    def get_centroid_clique_distances(self):
        self.updated = True
        if len(self.centroid_clique_distances) > 0:
            return self.centroid_clique_distances
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
        self.centroid_clique_distances = distances
        return self.centroid_clique_distances

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
        for clique in self.centroid_cliques:
            for i in range(len(clique)):
                for j in range(i + 1, len(clique)):
                    if ref[clique[i].get_name()] == ref[clique[j].get_name()]:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                    else:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                        arr[ref[clique[j].get_name()]][ref[clique[i].get_name()]] += 1
        return np.array(arr)



old = '''class CentroidProtein:
    def __init__(self, name, file_path, exclude_backbone=False, load_json=False, json_data_file_path=None, download_data=False, data_url=None, raw_data=None):
        self.updated = False
        AAs = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
               "SER", "THR", "TRP", "TYR", "VAL"]
        if load_json:
            file = open(json_data_file_path, "r")
            data = json.load(file)
            data = data["protein"]
            self.deserialize_res_json(data)
            self.centroid_cliques = []
            for clique in data["centroid_cliques"]:
                for i in range(len(clique)):
                    clique[i] = self.residues[clique[i]]
                    self.centroid_cliques.append(clique)
            self.centroid_clique_distances = list(data["centroid_clique_distances"])
            self.centroid_clique_frequency = list(data["centroid_clique_freq"])
        elif download_data:
            pdb_file = requests.get(data_url).text.split("\n")
            self.updated = True
            self.name = name
            self.exclude_backbone = exclude_backbone
            self.file_path = file_path
            self.residues = {}
            self.centroid_cliques = []
            # self.centroid_graph = None
            self.centroid_clique_distances = []
            self.centroid_clique_frequency = []
            atom_count = 0
            res_count = -1
            prev_res = -1
            for line in pdb_file:
                original_res_name = line[16:20].strip(" ")
                #print(original_res_name)
                if line[0:4] == "ATOM" and (original_res_name in AAs or original_res_name[1:] in AAs):
                    res_name = line[17:20].strip(" ")
                    res_id = int(line[22:26].strip(" "))
                    if prev_res != res_id:
                        prev_res = res_id
                        res_count += 1
                    # atom_id = int(line[6:11].strip(" "))
                    old_res_id = res_id
                    res_id = res_count
                    atom_id = atom_count
                    atom_count += 1
                    coordx = float(line[30:38].strip(" "))
                    coordy = float(line[38:46].strip(" "))
                    coordz = float(line[46:54].strip(" "))
                    bfactor = float(line[61:66].strip(" "))
                    symbol = line[76:78].strip(" ")
                    atom_name = line[12:16].strip(" ")
                    coords = (coordx, coordy, coordz)
                    chain = str(line[21].strip(" "))
                    atm = atom.Atom(symbol, atom_name, atom_id, coords, bfactor)
                    if self.residues.get(res_id) is None:
                        self.residues[res_id] = residue.Residue(res_name, res_id, [atm], chain, old_resid=old_res_id)
                    else:
                        self.residues[res_id].add_atom(atm)
        elif raw_data != None:
            self.updated = True
            self.name = name
            self.exclude_backbone = exclude_backbone
            self.file_path = file_path
            self.residues = {}
            self.centroid_cliques = []
            # self.centroid_graph = None
            self.centroid_clique_distances = []
            self.centroid_clique_frequency = []
            atom_count = 0
            res_count = -1
            prev_res = -1
            for line in raw_data:
                original_res_name = line[16:20].strip(" ")
                #print(original_res_name)
                if line[0:4] == "ATOM" and (original_res_name in AAs or original_res_name[1:] in AAs):
                    res_name = line[17:20].strip(" ")
                    #print(res_name)
                    res_id = int(line[22:26].strip(" "))
                    if prev_res != res_id:
                        prev_res = res_id
                        res_count += 1
                    # atom_id = int(line[6:11].strip(" "))
                    old_res_id = res_id
                    res_id = res_count
                    atom_id = atom_count
                    atom_count += 1
                    coordx = float(line[30:38].strip(" "))
                    coordy = float(line[38:46].strip(" "))
                    coordz = float(line[46:54].strip(" "))
                    bfactor = float(line[61:66].strip(" "))
                    symbol = line[76:78].strip(" ")
                    atom_name = line[12:16].strip(" ")
                    coords = (coordx, coordy, coordz)
                    chain = str(line[21].strip(" "))
                    atm = atom.Atom(symbol, atom_name, atom_id, coords, bfactor)
                    if self.residues.get(res_id) is None:
                        self.residues[res_id] = residue.Residue(res_name, res_id, [atm], chain, old_resid=old_res_id)
                    else:
                        self.residues[res_id].add_atom(atm)
        else:
            self.updated = True
            self.name = name
            self.exclude_backbone = exclude_backbone
            self.file_path = file_path
            self.residues = {}
            self.centroid_cliques = []
            # self.centroid_graph = None
            self.centroid_clique_distances = []
            self.centroid_clique_frequency = []
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
                        # atom_id = int(line[6:11].strip(" "))
                        old_res_id = res_id
                        res_id = res_count
                        atom_id = atom_count
                        atom_count += 1
                        coordx = float(line[30:38].strip(" "))
                        coordy = float(line[38:46].strip(" "))
                        coordz = float(line[46:54].strip(" "))
                        bfactor = float(line[61:66].strip(" "))
                        symbol = line[76:78].strip(" ")
                        atom_name = line[12:16].strip(" ")
                        coords = (coordx, coordy, coordz)
                        chain = str(line[21].strip(" "))
                        atm = atom.Atom(symbol, atom_name, atom_id, coords, bfactor)
                        if self.residues.get(res_id) is None:
                            self.residues[res_id] = residue.Residue(res_name, res_id, [atm], chain, old_resid=old_res_id)
                        else:
                            self.residues[res_id].add_atom(atm)


    def deserialize_res_json(self, data):
        self.name = data["name"]
        self.exclude_backbone = data["exclude_backbone"]
        self.file_path = data["file_path"]
        self.residues = {}
        for res in data["residues"]:
            res_name = res["name"]
            res_resid = res["resid"]
            res_atoms = []
            res_centroid = res["centroid"]
            for atm in res["atoms"]:
                atm_symbol = atm["symbol"]
                atm_name = atm["name"]
                atm_id = atm["atomid"]
                atm_coords = tuple(atm["coords"])
                atm_mc_sc = atm["mc_sc"]
                atm_mass = atm["atomic_mass"]
                A = atom.Atom("", "", "", "", load_json=True)
                A.symbol = atm_symbol
                A.name = atm_name
                A.atomid = atm_id
                A.coords = atm_coords
                A.mc_sc = atm_mc_sc
                A.atomic_mass = atm_mass
                res_atoms.append(A)
            R = residue.Residue("", "", "", load_json=True)
            R.name = res_name
            R.resid = res_resid
            R.atoms = res_atoms
            R.centroid = res_centroid
            self.residues[res_resid] = R

    def convert_clique_to_json(self, clique):
        return tuple((res.get_json_dict() for res in clique))

    def convert_clique_to_json_id(self, clique):
        return tuple((res.resid for res in clique))

    def get_json_dict(self):
        start_time = time.time()
        data = {
            "protein":
            {
                "name": self.name,
                "exclude_backbone": self.exclude_backbone,
                "file_path": self.file_path,
                "residues": tuple((self.residues[res_id].get_json_dict() for res_id in self.residues)),
                #"centroid_cliques": tuple((self.convert_clique_to_json(clique) for clique in self.centroid_cliques)),
                "centroid_cliques": tuple((self.convert_clique_to_json_id(clique) for clique in self.centroid_cliques)),
                "centroid_clique_distances": tuple(self.centroid_clique_distances),
                "centroid_clique_freq": tuple(self.centroid_clique_frequency)
            }
        }
        print("json generating time: {}".format(time.time()-start_time))
        return data

    def save_protein_data(self):   # DEPRECATED - USE TPP_Engine() PROJECT SYSTEM INSTEAD
        if not os.path.exists(r"{}\top_pro_pack_logs\{}".format(os.getcwd(), self.name)):
            os.makedirs(r"{}\top_pro_pack_logs\{}".format(os.getcwd(), self.name))
        file = open(r"{}\top_pro_pack_logs\{}\protein_data.json".format(os.getcwd(), self.name), "w")
        json.dump(self.get_json_dict(), file)
        file.close()

    def deserialize_json(self, data):
        pass


    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def generate_centroid_cliques(self, distance_cutoff=6):
        self.updated = True
        centroids = []
        centroid_res = {}
        for i in self.residues:
            centroids.append(self.residues[i].get_centroid())
            centroid_res[self.residues[i].get_centroid()] = self.residues[i]
        centroids = np.array(centroids)
        tri = scipy.spatial.qhull.Delaunay(centroids)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if get_dist(centroids[edge[0]], centroids[edge[1]]) <= distance_cutoff: edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        #self.centroid_graph = graph
        self.centroid_cliques = list(nx.find_cliques(graph))
        for i in range(len(self.centroid_cliques)):
            for j in range(len(self.centroid_cliques[i])):
                self.centroid_cliques[i][j] = centroid_res[tuple(list(centroids[self.centroid_cliques[i][j]]))]
        self.centroid_cliques = np.array(self.centroid_cliques)
        return self.centroid_cliques

    def get_clique_frequency(self):
        self.updated = True
        if len(self.centroid_clique_frequency) > 0:
            return self.centroid_clique_frequency
        if self.centroid_cliques is None:
            self.generate_centroid_cliques()
        freq_arr = [0, 0, 0, 0, 0, 0, 0]
        for i in self.centroid_cliques:
            freq_arr[len(i)] += 1
        self.centroid_clique_frequency = freq_arr
        return freq_arr

    def get_centroid_clique_distances(self):
        self.updated = True
        if len(self.centroid_clique_distances) > 0:
            return self.centroid_clique_distances
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
        self.centroid_clique_distances = distances
        return self.centroid_clique_distances

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
        for clique in self.centroid_cliques:
            for i in range(len(clique)):
                for j in range(i + 1, len(clique)):
                    if ref[clique[i].get_name()] == ref[clique[j].get_name()]:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                    else:
                        arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                        arr[ref[clique[j].get_name()]][ref[clique[i].get_name()]] += 1
        return np.array(arr)'''


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
    for clique in cliques:
        for i in range(len(clique)):
            for j in range(i + 1, len(clique)):
                if ref[clique[i].get_name()] == ref[clique[j].get_name()]:
                    arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                else:
                    arr[ref[clique[i].get_name()]][ref[clique[j].get_name()]] += 1
                    arr[ref[clique[j].get_name()]][ref[clique[i].get_name()]] += 1
    return np.array(arr)


