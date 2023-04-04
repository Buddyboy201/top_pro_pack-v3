import TPP.API.atom as atom
import TPP.API.residue as residue
import numpy as np
import scipy.spatial
import networkx as nx
import math
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

# python version=3.7.7+


AAs = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
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
    "CYS": 19,
}


class CentroidProtein:
    def __init__(
        self,
        name,
        file_path,
        exclude_backbone=False,
        distance_cutoff=6,
        filter_bfactor=60,
        filter_tmaf_conf=.75,
        tmaf=False
    ):

        self.name = name
        self.exclude_backbone = exclude_backbone
        self.distance_cutoff = distance_cutoff
        self.filter_bfactor = filter_bfactor
        self.filter_tmaf_conf = filter_tmaf_conf
        self.tmaf = tmaf
        self.file_path = file_path
        self.residues = {}
        self.centroid_cliques = []
        self._read_pdb()
        self.ss = None

    def update_ss(self):
        self.ss = dssp_dict_from_pdb_file(self.file_path, DSSP="mkdssp")

    def _read_pdb(self):
        atom_count = 0
        res_count = -1
        prev_res = -1
        with open(self.file_path) as pdb_file:

            for line in pdb_file:
                original_res_name = line[16:20].strip(" ")
                if line[0:4] == "ATOM" and (
                    original_res_name in AAs or original_res_name[1:] in AAs
                ):
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
                    atm = atom.Atom(symbol, atom_name, atom_id, coords)
                    if self.residues.get(res_id) is None:
                        self.residues[res_id] = residue.Residue(
                            res_name, res_id, [atm], chain, bfactor, old_resid=old_res_id
                        )
                    else:
                        self.residues[res_id].add_atom(atm)

    def _check_bfactor_threshold(self, res, bfactor_baseline):
        return res.get_bfactor() < bfactor_baseline

    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def generate_centroid_cliques(self):
        def _get_dist(coord1, coord2):
            return math.sqrt(
                (coord1[0] - coord2[0]) ** 2
                + (coord1[1] - coord2[1]) ** 2
                + (coord1[2] - coord2[2]) ** 2
            )

        centroids = []
        centroid_res = {}
        for res in self.residues:
            centroid = self.residues[res].get_centroid(
                exclude_backbone=self.exclude_backbone
            )
            if centroid is not None and self._check_bfactor_threshold(
                self.residues[res], bfactor_baseline=self.filter_bfactor
            ):
                centroids.append(centroid)
                centroid_res[centroid] = self.residues[res]
        centroids = np.array(centroids)
        tri = scipy.spatial.qhull.Delaunay(centroids)
        edges = []
        for n in tri.simplices:
            edge = sorted([n[0], n[1]])
            if _get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[2]])
            if _get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[0], n[3]])
            if _get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[2]])
            if _get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[1], n[3]])
            if _get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
            edge = sorted([n[2], n[3]])
            if _get_dist(centroids[edge[0]], centroids[edge[1]]) <= self.distance_cutoff:
                edges.append((edge[0], edge[1]))
        graph = nx.Graph(edges)
        self.centroid_cliques = list(nx.find_cliques(graph))
        for protein in range(len(self.centroid_cliques)):
            for res in range(len(self.centroid_cliques[protein])):
                self.centroid_cliques[protein][res] = centroid_res[
                    tuple(list(centroids[self.centroid_cliques[protein][res]]))
                ]
        self.centroid_cliques = np.array(self.centroid_cliques)
        return self.centroid_cliques


