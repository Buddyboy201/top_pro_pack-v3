import sys
sys.path.insert(1, r'C:\Users\aprak\PycharmProjects\TopProPack_v2_2\API')

import os
import API.centroid_protein as centroid_protein
import numpy as np
import pandas as pd
import statistics as stat
import math

class Energy:
    def __init__(self):
        self.STATIC_TOTAL_PAIRS_TABLE = pd.DataFrame(
            np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("int32"))
        self.STATIC_EPAIR_TABLE = pd.DataFrame(
            np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("float64"))
        self.up_to_date = False
        self.ref = {
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

    def initialize_static_total_pairs_table(self):
        pass

    def initialize_static_epair_table(self):
        pass

    def update_static_total_pairs_table(self, protein_pairs_matrix):
        self.STATIC_TOTAL_PAIRS_TABLE = pd.DataFrame(np.add(self.STATIC_TOTAL_PAIRS_TABLE.values.astype("int32"),
                                                            protein_pairs_matrix.astype("int32")).astype('int32'))

    def update_epair_values(self):
        epair_heat_map = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("float64")
        M_const = self.get_M(self.STATIC_TOTAL_PAIRS_TABLE)
        for aa_i in self.ref:
            for aa_j in self.ref:
                e_pair = self.get_epair(self.ref[aa_i], self.ref[aa_j], M_const, self.STATIC_TOTAL_PAIRS_TABLE)
                if aa_i == aa_j:
                    epair_heat_map[self.ref[aa_i]][self.ref[aa_j]] = e_pair
                else:
                    epair_heat_map[self.ref[aa_i]][self.ref[aa_j]] = e_pair
                    epair_heat_map[self.ref[aa_j]][self.ref[aa_i]] = e_pair
        self.STATIC_EPAIR_TABLE = pd.DataFrame(epair_heat_map)

    def get_static_total_pairs_table(self):
        return self.STATIC_TOTAL_PAIRS_TABLE

    def get_static_epair_table(self):
        return self.STATIC_EPAIR_TABLE

    def get_epair(self, aa_i, aa_j):
        return self.STATIC_EPAIR_TABLE.iloc[aa_i, aa_j]

    def get_M(self, aa_heat_map):
        total = 0
        for i in range(20):
            total += aa_heat_map.iloc[:, i].sum()
        return total

    def M_single(self, i, aa_heat_map):
        return aa_heat_map.iloc[:, i].sum()

    def M_E(self, i, j, M_const, aa_heat_map):
        if M_const == 0: return None
        return (self.M_single(i, aa_heat_map) * self.M_single(j, aa_heat_map)) / (M_const ** 2)

    def M_pair(self, i, j, M_const, aa_heat_map):
        if M_const == 0: return None
        return (aa_heat_map.iloc[j, i]) / M_const

    def get_epair(self, i, j, M_const, aa_heat_map):
        m_pair = self.M_pair(i, j, M_const, aa_heat_map)
        m_e = self.M_E(i, j, M_const, aa_heat_map)
        if M_const == 0 or m_pair == None or m_e == None or m_e == 0 or m_pair/m_e <= 0: return 0
        return -math.log( m_pair/ m_e , math.e)




