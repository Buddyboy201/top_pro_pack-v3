
import os
import TPP.API.centroid_protein as centroid_protein
import numpy as np
import pandas as pd
import statistics as stat
import math
from itertools import permutations

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
                e_pair = self._get_epair(self.ref[aa_i], self.ref[aa_j], M_const, self.STATIC_TOTAL_PAIRS_TABLE)
                if aa_i == aa_j:
                    epair_heat_map[self.ref[aa_i]][self.ref[aa_j]] = e_pair
                else:
                    epair_heat_map[self.ref[aa_i]][self.ref[aa_j]] = e_pair
                    epair_heat_map[self.ref[aa_j]][self.ref[aa_i]] = e_pair
        self.STATIC_EPAIR_TABLE = pd.DataFrame(epair_heat_map)

    def update_epair_values2(self, layer_map=None):
        epair_heat_map = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]] * 20).astype("float64")
        N_tot = self.get_M(self.STATIC_TOTAL_PAIRS_TABLE)
        layer_N_tot = None
        if isinstance(layer_map, pd.DataFrame):
            layer_N_tot = self.get_M(layer_map)
        print(N_tot, layer_N_tot)
        for aa_i in self.ref:
            for aa_j in self.ref:
                e_pair = self._get_epair2(self.ref[aa_i], self.ref[aa_j], N_tot, self.STATIC_TOTAL_PAIRS_TABLE, layer_map=layer_map, layer_N_tot=layer_N_tot)
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

    def _get_epair(self, i, j, M_const, aa_heat_map):
        m_pair = self.M_pair(i, j, M_const, aa_heat_map)
        m_e = self.M_E(i, j, M_const, aa_heat_map)
        if M_const == 0 or m_pair == None or m_e == None or m_e == 0 or m_pair/m_e <= 0: return 0
        return -math.log( m_pair/ m_e , math.e)




    def _P_i(self, i, N_tot, aa_heat_map, layer_map=None, layer_N_tot=None):
        if N_tot == 0: return None
        if not isinstance(layer_map, pd.DataFrame):
            #print("RIPPPP")
            return aa_heat_map.iloc[:, i].sum()/float(N_tot)
        #print(layer_N_tot)
        return layer_map.iloc[:, i].sum()/float(layer_N_tot)

    def _P_i_j(self, i, j, N_tot, aa_heat_map, M=100, layer_map=None, layer_N_tot=None):
        if N_tot == 0: return None
        if not isinstance(layer_map, pd.DataFrame): return (aa_heat_map.iloc[j, i]) / N_tot
        return (layer_map.iloc[j, i] + M*aa_heat_map.iloc[j, i]/N_tot)/(layer_N_tot + M)

    def _get_epair2(self, i, j, N_tot, aa_heat_map, layer_map=None, layer_N_tot=None):
        P_indv = (self._P_i(i, N_tot, aa_heat_map, layer_map=layer_map, layer_N_tot=layer_N_tot) * self._P_i(j, N_tot, aa_heat_map,
                                                                                    layer_map=layer_map, layer_N_tot=layer_N_tot))
        P_pair = self._P_i_j(i, j, N_tot, aa_heat_map, layer_map=layer_map, layer_N_tot=layer_N_tot)

        if N_tot == 0 or P_pair == None or P_indv == None or P_indv == 0 or P_pair/P_indv <= 0: return 0
        return -math.log(P_pair / P_indv)

class EnergyND:
    def __init__(self, dim=3):
        self.dim = dim

        #self.STATIC_TOTAL_GROUPS_TABLE = np.array(array).astype("int32")
        self.STATIC_TOTAL_GROUPS_TABLE = np.zeros(shape=(20,)*self.dim, dtype="int32")
        self.STATIC_EPAIR_TABLE = np.zeros(shape=(20,)*self.dim, dtype="float64")
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

    def update_static_total_groups_table(self, protein_groups_matrix):
        self.STATIC_TOTAL_GROUPS_TABLE = np.add(self.STATIC_TOTAL_GROUPS_TABLE,
                                                            protein_groups_matrix)

    def get_static_total_groups_table(self):
        return self.STATIC_TOTAL_GROUPS_TABLE

    def get_static_epair_table(self):
        return self.STATIC_EPAIR_TABLE

    def get_epair(self, residues):
        val = self.STATIC_EPAIR_TABLE
        for i in residues:
            val = val[i]
        return val

    def _get_counts_1(self, arr):
        count = 0
        for i in arr:
            if isinstance(i, type(np.array([]))):
                count += self._get_counts_1(i)
            else:
                count += i
        return count

    def _get_counts_2(self, arr):
        count = 0
        for i in arr:
            if isinstance(i[0], type(np.array([]))):
                count += self._get_counts_2(i)
            else:
                count += sum(i)
        return count

    def _P_A(self, A):
        arr = self.get_static_total_groups_table()[A]
        return self._get_counts_2(arr)

    def _P_2(self, X):
        arr = self.get_static_total_groups_table()
        for i in X:
            arr = arr[i]
        return self._get_counts_2(arr)

    def _P_3(self, X):
        arr = self.get_static_total_groups_table()
        for i in X:
            arr = arr[i]
        return self._get_counts_2(arr)

    def _P_4(self, X):
        arr = self.get_static_total_groups_table()
        for i in X:
            arr = arr[i]
        return self._get_counts_2(arr)

    def _get_epair(self, residues):
        if self.dim == 2:
            return -math.log(self._P_2(residues) / (self._P_A(residues[0]) * self._P_A(residues[1])),
                             math.e)
        elif self.dim == 3:
            return -math.log(self._P_3(residues) / (self._P_2(residues) * self._P_A(residues[1])),
                             math.e)
        elif self.dim == 4:
            pass



