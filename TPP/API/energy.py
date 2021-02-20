
import os
import TPP.API.centroid_protein as centroid_protein
import numpy as np
import pandas as pd
import statistics as stat
import math
from itertools import permutations, product

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
            #print(aa_heat_map.iloc[:, i])
            total += aa_heat_map.iloc[:, i].sum() + aa_heat_map.iloc[i, i]
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
    def __init__(self, M, L="ALL"):
        if M not in range(2, 5): raise Exception("Invalid clique dim: {}".format(M))
        self.M = M
        self.L = L
        self.STATIC_TOTAL_GROUPS_TABLE = np.zeros(shape=(20,) * M, dtype="int32")
        self.STATIC_EPAIR_TABLE = np.zeros(shape=(20,) * M, dtype="float64")
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

    def update_static_groups_table(self, counts_matrix):
        self.STATIC_TOTAL_GROUPS_TABLE = np.add(self.STATIC_TOTAL_GROUPS_TABLE, counts_matrix)

    def update_epair_table(self):
        unique_clique_combos = self._get_all_unique_cliques()
        total_clique_counts = self._get_total_group_counts(unique_clique_combos)
        print(total_clique_counts)
        for clique in unique_clique_combos:
            if len(list(set(clique))) == 1:
                self.STATIC_EPAIR_TABLE[self._get_res_codes(clique)] = self._compute_epair(clique, total_clique_counts)
            else:
                e_pair = self._compute_epair(clique, total_clique_counts)
                for combo in permutations(clique):
                    #print(combo, e_pair)
                    self.STATIC_EPAIR_TABLE[self._get_res_codes(combo)] = e_pair

    def _cartesian_product(self, residues):
        if self.M != len(residues):
            raise Exception("Mismatched clique dim with given indexing dim: {}".format(len(residues)))
        if self.M == 2:
            return product(residues, residues, repeat=1)
        elif self.M == 3:
            return product(residues, residues, residues, repeat=1)
        elif self.M == 4:
            return product(residues, residues, residues, residues, repeat=1)

    def get_static_groups_table(self):
        return self.STATIC_TOTAL_GROUPS_TABLE

    def get_counts(self, residues):
        res_codes = tuple(map(lambda x: self.ref[x], residues))
        return self.STATIC_TOTAL_GROUPS_TABLE[res_codes]

    def _get_res_codes(self, residues):
        return tuple(map(lambda x: self.ref[x], residues))

    def get_epair(self, residues):
        res_codes = self._get_res_codes(residues)
        return self.STATIC_EPAIR_TABLE[res_codes]

    def _compute_epair(self, residues, total_counts):
        if self.M == 2:
            return self._compute_epair_2(residues, total_counts)
        elif self.M == 3:
            return self._compute_epair_3(residues, total_counts)
        elif self.M == 4:
            return self._compute_epair_4(residues, total_counts)

    def _get_all_unique_cliques(self):
        residues = list(self.ref.keys())
        keys = list(self.ref.keys())
        for i in range(self.M - 1):
            #print(i)
            new_keys = []
            for res_i in keys:
                for res_j in residues:
                    val = None
                    if isinstance(res_i, type("")):
                        val = [res_i, res_j]
                    else:
                        val = [] + res_i
                        val.append(res_j)
                    val.sort()
                    # print(val)
                    new_keys.append(val)
            new_keys_str = list(set(list(map(lambda x: ";".join(sorted(x)), new_keys))))
            new_keys_str.sort()
            new_keys = list(map(lambda x: x.split(";"), new_keys_str))
            keys = new_keys
        return keys

    def _get_total_group_counts(self, unique_clique_combos):
        total = 0
        for combo in unique_clique_combos:
            total += self.get_counts(combo)
        return total

    def _compute_epair_2(self, residues, total_counts):
        if len(residues) != 2:
            raise Exception("Inavlid residues indexing of dim: {}".format(len(residues)))
        A, B = residues[0], residues[1]
        counts_A, counts_B, counts_AB = sum(self.get_counts([A])), sum(self.get_counts([B])), self.get_counts(residues)
        #total_counts = self._get_total_group_counts()
        P_A = counts_A / (2*total_counts)
        P_B = counts_B / (2*total_counts)
        P_AB = counts_AB / (2*total_counts)
        P_indv = P_A * P_B
        if P_indv == 0 or P_AB == 0:
            return 0
        return -np.log(P_AB / P_indv)


    def _compute_epair_3(self, residues, total_counts):
        if len(residues) != 3:
            raise Exception("Inavlid residues indexing of dim: {}".format(len(residues)))
        A, B, C = residues[0], residues[1], residues[2]
        counts_A, counts_BC, counts_ABC = sum(self.get_counts([A])), sum(self.get_counts([B, C])), self.get_counts(residues)
        P_A = counts_A / (3*total_counts)
        P_BC = counts_BC / (3*total_counts)
        P_ABC = counts_ABC / (3*total_counts)
        P_indv = P_A * P_BC
        if P_indv == 0 or P_ABC == 0:
            return 0
        return -np.log(P_ABC / P_indv)


    def _compute_epair_4(self, residues, total_counts):
        pass


class EnergyND2:
    def __init__(self, M, cliques, L="ALL"):
        self.M = M
        self.L = L
        if M not in range(2, 5): raise Exception("Invalid clique dim: {}".format(M))
        self.M = M
        self.L = L
        self.STATIC_EPAIR_TABLE = np.zeros(shape=(20,) * M, dtype="float64")
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
        if len(cliques[0]) != M:
            raise Exception("Mismatched set clique dim and input clique dim")
        self.total_res = 0
        self.res_hash = {1:{}, 2:{}, 3:{}, 4:{}}
        for clique in cliques:
            self.total_res += len(clique)
        for res in self.ref.keys():
            self.res_hash[1][res] = 0
        if self.M >= 2:

            for combo in product(list(self.ref.keys()), list(self.ref.keys())):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(2).get(val) is None: self.res_hash[2][val] = 0
        if self.M >= 3:
            for combo in product(list(self.ref.keys()), list(self.ref.keys()), list(self.ref.keys())):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(3).get(val) is None: self.res_hash[3][val] = 0
        if self.M == 4:
            for combo in product(list(self.ref.keys()), list(self.ref.keys()), list(self.ref.keys()), list(self.ref.keys())):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(4).get(val) is None: self.res_hash[4][val] = 0

        for clique in cliques:
            for res in clique:
                self.res_hash[1][res] += 1
            for i in range(len(clique)):
                for j in range(i+1, len(clique)):
                    pair = [clique[i], clique[j]]
                    pair.sort()
                    val = ";".join(pair)
                    self.res_hash[2][val] += 1
            if self.M > 3:
                for i in range(len(clique)):
                    for j in range(i + 1, len(clique)):
                        for k in range(j+1, len(clique)):
                            triplet = [clique[i], clique[j], clique[k]]
                            triplet.sort()
                            val = ";".join(triplet)
                            self.res_hash[3][val] += 1
            elif self.M == 3:
                clique.sort()
                val = ";".join(clique)
                self.res_hash[3][val] += 1

    def get_counts(self, residues):
        residues.sort()
        val = ";".join(residues)
        return self.res_hash[len(residues)][val]

    def _get_res_codes(self, residues):
        return tuple(map(lambda x: self.ref[x], residues))

    def update_epair_table(self):
        for combo_str in self.res_hash[self.M]:
            #print(combo_str)
            combo = combo_str.split(";")
            if len(list(set(combo))) == 1:
                self.STATIC_EPAIR_TABLE[self._get_res_codes(combo)] = self._compute_epair(combo)
            else:
                e_pair = self._compute_epair(combo)
                for clique_combo in permutations(combo):
                    self.STATIC_EPAIR_TABLE[self._get_res_codes(clique_combo)] = e_pair

    def _compute_epair_2(self, residues):
        if len(residues) != 2:
            raise Exception("Inavlid residues indexing of dim: {}".format(len(residues)))
        A, B = residues[0], residues[1]
        counts_A, counts_B, counts_AB = self.get_counts([A]), self.get_counts([B]), self.get_counts(residues)
        #total_counts = self._get_total_group_counts()
        P_A = counts_A / self.total_res
        P_B = counts_B / self.total_res
        P_AB = counts_AB / self.total_res
        P_indv = P_A * P_B
        if P_indv == 0 or P_AB == 0:
            return 0
        return -np.log(P_AB / P_indv)


    def _compute_epair_3(self, residues):
        if len(residues) != 3:
            raise Exception("Inavlid residues indexing of dim: {}".format(len(residues)))
        A, B, C = residues[0], residues[1], residues[2]
        counts_A, counts_BC, counts_ABC = self.get_counts([A]), self.get_counts([B, C]), self.get_counts(residues)
        P_A = counts_A / self.total_res
        P_BC = counts_BC / self.total_res
        P_ABC = counts_ABC / self.total_res
        P_indv = P_A * P_BC
        if P_indv == 0 or P_ABC == 0:
            return 0
        return -np.log(P_ABC / P_indv)


    def _compute_epair_4(self, residues, total_counts):
        pass

    def _compute_epair(self, residues):
        if self.M == 2:
            return self._compute_epair_2(residues)
        elif self.M == 3:
            return self._compute_epair_3(residues)
        elif self.M == 4:
            return self._compute_epair_4(residues)

