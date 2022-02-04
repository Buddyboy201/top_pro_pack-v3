import TPP.API.centroid_protein as centroid_protein
import numpy as np
import pandas as pd
import math
from itertools import permutations, product


class EnergyND2:
    def __init__(self, M, cliques, L="ALL"):
        self.M = M
        self.L = L
        if M not in range(2, 5):
            raise Exception("Invalid clique dim: {}".format(M))
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
            "CYS": 19,
        }
        if len(cliques[0]) != M:
            raise Exception("Mismatched set clique dim and input clique dim")
        self.total_res = 0
        self.res_hash = {1: {}, 2: {}, 3: {}, 4: {}}
        for clique in cliques:
            self.total_res += len(clique)
        for res in self.ref.keys():
            self.res_hash[1][res] = 0
        if self.M >= 2:

            for combo in product(list(self.ref.keys()), list(self.ref.keys())):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(2).get(val) is None:
                    self.res_hash[2][val] = 0
        if self.M >= 3:
            for combo in product(
                list(self.ref.keys()), list(self.ref.keys()), list(self.ref.keys())
            ):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(3).get(val) is None:
                    self.res_hash[3][val] = 0
        if self.M == 4:
            for combo in product(
                list(self.ref.keys()),
                list(self.ref.keys()),
                list(self.ref.keys()),
                list(self.ref.keys()),
            ):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(4).get(val) is None:
                    self.res_hash[4][val] = 0

        for clique in cliques:
            for res in clique:
                self.res_hash[1][res] += 1
            for i in range(len(clique)):
                for j in range(i + 1, len(clique)):
                    pair = [clique[i], clique[j]]
                    pair.sort()
                    val = ";".join(pair)
                    self.res_hash[2][val] += 1
            if self.M > 3:
                for i in range(len(clique)):
                    for j in range(i + 1, len(clique)):
                        for k in range(j + 1, len(clique)):
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
            # print(combo_str)
            combo = combo_str.split(";")
            if len(list(set(combo))) == 1:
                self.STATIC_EPAIR_TABLE[
                    self._get_res_codes(combo)
                ] = self._compute_epair(combo)
            else:
                e_pair = self._compute_epair(combo)
                for clique_combo in permutations(combo):
                    self.STATIC_EPAIR_TABLE[self._get_res_codes(clique_combo)] = e_pair

    def _compute_epair_2(self, residues):
        if len(residues) != 2:
            raise Exception(
                "Inavlid residues indexing of dim: {}".format(len(residues))
            )
        A, B = residues[0], residues[1]
        counts_A, counts_B, counts_AB = (
            self.get_counts([A]),
            self.get_counts([B]),
            self.get_counts(residues),
        )
        # total_counts = self._get_total_group_counts()
        P_A = counts_A / self.total_res
        P_B = counts_B / self.total_res
        P_AB = counts_AB / self.total_res
        P_indv = P_A * P_B
        if P_indv == 0 or P_AB == 0:
            return 0
        return -np.log(P_AB / P_indv)

    def _compute_epair_3(self, residues):
        if len(residues) != 3:
            raise Exception(
                "Inavlid residues indexing of dim: {}".format(len(residues))
            )
        A, B, C = residues[0], residues[1], residues[2]
        counts_A, counts_BC, counts_ABC = (
            self.get_counts([A]),
            self.get_counts([B, C]),
            self.get_counts(residues),
        )
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
