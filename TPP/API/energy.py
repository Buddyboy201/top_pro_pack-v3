import numpy as np
from itertools import permutations, product
from scipy.spatial import KDTree
from TPP.API.constants import AA_REF, L_MAP


def get_cliques(project, M=2, L="ALL"):
    if M < 2:
        raise Exception(f"Invalid clique dim: {M}")
    cliques = list()
    for structure_id, P in project.proteins.items():
        for clique in P.centroid_cliques:
            if len(clique) != M:
                continue
            layers = [res.layerinfo for res in clique if res.layerinfo in L_MAP[L]]
            if len(layers) != len(clique):
                continue
            cliques.append(tuple(sorted([res.name for res in clique])))
    return cliques


class EnergyMD:
    def __init__(self, project, M=2, L="ALL"):
        if M not in range(2, 5):
            raise Exception("Invalid clique dim: {}".format(M))
        self.M = M
        self.L = L
        self.STATIC_EPAIR_TABLE = np.zeros(shape=(20,) * M, dtype="float64")
        self.result = None
        cliques = get_cliques(project, M=M, L=L)
        if len(cliques[0]) != M:
            raise Exception("Mismatched set clique dim and input clique dim")

        self.total_res = 0
        self.res_hash = {1: {}, 2: {}, 3: {}, 4: {}}
        for clique in cliques:
            self.total_res += len(clique)
            # NOTE: since the input will be changing to a project, would it be better to use the total num of residues in each structure in the project? Keep implementation consistent for now.
        for res in AA_REF.keys():
            self.res_hash[1][res] = 0
        if self.M >= 2:
            for combo in product(list(AA_REF.keys()), list(AA_REF.keys())):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(2).get(val) is None:
                    self.res_hash[2][val] = 0
        if self.M >= 3:
            for combo in product(
                    list(AA_REF.keys()), list(AA_REF.keys()), list(AA_REF.keys())
            ):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(3).get(val) is None:
                    self.res_hash[3][val] = 0
        if self.M == 4:
            for combo in product(
                    list(AA_REF.keys()),
                    list(AA_REF.keys()),
                    list(AA_REF.keys()),
                    list(AA_REF.keys()),
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
            if self.M > 3:  # upper-bound currently at M=4
                for i in range(len(clique)):
                    for j in range(i + 1, len(clique)):
                        for k in range(j + 1, len(clique)):
                            triplet = [clique[i], clique[j], clique[k]]
                            triplet.sort()
                            val = ";".join(triplet)
                            self.res_hash[3][val] += 1
                # clique.sort() - NOTE: clique already sorted per get_rows_project impl
                val = ";".join(clique)
                self.res_hash[4][val] += 1
            elif self.M == 3:
                # clique.sort() - NOTE: clique already sorted per get_rows_project impl
                val = ";".join(clique)
                self.res_hash[3][val] += 1

    def get_counts(self, residues):
        residues.sort()
        val = ";".join(residues)
        return self.res_hash[len(residues)][val]

    def _get_res_codes(self, residues):
        return tuple(map(lambda x: AA_REF[x], residues))

    def update_epair_table(self):
        for combo_str in self.res_hash[self.M]:
            combo = combo_str.split(";")
            if len(list(set(combo))) == 1:
                self.STATIC_EPAIR_TABLE[
                    self._get_res_codes(combo)
                ] = self._compute_epair(combo)
            else:
                e_pair = self._compute_epair(combo)
                for clique_combo in permutations(combo):
                    self.STATIC_EPAIR_TABLE[self._get_res_codes(clique_combo)] = e_pair
        self.result = self.STATIC_EPAIR_TABLE

    def get_result(self):
        if self.result is None:
            self.update_epair_table()
        return self.result

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

    def _compute_epair_4(self, residues):
        if len(residues) != 4:
            raise Exception("Invalid residues indexing of dim: {}".format(len(residues)))
        A, B, C, D = residues
        counts_A, counts_BCD, counts_ABCD = self.get_counts([A]), self.get_counts([B, C, D]), self.get_counts(residues)
        P_A = counts_A / self.total_res
        P_BCD = counts_BCD / self.total_res
        P_ABCD = counts_ABCD / self.total_res
        P_indv = P_A * P_BCD
        if P_indv == 0 or P_ABCD == 0:
            return 0
        return -np.log(P_ABCD / P_indv)

    def _compute_epair(self, residues):
        if self.M == 2:
            return self._compute_epair_2(residues)
        elif self.M == 3:
            return self._compute_epair_3(residues)
        elif self.M == 4:
            return self._compute_epair_4(residues)


# POTENTIALLY DEPRECATED - NEED TO VALIDATE EnergyMD CLASS
class EnergyND2:
    def __init__(self, M, cliques, L="ALL"):
        self.M = M
        self.L = L
        if M not in range(2, 5):
            raise Exception("Invalid clique dim: {}".format(M))
        self.M = M
        self.L = L
        self.STATIC_EPAIR_TABLE = np.zeros(shape=(20,) * M, dtype="float64")
        self.result = None

        if len(cliques[0]) != M:
            raise Exception("Mismatched set clique dim and input clique dim")
        self.total_res = 0
        self.res_hash = {1: {}, 2: {}, 3: {}, 4: {}}
        for clique in cliques:
            self.total_res += len(clique)
            # NOTE: since the input will be changing to a project, would it be better to use the total num of residues in each structure in the project? Keep implementation consistent for now.
        for res in AA_REF.keys():
            self.res_hash[1][res] = 0
        if self.M >= 2:
            for combo in product(list(AA_REF.keys()), list(AA_REF.keys())):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(2).get(val) is None:
                    self.res_hash[2][val] = 0
        if self.M >= 3:
            for combo in product(
                list(AA_REF.keys()), list(AA_REF.keys()), list(AA_REF.keys())
            ):
                vals = list(combo)
                vals.sort()
                val = ";".join(vals)
                if self.res_hash.get(3).get(val) is None:
                    self.res_hash[3][val] = 0
        if self.M == 4:
            for combo in product(
                list(AA_REF.keys()),
                list(AA_REF.keys()),
                list(AA_REF.keys()),
                list(AA_REF.keys()),
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
            if self.M > 3: # upper-bound currently at M=4
                for i in range(len(clique)):
                    for j in range(i + 1, len(clique)):
                        for k in range(j + 1, len(clique)):
                            triplet = [clique[i], clique[j], clique[k]]
                            triplet.sort()
                            val = ";".join(triplet)
                            self.res_hash[3][val] += 1
                clique.sort()
                val = ";".join(clique)
                self.res_hash[4][val] += 1
            elif self.M == 3:
                clique.sort()
                val = ";".join(clique)
                self.res_hash[3][val] += 1

    def get_counts(self, residues):
        residues.sort()
        val = ";".join(residues)
        return self.res_hash[len(residues)][val]

    def _get_res_codes(self, residues):
        return tuple(map(lambda x: AA_REF[x], residues))

    def update_epair_table(self):
        for combo_str in self.res_hash[self.M]:
            combo = combo_str.split(";")
            if len(list(set(combo))) == 1:
                self.STATIC_EPAIR_TABLE[
                    self._get_res_codes(combo)
                ] = self._compute_epair(combo)
            else:
                e_pair = self._compute_epair(combo)
                for clique_combo in permutations(combo):
                    self.STATIC_EPAIR_TABLE[self._get_res_codes(clique_combo)] = e_pair
        self.result = self.STATIC_EPAIR_TABLE

    def get_result(self):
        if self.result is None:
            self.update_epair_table()
        return self.result

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

    def _compute_epair_4(self, residues):
        if len(residues) != 4:
            raise Exception("Invalid residues indexing of dim: {}".format(len(residues)))
        A, B, C, D = residues
        counts_A, counts_BCD, counts_ABCD = self.get_counts([A]), self.get_counts([B, C, D]), self.get_counts(residues)
        P_A = counts_A / self.total_res
        P_BCD = counts_BCD / self.total_res
        P_ABCD = counts_ABCD / self.total_res
        P_indv = P_A * P_BCD
        if P_indv == 0 or P_ABCD == 0:
            return 0
        return -np.log(P_ABCD / P_indv)

    def _compute_epair(self, residues):
        if self.M == 2:
            return self._compute_epair_2(residues)
        elif self.M == 3:
            return self._compute_epair_3(residues)
        elif self.M == 4:
            return self._compute_epair_4(residues)


def get_new_cen6(cen6):
    return min(max(int(cen6), 1), 8) if cen6 <= 8 else 9

# POTENTIALLY DEPRECATED - NEED TO VALIDATE NEW get_new_cen6 FUNCTION
def get_new_cen6_old(cen6):
    new_cen6 = None
    if cen6 <= 1:
        new_cen6 = 1
    elif cen6 <= 2:
        new_cen6 = 2
    elif cen6 <= 3:
        new_cen6 = 3
    elif cen6 <= 4:
        new_cen6 = 4
    elif cen6 <= 5:
        new_cen6 = 5
    elif cen6 <= 6:
        new_cen6 = 6
    elif cen6 <= 7:
        new_cen6 = 7
    elif cen6 <= 8:
        new_cen6 = 8
    else:
        new_cen6 = 9
    return new_cen6


def get_total_count(all_cen6):
    total_count = sum([len(all_cen6[structure_id]) for structure_id in all_cen6])
    return total_count


def get_counts_res(all_cen6):
    counts_res = {type_: 0 for type_ in AA_REF}
    for structure_id in all_cen6:
        for resid in all_cen6[structure_id]:
            type_ = all_cen6[structure_id][resid]["type"]
            counts_res[type_] += 1
    return counts_res


def get_counts_layer_cen6(all_cen6):
    counts = {layer: {cen6: 0 for cen6 in range(0, 10)} for layer in range(1, 8)}
    for structure_id in all_cen6:
        for resid in all_cen6[structure_id]:
            layer = all_cen6[structure_id][resid]["layer"]
            cen6 = all_cen6[structure_id][resid]["cen6"]
            new_cen6 = get_new_cen6(cen6)
            #if layer is None: # due to lack of layerinfo
            #    layer = 7 # 'ALL' layer default fallback
            counts[layer][new_cen6] += 1
    return counts


def get_counts_res_layer_cen6(all_cen6):
    counts = {type_: {layer: {cen6: 0 for cen6 in range(0, 10)} for layer in range(1, 8)} for type_ in AA_REF}
    for structure_id in all_cen6:
        for resid in all_cen6[structure_id]:
            type_ = all_cen6[structure_id][resid]["type"]
            layer = all_cen6[structure_id][resid]["layer"]
            cen6 = all_cen6[structure_id][resid]["cen6"]
            new_cen6 = get_new_cen6(cen6)
            counts[type_][layer][new_cen6] += 1
    return counts


def P_aa(aa_i, total_count, counts_res):
    # print(f"P({aa_i}) = {counts_res[aa_i]} / {total_count}")
    return counts_res[aa_i] / total_count


def P_aa_L_B(aa_i, L, B, counts_layer_cen6, counts_res_layer_cen6):
    # print(f"P({aa_i} | {L}, {B}) = {counts_res_layer_cen6[aa_i][L][B]} / {counts_layer_cen6[L][B]}")
    return counts_res_layer_cen6[aa_i][L][B] / counts_layer_cen6[L][B]


def get_E_env(aa_i, L, B, total_count, counts_res, counts_layer_cen6, counts_res_layer_cen6):
    if counts_res_layer_cen6[aa_i][L][B] == 0:
        return 5.0
    top = P_aa_L_B(aa_i, L, B, counts_layer_cen6, counts_res_layer_cen6)
    bottom = P_aa(aa_i, total_count, counts_res)
    print(f"E_env_{aa_i}_{L}_{B} = -log({top} / {bottom}) = -log({top / bottom}) = {-np.log(top / bottom)}")
    return -np.log(top / bottom)


def get_all_E_env(all_cen6):
    total_count = get_total_count(all_cen6)
    counts_res = get_counts_res(all_cen6)
    counts_layer_cen6 = get_counts_layer_cen6(all_cen6)
    counts_res_layer_cen6 = get_counts_res_layer_cen6(all_cen6)
    result = {type_: {layer: {new_cen6: 0 for new_cen6 in range(1, 10)} for layer in range(1, 8)} for type_ in AA_REF}
    for type_ in AA_REF:
        for layer in range(1, 8):
            for new_cen6 in range(1, 10):
                PP_bur = get_E_env(type_, layer, new_cen6, total_count, counts_res, counts_layer_cen6,
                                   counts_res_layer_cen6)
                result[type_][layer][new_cen6] = PP_bur
                print(f"MEM_ENV_CEN6 {type_} {layer} {new_cen6} {PP_bur}")
    return result


class EnergyEnvKDTree:
    def __init__(self, project):
        self.project = project
        self.result = None

    def _find_res_cen6(self, centroid, tree, radius=6):
        cen6 = tree.query_ball_point(centroid, radius, return_length=True)
        return max(cen6 - 1, 0)  # subtract out identity point

    def _get_all_cen6(self):
        structures = dict()
        for structure_id in self.project.proteins:
            structures[structure_id] = dict()
            structure = self.project.proteins[structure_id]
            centroid_resids = structure.get_centroid_resids_nonetype_check_only()
            centroids = list(centroid_resids.values())
            tree = KDTree(centroids)
            for resid, centroid in centroid_resids.items():
                structures[structure_id][resid] = {
                    "resid": resid,
                    "cen6": self._find_res_cen6(centroid, tree),
                    "type": structure.residues[resid].name,
                    "layer": structure.residues[resid].layerinfo
                        if structure.residues[resid].layerinfo is not None else 7
                }
        return structures

    def update_all_e_env(self):
        self.result = get_all_E_env(self._get_all_cen6())

    def get_result(self):
        if self.result is None:
            self.update_all_e_env()
        return self.result


def update_resid_centroid_ref(structure, res_centroid_ref):
    resid_centroid_cliques = [[r.resid for r in clique] for clique in structure.centroid_cliques]
    for clique in resid_centroid_cliques:
        for resid in clique:
            res_centroid_ref[resid].update(clique)
    return res_centroid_ref


class EnergyEnvClique:
    def __init__(self, project):
        self.project = project
        self.result = None

    def _find_res_cen6(self, res, res_centroid_ref, radius=6):
        cen6 = len(res_centroid_ref[res.resid])
        return max(cen6 - 1, 0)

    def _get_all_cen6(self):
        structures = dict()
        for structure_id in self.project.proteins:
            structures[structure_id] = dict()
            structure = self.project.proteins[structure_id]
            centroid_resids = structure.get_centroid_resids_nonetype_check_only()
            resid_centroid_ref = dict()
            for resid in centroid_resids:
                resid_centroid_ref[resid] = set()
            resid_centroid_ref = update_resid_centroid_ref(structure, resid_centroid_ref)
            # resid_centroid_cliques = [[r.resid for r in clique] for clique in structure.centroid_cliques]
            for resid in centroid_resids:
                # for clique in resid_centroid_cliques:

                structures[structure_id][resid] = {
                    "resid": resid,
                    "cen6": self._find_res_cen6(structure.residues[resid], resid_centroid_ref),
                    "type": structure.residues[resid].name,
                    "layer": structure.residues[resid].layerinfo
                        if structure.residues[resid].layerinfo is not None else 7
                }
        return structures

    def update_all_e_env(self):
        self.result = get_all_E_env(self._get_all_cen6())

    def get_result(self):
        if self.result is None:
            self.update_all_e_env()
        return self.result
