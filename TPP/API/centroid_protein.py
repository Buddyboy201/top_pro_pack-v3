import TPP.API.atom as atom
import TPP.API.residue as residue
import scipy.spatial
import networkx as nx
import math
from TPP.API.constants import AAs, L_MAP

# TODO: remove filter_bfactor parameter?
# TODO: remove self._check_bfactor_threshold?
# TODO: remove any references to bfactor checks?
# TODO: make self.generate_centroid_cliques into self.update_centroid_cliques()
# TODO: add self.generate_centroid_cliques method and add "caching" check for self.centroid_cliques at generation to
#  avoid recomputing accidentally


class CentroidProtein:
    def __init__(
        self,
        structure_id,
        file_path,
        exclude_backbone=False,
        distance_cutoff=6,
        filter_bfactor=60,  # TODO: change name to something better
        filter_pLDDT=70, # formerly 75 # TODO: change name to something better
        is_alphafold=False
    ):

        self.structure_id = structure_id
        self.exclude_backbone = exclude_backbone
        self.distance_cutoff = distance_cutoff
        self.filter_bfactor = filter_bfactor
        self.filter_pLDDT = filter_pLDDT
        self.is_alphafold = is_alphafold
        self.file_path = file_path
        self.residues = {}
        self.centroid_cliques = None
        self._read_pdb()
        self._update_centroids()

    def _update_centroids(self):
        for res_id in self.residues:
            self.residues[res_id].get_centroid(exclude_backbone=self.exclude_backbone)

    def get_centroid_resids(self, en_nonecheck=True, en_thresholdcheck=True, en_layercheck=True, L="ALL"):
        result = {}
        for resid, res in self.residues.items():
            if en_nonecheck and res.get_centroid(exclude_backbone=self.exclude_backbone) is None:
                continue
            elif en_thresholdcheck and not self._check_threshold(res):
                continue
            elif en_layercheck and not (L_MAP[L][0] <= res.layerinfo <= L_MAP[L][-1]):
                continue
            else:
                result[resid] = res.get_centroid(exclude_backbone=self.exclude_backbone)
        return result

    # DEPRECATED
    # TODO: last skip condition is weird - and Nonecheck should always be on s.t. an exception is handled later for unviable centroids
    def get_centroid_resids_old(self, enable_nonecheck=True, enable_bfactorcheck=True, enable_layercheck=False, enable_pLDDTcheck=True, L="ALL"):
        result = {}
        for resid in self.residues:
            if (enable_nonecheck and not self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone) is not None) or \
                    not self.is_alphafold and enable_bfactorcheck and not self._check_bfactor_threshold(self.residues[resid], bfactor_baseline=self.filter_bfactor) or \
                    enable_layercheck and not L_MAP[L][0] <= self.residues[resid].layerinfo <= L_MAP[L][-1] or \
                    self.is_alphafold and enable_pLDDTcheck and not self._check_pLDDT_threshold(self.residues[resid], pLDDT_baseline=self.filter_pLDDT) or \
                    not (enable_nonecheck or enable_bfactorcheck or enable_layercheck):
                continue
            else:
                result[resid] = self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone)
        return result

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

    def _check_threshold(self, res):
        if self.is_alphafold:
            return res.conf_score > self.filter_pLDDT
        return res.conf_score < self.filter_bfactor

    # DEPRECATED
    def _check_bfactor_threshold(self, res, bfactor_baseline):
        return res.get_bfactor() < bfactor_baseline

    # DEPRECATED
    def _check_pLDDT_threshold(self, res, pLDDT_baseline):
        return res.get_bfactor() > pLDDT_baseline

    def generate_centroid_cliques(self, skip_bfactor_check=False, skip_pLDDT_check=False): # TODO: Need to update with newer centroid compute code?
        def _get_dist(coord1, coord2):
            return math.sqrt(
                (coord1[0] - coord2[0]) ** 2
                + (coord1[1] - coord2[1]) ** 2
                + (coord1[2] - coord2[2]) ** 2
            )

        resid_centroids_map = self.get_centroid_resids(enable_bfactorcheck=not skip_bfactor_check,
                                                       enable_pLDDTcheck=not skip_pLDDT_check)
        centroids = list(resid_centroids_map.values())
        centroid_res = {centroid: self.residues[resid] for resid, centroid in resid_centroids_map.items()}

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
        for res in range(len(self.centroid_cliques)):
            for clique_res in range(len(self.centroid_cliques[res])):
                self.centroid_cliques[res][clique_res] = centroid_res[
                    tuple(centroids[self.centroid_cliques[res][clique_res]])
                ]
        # self.centroid_cliques = self.centroid_cliques  # what's this line even do?
        return self.centroid_cliques


