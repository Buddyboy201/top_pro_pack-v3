import TPP.API.atom as atom
import TPP.API.residue as residue
import scipy.spatial
import networkx as nx
import math
from TPP.API.constants import AAs, L_MAP

# TODO: remove filter_bfactor parameter?
# TODO: replace self.name with self.structure_id
# TODO: remove deprecated functions
# TODO: remove self._check_bfactor_threshold?
# TODO: remove unnecessary getters
# TODO: remove any references to bfactor checks?
# TODO: make self.generate_centroid_cliques and internal method
# TODO: add self.generate_centroid_cliques method and add "caching" check for self.centroid_cliques at generation to
#  avoid recomputing accidentally
# TODO: rename tmaf parameter to something more intelligible
# TODO: reformat using formatter

class CentroidProtein:
    def __init__(
        self,
        name,
        file_path,
        exclude_backbone=False,
        distance_cutoff=6,
        filter_bfactor=60,
        filter_pLDDT=70, # formerly 75
        tmaf=False
    ):

        self.name = name
        self.exclude_backbone = exclude_backbone
        self.distance_cutoff = distance_cutoff
        self.filter_bfactor = filter_bfactor
        self.filter_pLDDT = filter_pLDDT
        self.tmaf = tmaf
        self.file_path = file_path
        self.residues = {}
        self.centroid_cliques = None
        self._read_pdb()
        self._update_centroids()

    def _update_centroids(self):
        for res_id in self.residues:
            self.residues[res_id].get_centroid(exclude_backbone=self.exclude_backbone)

    # TODO: last skip condition is weird - and Nonecheck should always be on s.t. an exception is handled later for unviable centroids
    def get_centroid_resids(self, enable_nonecheck=True, enable_bfactorcheck=True, enable_layercheck=False, enable_pLDDTcheck=True, L="ALL"):
        result = {}
        for resid in self.residues:
            if (enable_nonecheck and not self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone) is not None) or \
                    not self.tmaf and enable_bfactorcheck and not self._check_bfactor_threshold(self.residues[resid], bfactor_baseline=self.filter_bfactor) or \
                    enable_layercheck and not L_MAP[L][0] <= self.residues[resid].layerinfo <= L_MAP[L][-1] or \
                    self.tmaf and enable_pLDDTcheck and not self._check_pLDDT_threshold(self.residues[resid], pLDDT_baseline=self.filter_pLDDT) or \
                    not (enable_nonecheck or enable_bfactorcheck or enable_layercheck):
                continue
            else:
                result[resid] = self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone)
        return result

    # DEPRECATED?
    def old_get_centroid_resids(self, L="ALL"):
        return {resid: self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone) for resid in self.residues
                         if self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone) is not None
                         and self._check_bfactor_threshold(self.residues[resid], bfactor_baseline=self.filter_bfactor)
                         and L_MAP[L][0] <= self.residues[resid].layerinfo <= L_MAP[L][-1]}

    # DEPRECATED?
    def get_centroid_resids_nonetype_check_only(self):
        return {resid: self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone) for resid in self.residues
                         if self.residues[resid].get_centroid(exclude_backbone=self.exclude_backbone) is not None}

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

    def _check_pLDDT_threshold(self, res, pLDDT_baseline):
        return res.get_bfactor() > pLDDT_baseline

    def get_name(self):
        return self.name

    def get_file_path(self):
        return self.file_path

    def get_residues(self):
        return self.residues

    def generate_centroid_cliques(self, skip_bfactor_check=False, skip_pLDDT_check=False): # TODO: Need to update with newer centroid compute code
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
        self.centroid_cliques = self.centroid_cliques
        return self.centroid_cliques


