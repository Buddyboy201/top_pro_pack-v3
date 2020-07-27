

import TPP.API.centroid_protein as centroid_protein
import TPP.API.atom as atom
import TPP.API.residue as residue
import TPP.API.energy as energy
import os
import sys
import sqlalchemy
from pathlib import Path


class TPP_Engine:
    def __init__(self):
        try:
            os.makedirs(Path.home() / Path("top_pro_pack/bin"))
            print("Initializing top_pro_pack data folder at {}".format(Path.home() / Path("top_pro_pack/bin")))
        except:
            print("top_pro_pack data files located at {}".format(Path.home() / Path("top_pro_pack/bin")))
        self.base_lib_path = Path.home() / Path("top_pro_pack")
        self.projects = {}


    def add_protein(self, project_name, name, file_path, json_load=True, data_load=True, data_url="https://files.rcsb.org/download/{}", raw_data=None):
        out = self.init_protein(name, file_path, json_load=json_load, data_load=data_load, data_url=data_url.format(name), raw_data=raw_data)
        if type(out) is Exception: print(out)
        else: self.projects[project_name].append(out)

    def add_dataset(self, project_name, proteins, modifers={"json_load": True, "data_load": True, "data_url": "https://files.rcsb.org/download/{}", "raw_data": None}):
        prev_pdb = ""
        for pdb, file_path in proteins:
            if prev_pdb != pdb: self.add_protein(project_name, pdb, file_path, json_load=modifers["json_load"], data_load=modifers["data_load"], data_url=modifers["data_url"], raw_data=modifers["raw_data"])
            prev_pdb = pdb

    def create_new_project(self, name="project_{}", exclude_backbone=False, distance_cutoff=6, proteins=None):
        if name == "project_{}":
            name = name.format(len(self.projects)+1)
        print("Attempting to create new project: {}".format(name))
        try:
            project_path = self.base_lib_path / Path("bin/{}".format(name))
            os.makedirs(project_path)
            self.projects[name] = []
            if proteins is not None:
                self.add_dataset(name, proteins)
            print("Project {} created!".format(name))
        except:
            print("Project {} already exists, cancelling operation".format(name))


    def load_protein_json(self, project, name):
        file_path = self.base_lib_path / Path("bin/{}/{}/data.json".format(project, name))
        P = centroid_protein.CentroidProtein("", "", load_json=True, json_data_file_path=file_path)
        #self.proteins.append(P)
        self.E.update_static_total_pairs_table(P.get_heatmap_data_centroid())
        return P

    def init_protein(self, project, name, file_path, json_load=True, data_load=True, data_url="https://files.rcsb.org/download/{}", raw_data=None):
        if name in os.listdir(os.getcwd() + "\\top_pro_pack_logs") and json_load:
            print("Attempting to load {} from JSON".format(name))
            return self.load_protein_json(project, name)
        elif len(file_path) > 0:
            print("Atempting to process {} from directly from pdb file".format(name))
            try: P = centroid_protein.CentroidProtein(name, file_path, exclude_backbone=self.exclude_backbone)
            except:
                e = sys.exc_info()[0]
                return Exception(e)
            if len(P.residues) > 0:
                P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
                # self.proteins.append(P)
                self.E.update_static_total_pairs_table(P.get_heatmap_data_centroid())
                return P
            else:
                return Exception("{} is empty".format(P.name))
        elif data_load and data_url is not None:
            data_url = data_url.format(name[:4] + ".pdb")
            print("Attempting to download/process {} from RCSB".format(name))
            try:
                P = centroid_protein.CentroidProtein(name, "", exclude_backbone=self.exclude_backbone,
                                                     download_data=data_load, data_url=data_url)
            except sqlalchemy.orm.exc.NoResultFound:
                return Exception("{} does not exist in RCSB database".format(name))
            if len(P.residues) > 0:
                P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
                # self.proteins.append(P)
                self.E.update_static_total_pairs_table((P.get_heatmap_data_centroid()))
                return P
            else:
                return Exception("{} is empty".format(P.name))
        elif raw_data != None:
            print("Atempting to process {} from raw text".format(name))
            try: P = centroid_protein.CentroidProtein(name, "", exclude_backbone=self.exclude_backbone, download_data=data_load, data_url=data_url, raw_data=raw_data)
            except:
                e = sys.exc_info()[0]
                return Exception(e)
            if len(P.residues) > 0:
                P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
                # self.proteins.append(P)
                self.E.update_static_total_pairs_table((P.get_heatmap_data_centroid()))
                return P
            else:
                return Exception("{} is empty".format(P.name))
        else:
            print("All processing attempts failed for {}, check provided info and try again".format(name))