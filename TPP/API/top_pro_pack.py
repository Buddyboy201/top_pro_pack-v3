

import TPP.API.centroid_protein as centroid_protein
import TPP.API.atom as atom
import TPP.API.residue as residue
import TPP.API.energy as energy
import os
import sys
import json
import time
import TPP.API.visualizer as visualizer
import requests
import sqlalchemy
import numpy as np
import Bio
#TODO: gather 100 sample proteins, do distance distribution analysis, send report
#TODO: use same 100 proteins to generate benchmark runtime speeds for processing proteins manually v. loading from indv. json files/loading from bulk json file

class TPP_Engine:
    def __init__(self, exclude_backbone=False, distance_cutoff=6):
        #os.chdir(r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2")
        #self.STATIC_TOTAL_PAIRS_TABLE = None
        #self.STATIC_EPAIR_TABLE = None
        self.exclude_backbone = exclude_backbone
        self.distance_cutoff = distance_cutoff
        self.proteins = []
        self.E = energy.Energy()
        self.projects = {}
        if "top_pro_pack_logs" not in os.listdir(os.getcwd()):
            os.mkdir("top_pro_pack_logs")

    def load_all_saved_data(self):
        files = os.listdir(os.getcwd()+"\\top_pro_pack_logs")
        #times = []
        #sizes = []
        #start_time = time.time()
        #last_time = start_time
        for name in files:
            #curr_time = time.time()
            self.load_protein(name)
            #diff = curr_time-last_time
            #times.append(diff)
            #last_time = curr_time
            #size = os.stat(os.getcwd()+"\\top_pro_pack_logs\\{}\\protein_data.json".format(name)).st_size
            #sizes.append(size)
            #print(name, curr_time-start_time, diff, size)
        #print("avg. rate: {} | avg. file size: {}".format(len(self.proteins)/(curr_time - start_time), sum(sizes)/float(len(self.proteins))))

    def engine_shutdown(self):
        for i in range(len(self.proteins)):
            if self.proteins[i].updated:
                start_time = time.time()
                self.proteins[i].save_protein_data()
                print("total saving time: {}".format(time.time()-start_time))

    def save_project(self, project_name):
        pass

    def add_protein(self, project_name, name, file_path, json_load=True, data_load=True, data_url="https://files.rcsb.org/download/{}", raw_data=None):
        out = self.init_protein(name, file_path, json_load=json_load, data_load=data_load, data_url=data_url.format(name), raw_data=raw_data)
        if type(out) is Exception: print(out)
        else: self.projects[project_name].append(out)

    def add_dataset(self, project_name, proteins, modifers={"json_load": True, "data_load": True, "data_url": "https://files.rcsb.org/download/{}", "raw_data": None}):
        prev_pdb = ""
        for pdb, file_path in proteins:
            if prev_pdb != pdb: self.add_protein(project_name, pdb, file_path, json_load=modifers["json_load"], data_load=modifers["data_load"], data_url=modifers["data_url"], raw_data=modifers["raw_data"])
            prev_pdb = pdb

    def create_new_project(self, name = "project_{}", proteins=None, modifers={"json_load": True, "data_load": True, "data_url": "https://files.rcsb.org/download/{}", "raw_data": None}):
        if name == "project_{}":
            name = name.format(len(self.projects)+1)
        print("Creating new project {}".format(name))
        self.projects[name] = []
        if proteins is not None:
            self.add_dataset(name, proteins)
        print("Project {} created!".format(name))

    def load_protein(self, name):
        file_path = os.getcwd() +"\\top_pro_pack_logs\\{}\\protein_data.json".format(name)
        P = centroid_protein.CentroidProtein("", "", load_json=True, json_data_file_path=file_path)
        #self.proteins.append(P)
        self.E.update_static_total_pairs_table(P.get_heatmap_data_centroid())
        return P

    def init_protein(self, name, file_path, json_load=True, data_load=True, data_url="https://files.rcsb.org/download/{}", raw_data=None):
        if name in os.listdir(os.getcwd() + "\\top_pro_pack_logs") and json_load:
            print("Attempting to load {} from JSON".format(name))
            return self.load_protein(name)
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

    def initialize_static_total_pairs_table(self):
        pass

    def initialize_static_epair_table(self):
        pass


