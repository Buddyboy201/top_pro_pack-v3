import centroid_protein
import atom
import residue
import energy

class Project:
    def __init__(self, exclude_backbone=False, distance_cutoff=6):
        self.exclude_backbone = exclude_backbone
        self.distance_cutoff = distance_cutoff
        self.proteins = []
        self.E = energy.Energy()
        self.projects = {}

    def save_project(self, project_name):
        pass

    def add_dataset(self, project_name, proteins, modifers={"json_load": True, "data_load": True, "data_url": "https://files.rcsb.org/download/{}"}):
        prev_pdb = ""
        for pdb, file_path in proteins:
            if prev_pdb != pdb: self.add_protein(project_name, pdb, file_path, json_load=modifers["json_load"], data_load=modifers["data_load"], data_url=modifers["data_url"])
            prev_pdb = pdb