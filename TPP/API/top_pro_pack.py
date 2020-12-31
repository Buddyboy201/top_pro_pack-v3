
import os
import sys
import sqlalchemy
from pathlib import Path
from TPP.API.centroid_protein import CentroidProtein
import json
from shutil import copyfile
from time import perf_counter

def get_config(name, pdb_path, json_path, exclude_backbone, distance_cutoff, ignored_paths):
    config = {
        "name": name,
        "pdb_path": Path(pdb_path).__str__(),
        "json_path": Path(json_path).__str__(),
        "exclude_backbone": exclude_backbone,
        "distance_cutoff": distance_cutoff,
        "ignored_paths": [ Path(file).__str__() for file in ignored_paths ]
    }
    return config

def create_project(config_path, name, pdb_path, json_path, exclude_backbone=False, distance_cutoff=6, ignored_paths=[]):
    #config_path = Path.cwd() / Path("{}_config.json".format(name))
    config = get_config(name=name, pdb_path=pdb_path, json_path=json_path, exclude_backbone=exclude_backbone, distance_cutoff=distance_cutoff, ignored_paths=ignored_paths)
    with open(config_path, "wt") as file:
        json.dump(config, file)


class Project:
    def __init__(self, config_path):
        self._init_project(config_path)
        self.config_path = Path(config_path)
        self.proteins = {}

    def _get_function_perf_decorator(func):
        def inner(self, id, filename):
            start = perf_counter()
            out = func(self, id, filename)
            end = perf_counter()
            print(end - start)
            return out

        return inner

    def _init_project(self, config_path):
        if not Path(config_path).is_file():
            raise Exception("invalid config path: {}".format(Path(config_path)))
        with open(config_path, "rt") as config_file:
            config = json.load(config_file)
            self.distance_cutoff = config["distance_cutoff"]
            self.exclude_backbone = config["exclude_backbone"]
            self.name = config["name"]
            self.pdb_path = Path(config["pdb_path"])
            self.json_path = Path(config["json_path"])
            self.ignored_paths = [ Path(file) for file in config["ignored_paths"] ]
            self.ignore_links = {}
            if not self.pdb_path.is_dir():
                self.pdb_path.mkdir(parents=True)
            if not self.json_path.is_dir():
                self.json_path.mkdir(parents=True)

    def get_protein(self, id):
        try:
            if not self.ignore_links.get(id):
                return self.proteins[id]
            else:
                return None
        except:
            raise Exception("{} is invalid/ignored".format(id))

    @_get_function_perf_decorator
    def load_protein(self, id, file_name):
        file_path = None
        if Path(file_name).suffix == ".json":
            file_path = self.json_path / Path(file_name)
        else:
            file_path = self.pdb_path / Path(file_name)
        if file_path.is_file():
            if Path(file_path) not in self.ignored_paths:
                val = self._init_protein(id, file_path)
                if isinstance(val, Exception):
                    return val
                self.proteins[id] = val
                self.ignore_links[id] = False
                return val
            else:
                self.ignore_links[id] = True
                return None
        else:
            raise Exception("Not a valid {} file".format(file_path.suffix))


    def add_protein(self, file_path):
        if Path(file_path).is_file():
            if Path(file_path).suffix == ".json":
                new_file_path = self.json_path / Path(file_path).name
                if self.json_links.get(name) is not None:
                    raise Exception("{} already taken by {}".format(name, self.json_links.get(name)))
                else:
                    self.json_links[name] = Path(file_path).suffix
            else:
                new_file_path = self.pdb_path / Path(file_path).name

            copyfile(Path(file_path), new_file_path)
        else:
            raise Exception("Not a valid {} file".format(file_path.suffix))

    def add_ignored_path(self, file_path):
        if Path(file_path).is_file():
            self.ignored_paths.append(Path(file_path))
        else:
            raise Exception("{} does not exist".format(Path(file_path)))

    def remove_ignored_path(self, file_path):
        if Path(file_path).is_file():
            self.ignored_paths.remove(Path(file_path))
        else:
            raise Exception("{} does not exist".format(Path(file_path)))

    def load_all_pdbs(self, ids, pdb_filter=None):
        try:
            for pdb_file, id in zip(self.pdb_path.iterdir(), ids):
                print("loading {} as {} ...".format(Path(pdb_file), id))
                val = self.load_protein(id, Path(pdb_file))
                if isinstance(val, Exception):
                    print(val)
                elif isinstance(val, type(None)):
                    print("{} is ignored".format(pdb_file))
                else:
                    print("{} loaded as {}".format(pdb_file, id))
        except:
            raise Exception("All pdbs could not be loaded or handled")

    def load_all_json(self, ids):
        try:
            for json_file, id in zip(self.json_path.iterdir(), ids):
                self.load_protein(id, Path(json_file))
        except:
            raise Exception("All jsons could not be loaded")

    def get_config(self):
        config = {
            "name": self.name,
            "pdb_path": Path(self.pdb_path).__str__(),
            "json_path": Path(self.json_path).__str__(),
            "exclude_backbone": self.exclude_backbone,
            "distance_cutoff": self.distance_cutoff,
            "ignored_paths": self.ignored_paths
        }
        return config

    def list_pdb_files(self):
        return self.pdb_path.glob("*.pdb")

    def list_json_files(self):
        return self.json_path.glob("*.json")

    def list_ignored(self):
        return self.ignored_paths

    def _init_protein(self, id, file_path):
        try:
            P = CentroidProtein(id, file_path, exclude_backbone=self.exclude_backbone)
        except:
            e = sys.exc_info()[0]
            return Exception(e)
        if len(P.residues) > 0 and Path(file_path).suffix != ".json":
            P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
        else:
            return Exception("{} is empty".format(P.name))
        return P










test_code = '''class Project:
    def __init__(self, config_path):
        self._init_project(config_path)
        self.config_path = Path(config_path)
        self.loaded_proteins = {}

    def _update_links(self):
        for pdb_name in self.pdb_links:
            if self.pdb_links.get(pdb_name) not in list(self.pdb_path.iterdir()):
                self.pdb_links.pop(pdb_name, None)
        for json_name in self.json_links:
            if self.json_links.get(json_name) not in list(self.json_path.iterdir()):
                self.json_links.pop(json_name, None)
        self.ignore_links = [Path(file) for file in self.ignore_links if file in list(self.pdb_path.iterdir())+list(self.json_path.iterdir())]
        #if len([Path(file) for file in list(self.pdb_path.iterdir())+list(self.json_path.iterdir()) if Path(file) not in list(self.pdb_links.values())+list(self.json_links.values())+list(self.ignore_links)]) > 0:
            #raise Exception("Not all files have valid identifier providided")


    def _init_project(self, config_path):
        if not Path(config_path).is_file():
            raise Exception("invalid config path: {}".format(Path(config_path)))
        with open(config_path, "rt") as config_file:
            config = json.load(config_file)
            self.distance_cutoff = config["distance_cutoff"]
            self.exclude_backbone = config["exclude_backbone"]
            self.name = config["name"]
            self.pdb_path = Path(config["pdb_path"])
            self.json_path = Path(config["json_path"])
            self.pdb_links = config["pdb_links"]
            self.json_links = config["json_links"]
            self.ignore_links = config["ignore_links"]
            if not self.pdb_path.is_dir():
                self.pdb_path.mkdir(parents=True)
            if not self.json_path.is_dir():
                self.json_path.mkdir(parents=True)
            self._update_links()

    def _init_protein(self, name, file_path):
        try:
            P = CentroidProtein(name, file_path, exclude_backbone=self.exclude_backbone)
        except:
            e = sys.exc_info()[0]
            return Exception(e)
        if len(P.residues) > 0 and Path(file_path).suffix != ".json":
            P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
        else:
            return Exception("{} is empty".format(P.name))
        return P


    def get_config(self):
        config = {
            "name": self.name,
            "pdb_path": self.pdb_path.__str__(),
            "json_path": self.json_path.__str__(),
            "exclude_backbone": self.exclude_backbone,
            "distance_cutoff": self.distance_cutoff,
            "pdb_links": self.pdb_links,
            "json_links": self.json_links,
            "ignore_links": self.ignore_links
        }
        return config

    def add_protein(self, name, file_path):
        if Path(file_path).is_file():
            if Path(file_path).suffix == ".json":
                new_file_path = self.json_path / Path(file_path).name
                if self.json_links.get(name) is not None:
                    raise Exception("{} already taken by {}".format(name, self.json_links.get(name)))
                else:
                    self.json_links[name] = Path(new_file_path).name
            else:
                new_file_path = self.pdb_path / Path(file_path).name
                self.pdb_links[name] = new_file_path

            copyfile(Path(file_path), new_file_path)
        else:
            raise Exception("Not a valid {} file".format(file_path.suffix))

    def load_protein(self, name, priority="pdb"):
        if priority == "pdb":
            if self.pdb_links.get(name) is not None and self.pdb_links.get(name).is_file() and self.pdb_links.get(name) not in self.ignore_links:
                self._init_protein(name, self.pdb_path / Path(self.pdb_links.get(name)))
            else:



class Project:
    def __init__(self, config_path):
        self._init_project(config_path)
        self.config_path = Path(config_path)
        self.proteins = {}

    def get_protein(self, name):
        try:
            if self.ignore_links.get(name) is None:
                return self.proteins[name]
            else:
                return None
        except:
            return Exception("{} not loaded yet".format(name))

    def get_filename_from_name(self, name, priority="pdb"):
        if priority == "json":
            return self.json_links.get(name)
        elif priority == "pdb":
            return self.pdb_links.get(name)
        else:
            return self.ignore_links.get(name)

    def load_protein(self, name, file_name):
        file_path = None
        if Path(file_name).suffix == ".json":
            file_path = self.json_path / Path(file_name)
        else:
            file_path = self.pdb_path / Path(file_name)
        if file_path.is_file():
            self.proteins[name] = self._init_protein(name, file_path)
        else:
            raise Exception("Not a valid {} file".format(file_path.suffix))

    def add_protein(self, name, file_path):
        if Path(file_path).is_file():
            if Path(file_path).suffix == ".json":
                new_file_path = self.json_path / Path(file_path).name
                if self.json_links.get(name) is not None:
                    raise Exception("{} already taken by {}".format(name, self.json_links.get(name)))
                else:
                    self.json_links[name] = Path(file_path).suffix
            else:
                new_file_path = self.pdb_path / Path(file_path).name

            copyfile(Path(file_path), new_file_path)
        else:
            raise Exception("Not a valid {} file".format(file_path.suffix))

    def ignore_protein(self, name):


    def get_config(self):
        config = {
            "name": self.name,
            "pdb_path": self.pdb_path.__str__(),
            "json_path": self.json_path.__str__(),
            "ignore_path": self.ignore_path.__str__(),
            "exclude_backbone": self.exclude_backbone,
            "distance_cutoff": self.distance_cutoff,
            "pdb_links": self.pdb_links,
            "json_links": self.json_links,
            "ignore_links": self.ignore_links
        }
        return config

    def _update_links(self):
        for pdb_name in self.pdb_links:
            if self.pdb_links.get(pdb_name) not in list(self.list_pdbs()):
                self.pdb_links.pop(pdb_name, None)
        for json_name in self.json_links:
            if self.json_links.get(json_name) not in list(self.list_json()):
                self.json_links.pop(json_name, None)
        for ignore_name in self.ignore_links:
            if self.ignore_links.get(ignore_name) not in list(self.list_ignored()):
                self.ignore_links.pop(ignore_name, None)

    def _init_project(self, config_path):
        if not Path(config_path).is_file():
            raise Exception("invalid config path: {}".format(Path(config_path)))
        with open(config_path, "rt") as config_file:
            config = json.load(config_file)
            self.distance_cutoff = config["distance_cutoff"]
            self.exclude_backbone = config["exclude_backbone"]
            self.name = config["name"]
            self.pdb_path = Path(config["pdb_path"])
            self.json_path = Path(config["json_path"])
            self.ignore_path = Path(config["ignore_path"])
            self.pdb_links = config["pdb_links"]
            self.json_links = config["json_links"]
            self.ignore_links = config["ignore_links"]
            if not self.pdb_path.is_dir():
                self.pdb_path.mkdir(parents=True)
            if not self.ignore_path.is_dir():
                self.ignore_path.mkdir(parents=True)
            if not self.json_path.is_dir():
                self.json_path.mkdir(parents=True)
            self._update_links()

    def list_loaded_proteins(self):
        return self.proteins.keys()

    def list_pdbs(self):
        return self.pdb_path.iterdir()

    def list_json(self):
        return self.json_path.iterdir()

    def list_ignored(self):
        return self.ignore_path.iterdir()

    def get_name(self):
        return self.name

    def get_pdb_path(self):
        return self.pdb_path

    def get_json_path(self):
        return self.json_path

    def get_ignore_path(self):
        return self.ignore_path

    def is_mc(self):
        return not self.exclude_backbone

    def get_cutoff(self):
        return self.distance_cutoff

    def _init_protein(self, name, file_path):
        try:
            P = CentroidProtein(name, file_path, exclude_backbone=self.exclude_backbone)
        except:
            e = sys.exc_info()[0]
            return Exception(e)
        if len(P.residues) > 0 and Path(file_path).suffix != ".json":
            P.generate_centroid_cliques(distance_cutoff=self.distance_cutoff)
        else:
            return Exception("{} is empty".format(P.name))
        return P


class TPP_Engine:
    def __init__(self):
        self.projects = {}

    def load_project(self, config_path):
        proj = Project(Path(config_path))
        self.projects[proj.get_name()] = proj'''











































old = '''class TPP_Engine:
    def __init__(self):
        try:
            os.makedirs(Path.home() / Path("top_pro_pack/bin"))
            print("Initializing top_pro_pack data folder at {}".format(Path.home() / Path("top_pro_pack/bin")))
        except:
            print("top_pro_pack data files located at {}".format(Path.home() / Path("top_pro_pack/bin")))
        self.base_lib_path = Path.home() / Path("top_pro_pack")
        self.projects = {}
        self.exclude_backbone = False
        self.distance_cutoff = 6


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
                #self.E.update_static_total_pairs_table(P.get_heatmap_data_centroid())
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
                #self.E.update_static_total_pairs_table((P.get_heatmap_data_centroid()))
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
                #self.E.update_static_total_pairs_table((P.get_heatmap_data_centroid()))
                return P
            else:
                return Exception("{} is empty".format(P.name))
        else:
            print("All processing attempts failed for {}, check provided info and try again".format(name))'''