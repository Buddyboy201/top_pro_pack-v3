import sys
from pathlib import Path
from TPP.API.centroid_protein import CentroidProtein
import json
from shutil import copyfile
from time import perf_counter
from TPP.API.verbose import handle_debug


# filter_bfactor: <default baseline, but can set custom value>


def get_config(
    name, pdb_path, exclude_backbone, distance_cutoff, filter_bfactor, ignored_paths, tmaf
):
    config = {
        "name": name,
        "pdb_path": Path(pdb_path).__str__(),
        "exclude_backbone": exclude_backbone,
        "distance_cutoff": distance_cutoff,
        "filter_bfactor": filter_bfactor,  # remove res if any atms fail baseline
        "ignored_paths": [Path(file).__str__() for file in ignored_paths],
        "tmaf": tmaf
    }
    return config


# filter_bfactor default baseline currently temp, will be changed later to more ideal value
def create_project(
    config_path,
    name,
    pdb_path,
    exclude_backbone=False,
    distance_cutoff=6,
    filter_bfactor=60,
    tmaf=False,
    ignored_paths=[],
):
    config = get_config(
        name=name,
        pdb_path=pdb_path,
        exclude_backbone=exclude_backbone,
        distance_cutoff=distance_cutoff,
        filter_bfactor=filter_bfactor,
        ignored_paths=ignored_paths,
        tmaf=tmaf
    )

    with open(config_path, "wt") as file:
        json.dump(config, file)


class Project:
    def __init__(self, config_path):
        self._init_project(config_path)
        self.config_path = Path(config_path)
        self.proteins = {}

    def generate_default_ids(self):
        return [
            f.stem if f not in self.list_ignored() else ""
            for f in self.list_pdb_files()
        ]

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
            self.filter_bfactor = config["filter_bfactor"]
            self.name = config["name"]
            self.pdb_path = Path(config["pdb_path"])
            self.ignored_paths = [Path(file) for file in config["ignored_paths"]]
            self.tmaf = config["tmaf"]
            self.ignore_links = {}
            if not self.pdb_path.is_dir():
                self.pdb_path.mkdir(parents=True)

    def get_protein(self, id):
        try:
            if not self.ignore_links.get(id):
                return self.proteins[id]
            else:
                return None
        except:
            raise Exception("{} is invalid/ignored".format(id))

    # @_get_function_perf_decorator # debugging funtion !!!
    def load_protein(self, id, file_name):
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
            for pdb_file, id in zip(self.list_pdb_files(), ids):
                handle_debug(print, "loading {} as {} ...".format(Path(pdb_file), id))
                try:
                    val = self.load_protein(id, Path(pdb_file))
                    if isinstance(val, Exception):
                        handle_debug(print, val)
                    elif isinstance(val, type(None)):
                        handle_debug(print, "{} is ignored".format(pdb_file))
                    else:
                        handle_debug(print, "{} loaded as {}".format(pdb_file, id))
                except:
                    print("{} could not be loaded".format(pdb_file))
        except:
            raise Exception("All pdbs could not be loaded or handled")

    def get_config(self):
        config = {
            "name": self.name,
            "pdb_path": Path(self.pdb_path).__str__(),
            "exclude_backbone": self.exclude_backbone,
            "distance_cutoff": self.distance_cutoff,
            "filter_bfactor": self.filter_bfactor,
            "ignored_paths": self.ignored_paths,
            "tmaf": self.tmaf
        }
        return config

    def list_pdb_files(self):
        return self.pdb_path.glob("*.pdb")

    def list_ignored(self):
        return self.ignored_paths

    def _init_protein(self, id, file_path):
        try:
            P = CentroidProtein(
                id,
                file_path,
                exclude_backbone=self.exclude_backbone,
                distance_cutoff=self.distance_cutoff,
                filter_bfactor=self.filter_bfactor,
                tmaf=self.tmaf
            )
        except:
            e = sys.exc_info()[0]
            return Exception(e)
        if len(P.residues) > 0:
            P.generate_centroid_cliques()
        else:
            return Exception("{} is empty".format(P.name))
        return P
