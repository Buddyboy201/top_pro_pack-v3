import sys
from pathlib import Path
from TPP.API.centroid_protein import CentroidProtein
import json
from shutil import copyfile
from TPP.API.verbose import handle_debug


# filter_bfactor: <default baseline, but can set custom value>


# TODO: remove filter_bfactor parameter and all references and checks to it?


def get_config(
    name, pdb_path, exclude_backbone, distance_cutoff, filter_bfactor, ignored_paths, is_alphafold, filter_pLDDT
):
    config = {
        "name": name,
        "pdb_path": Path(pdb_path).__str__(),
        "exclude_backbone": exclude_backbone,
        "distance_cutoff": distance_cutoff,
        "filter_bfactor": filter_bfactor,  # remove res if any atms fail baseline
        "ignored_paths": [Path(file).__str__() for file in ignored_paths],
        "is_alphafold": is_alphafold,
        "filter_pLDDT": filter_pLDDT
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
    is_alphafold=False,
    filter_pLDDT=70,
    ignored_paths=tuple(),
):
    config = get_config(
        name=name,
        pdb_path=pdb_path,
        exclude_backbone=exclude_backbone,
        distance_cutoff=distance_cutoff,
        filter_bfactor=filter_bfactor,
        ignored_paths=ignored_paths,
        is_alphafold=is_alphafold,
        filter_pLDDT=filter_pLDDT
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
            self.is_alphafold = config["is_alphafold"]
            self.filter_pLDDT = config["filter_pLDDT"]
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

    def load_protein(self, id, file_path, skip_clique_gen=False, skip_layer_info=True, out_dir=None,
                     skip_bfactor_check=False, skip_pLDDT_check=False):
        # file_path = self.pdb_path / Path(file_name)
        out_path = None
        if not skip_layer_info:
            out_path = out_dir / Path(f"{id}.out")
        if file_path.is_file():
            if Path(file_path) not in self.ignored_paths:
                val = self._init_protein(id, file_path, skip_clique_gen=skip_clique_gen,
                                         skip_layer_info=skip_layer_info, out_path=out_path,
                                         skip_bfactor_check=skip_bfactor_check, skip_pLDDT_check=skip_pLDDT_check)
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

    def load_all_pdbs(self, ids, skip_clique_gen=False, skip_layer_info=True, out_dir=None, pdb_filter=None,
                      skip_bfactor_check=False, skip_pLDDT_check=False):
        if out_dir is None:
            skip_layer_info = True
        try:
            for pdb_file, id in zip(self.list_pdb_files(), ids):
                handle_debug(print, "loading {} as {} ...".format(Path(pdb_file), id))
                try:
                    val = self.load_protein(id, Path(pdb_file), skip_clique_gen=skip_clique_gen,
                                            skip_layer_info=skip_layer_info, out_dir=out_dir,
                                            skip_bfactor_check=skip_bfactor_check, skip_pLDDT_check=skip_pLDDT_check)
                    if isinstance(val, Exception):
                        handle_debug(print, val)
                    elif isinstance(val, type(None)):
                        handle_debug(print, "{} is ignored".format(pdb_file))
                    else:
                        handle_debug(print, "{} loaded as {}".format(pdb_file, id))
                except:
                    e = sys.exc_info()
                    print(e)
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
            "is_alphafold": self.is_alphafold
        }
        return config

    def list_pdb_files(self):
        return self.pdb_path.glob("*.pdb")

    def list_ignored(self):
        return self.ignored_paths

    def _process_out_file(self, P, out_path, min_hydrophobic_residues=34, residue_baseline=30):
        def _get_layer_resid(resid, ref):
            return ref[resid + 1]

        def _get_cen6_resid(resid, ref):
            return ref[resid + 1]

        def _get_filtered_out_lines(out_file):
            with open(out_file, "rt") as file:
                lines = file.readlines()
                return [
                    [i for i in line.split(" ") if i != ""]
                    for line in lines
                    if line.split(" ")[0].strip(" ") == "2016Menv"
                ]
        if out_path.is_file():
            flags = [
                P.structure_id,
                out_path.__str__(),
            ]
            handle_debug(print, "out file found for {}".format(P.name))
            hydrophobic_count = 0
            layer_ref = {}
            cen6_ref = {}
            content = _get_filtered_out_lines(
                Path(out_path)
            )
            for line in content:
                res = line[2].strip(" ")
                id = int(line[1].strip(" "))
                layer = int(line[4].strip(" "))
                cen6 = float(line[5].strip(" "))
                layer_ref[id] = layer
                cen6_ref[id] = cen6
                if layer == 3 or layer == 4:
                    hydrophobic_count += 1

            if hydrophobic_count < min_hydrophobic_residues:
                flags.append("below hydrophobicity baseline")
            if len(P.residues) < residue_baseline:
                flags.append("below residue baseline")
            if len(layer_ref) != len(P.residues):
                flags.append("out file / pdb residue count mismatch")

            if len(flags) > 2:
                return Exception(f"{P.structure_id} out file indicates bad structure with flags {', '.join(flags)}")
            else:
                for res in P.residues:
                    P.residues[res].layerinfo = _get_layer_resid(res, layer_ref)
                    # P.residues[res].tmpcen6info = _get_cen6_resid(res, cen6_ref)  # used to validate Eenv calc alg
            return "SUCCESS"
        else:
            return Exception(f"out file for {P.structure_id} does not exist in {Path(out_path).parent}")

    def _init_protein(self, id, file_path, skip_clique_gen=False, skip_layer_info=True, out_path=None,
                      skip_bfactor_check=False, skip_pLDDT_check=False):
        try:
            P = CentroidProtein(
                id,
                file_path,
                exclude_backbone=self.exclude_backbone,
                distance_cutoff=self.distance_cutoff,
                filter_bfactor=self.filter_bfactor,
                is_alphafold=self.is_alphafold,
                filter_pLDDT=self.filter_pLDDT
            )
        except:
            e = sys.exc_info()[0]
            return Exception(e)
        if len(P.residues) > 0:
            if not skip_layer_info:
                res = self._process_out_file(P, out_path)
                if isinstance(res, Exception):
                    return res
            else:
                for res in P.residues:
                    P.residues[res].layerinfo = 7
                handle_debug(print, "{} skipped layer info merging".format(P.structure_id))

            if not skip_clique_gen:
                P.generate_centroid_cliques(skip_bfactor_check=skip_bfactor_check, skip_pLDDT_check=skip_pLDDT_check)
            else:
                handle_debug(print, "{} skipped clique gen".format(P.structure_id))
        else:
            return Exception("{} is empty".format(P.structure_id))
        return P
