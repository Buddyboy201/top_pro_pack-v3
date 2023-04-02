import csv
from TPP.API.verbose import handle_debug
from TPP.API import verbose
from pathlib import Path
from bs4 import BeautifulSoup


def _get_filtered_out_lines(out_file):
    with open(out_file, "rt") as file:
        lines = file.readlines()
        return [
            [i for i in line.split(" ") if i != ""]
            for line in lines
            if line.split(" ")[0].strip(" ") == "2016Menv"
        ]

def _get_tmdet_tm_count(tmdet_file):
    with open(tmdet_file, "rt") as tmdet_file:
        region_data = BeautifulSoup(tmdet_file, "xml").find_all("REGION", {'type': ['m', 'M']})
        tm_count = 0
        for region in region_data:
            tm_count += int(region.get("pdb_end")) - int(region.get("pdb_beg")) + 1
        return tm_count

def _get_layer_resid(resid, layer_ref):
    return layer_ref[resid + 1]


def _get_clique_with_names_only(sorted_clique):
    return ";".join([res.name for res in sorted_clique])


def _get_clique_with_resid_only(sorted_clique):
    return ";".join([str(i.resid) for i in sorted_clique])


def _get_clique_with_old_resid_only(sorted_clique):
    return ";".join([str(i.old_resid) for i in sorted_clique])


def _get_clique_layer_info_only(sorted_clique, layer_ref=None):
    if layer_ref is not None:
        resids = [res.resid for res in sorted_clique]
        return ";".join([str(_get_layer_resid(resid, layer_ref)) for resid in resids])
    else:
        return "NULL"


def _push_clique_to_buffer(clique, pdb_name, layer_ref, row_counter, buffer):
    sorted_clique = sorted([res for res in clique], key=lambda r: r.name)
    buffer.append(
        {
            "id": row_counter,
            "size": len(sorted_clique),
            "clique": _get_clique_with_names_only(sorted_clique),
            "resid": _get_clique_with_resid_only(sorted_clique),
            "oldresid": _get_clique_with_old_resid_only(sorted_clique),
            "layerinfo": _get_clique_layer_info_only(sorted_clique, layer_ref),
            "pdbname": pdb_name,
        }
    )
    return row_counter + 1

def gen_clique_db_TMAF(proj, db_path, tmdet_dir, min_tm_residues=34, residue_baseline=30, bad_proteins_file_path="bad_proteins_file.txt"):
    with open(db_path, 'wt', newline='') as db_file:
        headers = ["id", "size", "clique", "resid", "oldresid", "layerinfo", "pdbname"]
        row_counter = 0
        writer = csv.DictWriter(db_file, fieldnames=headers, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        bad_proteins = []
        buffer = []
        for pdb_id_clean in proj.proteins:
            # check for tmdet_file
            if Path(tmdet_dir / Path(f"{pdb_id_clean}_tmdet.xml")).is_file():
                # init flags
                flags = [
                    pdb_id_clean,
                    Path(tmdet_dir / Path(f"{pdb_id_clean}_tmdet.xml")).__str__(),
                ]
                handle_debug(print, "tmdet file found for {}".format(pdb_id_clean))
                P = proj.get_protein(pdb_id_clean)
                hydrophobic_count = _get_tmdet_tm_count(Path(tmdet_dir / Path(f"{pdb_id_clean}_tmdet.xml")))
                layer_ref = None
                if hydrophobic_count < min_tm_residues:
                    flags.append("below tm residue baseline")
                if len(P.residues) < residue_baseline:
                    flags.append("below residue baseline")

                if len(flags) > 2:
                    bad_proteins.append(",".join(flags) + "\n")
                else:
                    # handle_debug(print, None, len(P.residues), pdb_id_clean)
                    cliques = P.centroid_cliques
                    # buffer = []
                    for clique in cliques:
                        row_counter = _push_clique_to_buffer(clique, P.name, layer_ref, row_counter, buffer)
            else:
                handle_debug(
                    print,
                    "tmdet file for {} does not exist in {}".format(pdb_id_clean, tmdet_dir),
                )
                bad_proteins.append(",".join((pdb_id_clean, "", "missing tmdet file")) + "\n")
        if verbose.VERBOSE:
            with open(bad_proteins_file_path, "wt") as bp_file:
                bp_file.writelines(bad_proteins)

        writer.writerows(buffer)


def gen_clique_db(proj, db_path, out_dir, min_hydrophobic_residues=34, residue_baseline=30, bad_proteins_file_path="bad_proteins_file.txt"):
    with open(db_path, 'wt', newline='') as db_file:
        headers = ["id", "size", "clique", "resid", "oldresid", "layerinfo", "pdbname"]
        row_counter = 0
        writer = csv.DictWriter(db_file, fieldnames=headers, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        bad_proteins = []
        buffer = []
        for pdb_id_clean in proj.proteins:
            # pdb_id_clean = pdb_id[len("Menv_color_memb_cen_nor_"):len("Menv_color_memb_cen_nor_")+4]
            if Path(out_dir / Path("{}.out".format(pdb_id_clean))).is_file():
                flags = [
                    pdb_id_clean,
                    Path(out_dir / Path("{}.out".format(pdb_id_clean))).__str__(),
                ]
                handle_debug(print, "out file found for {}".format(pdb_id_clean))
                P = proj.get_protein(pdb_id_clean)
                hydrophobic_count = 0
                layer_ref = {}
                content = _get_filtered_out_lines(
                    Path(out_dir / Path("{}.out".format(pdb_id_clean)))
                )
                for line in content:
                    res = line[2].strip(" ")
                    id = int(line[1].strip(" "))
                    layer = int(line[4].strip(" "))
                    layer_ref[id] = layer
                    if layer == 3 or layer == 4:
                        hydrophobic_count += 1

                if hydrophobic_count < min_hydrophobic_residues:
                    flags.append("below hydrophobicity baseline")
                if len(P.residues) < residue_baseline:
                    flags.append("below residue baseline")
                if len(layer_ref) != len(P.residues):
                    flags.append("out file / pdb residue count mismatch")

                if len(flags) > 2:
                    bad_proteins.append(",".join(flags) + "\n")
                else:
                    # handle_debug(print, len(layer_ref), len(P.residues), pdb_id_clean)
                    cliques = P.centroid_cliques
                    # buffer = []
                    for clique in cliques:
                        row_counter = _push_clique_to_buffer(clique, P.name, layer_ref, row_counter, buffer)
            else:
                handle_debug(
                    print,
                    "out file for {} does not exist in {}".format(pdb_id_clean, out_dir),
                )
                bad_proteins.append(",".join((pdb_id_clean, "", "missing out file")) + "\n")

        if verbose.VERBOSE:
            with open(bad_proteins_file_path, "wt") as bp_file:
                bp_file.writelines(bad_proteins)

        writer.writerows(buffer)
