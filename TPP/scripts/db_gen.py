from TPP.API.top_pro_pack import Project
from pathlib import Path
from sqlalchemy import create_engine, MetaData, Table, Column, Integer, String
from TPP.API.verbose import handle_debug
from TPP.API import verbose

def gen_table_and_connect_engine(db_path):
    engine = create_engine(db_path, echo=True)
    meta = MetaData()
    cliques_table = Table(
        "cliques", meta,
        Column("id", Integer, primary_key=True),
        Column("size", Integer),
        Column("clique", String),
        Column("resid", String),
        Column("oldresid", String),
        Column("layerinfo", String),
        Column("pdbname", String)
    )
    meta.create_all(engine)
    conn = engine.connect()
    return conn, cliques_table, meta, engine


def _get_filtered_out_lines(out_file):
    with open(out_file, "rt") as file:
        lines = file.readlines()
        return [[i for i in line.split(" ") if i != ""]
                for line in lines if line.split(" ")[0].strip(" ") == "2016Menv"]


def _get_layer_resid(resid, layer_ref):
    return layer_ref[resid + 1]


def _get_clique_with_names_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([i.name for i in clique])


def _get_clique_with_resid_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([str(i.resid) for i in clique])


def _get_clique_with_old_resid_only(clique):
    clique.sort(key=lambda x: x.name)
    return ";".join([str(i.old_resid) for i in clique])


def _get_clique_layer_info_only(clique, layer_ref):
    clique.sort(key=lambda x: x.name)
    resids = [res.resid for res in clique]
    return ";".join([str(_get_layer_resid(resid, layer_ref)) for resid in resids])


def _push_clique_to_buffer(clique, pdb_name, layer_ref, buffer):
    buffer.append({
        "size": len(clique), "clique": _get_clique_with_names_only(clique),
        "resid": _get_clique_with_resid_only(clique),
        "oldresid": _get_clique_with_old_resid_only(clique),
        "layerinfo": _get_clique_layer_info_only(clique, layer_ref), "pdbname": pdb_name
        }
    )


def _bulk_insert_cliques_into_db(buffer, conn, table):
    return conn.execute(table.insert(), buffer)

# DEPRECATED - bulk insert should be used instead for faster db updates
def _insert_clique_into_db(clique, pdb_name, layer_ref, conn, table):
    ins = table.insert().values(size=len(clique), clique=_get_clique_with_names_only(clique),
                                resid=_get_clique_with_resid_only(clique),
                                oldresid=_get_clique_with_old_resid_only(clique),
                                layerinfo=_get_clique_layer_info_only(clique, layer_ref), pdbname=pdb_name)
    result = conn.execute(ins)
    return result


'''def generate_clique_db(proj, cliques_table, conn, out_dir, min_hydrophobic_residues=34, residue_baseline=30):
    bad_proteins = []
    for pdb_id in proj.proteins:
        pdb_id_clean = pdb_id[len("Menv_color_memb_cen_nor_"):len("Menv_color_memb_cen_nor_")+4]
        if Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))).is_file():
            flags = [pdb_id, Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))).__str__()]
            handle_debug(print, "out file found for {}".format(pdb_id))
            P = proj.get_protein(pdb_id)
            hydrophobic_count = 0
            layer_ref = {}
            content = _get_filtered_out_lines(Path(out_dir / Path("{}_Menv.out".format(pdb_id_clean))))
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
                handle_debug(print, len(layer_ref), len(P.residues), pdb_id_clean)
                cliques = P.centroid_cliques
                buffer = []
                for clique in cliques:
                    # insert_clique_into_db(clique, P.name, layer_ref, conn, cliques_table)
                    _push_clique_to_buffer(clique, P.name, layer_ref, buffer)
                _bulk_insert_cliques_into_db(buffer, conn, cliques_table)

            # insert comment 3 here

        else:
            handle_debug(print, "out file for {} does not exist in {}".format(pdb_id, out_dir))
            bad_proteins.append(",".join((pdb_id, "", "missing out file")) + "\n")

    if verbose.VERBOSE:
        with open("bad_proteins_file.txt", "wt") as bp_file:
            bp_file.writelines(bad_proteins)'''



def generate_clique_db(proj, cliques_table, conn, out_dir, min_hydrophobic_residues=34, residue_baseline=30):
    bad_proteins = []
    buffer = []
    for pdb_id_clean in proj.proteins:
        # pdb_id_clean = pdb_id[len("Menv_color_memb_cen_nor_"):len("Menv_color_memb_cen_nor_")+4]
        if Path(out_dir / Path("{}.out".format(pdb_id_clean))).is_file():
            flags = [pdb_id_clean, Path(out_dir / Path("{}.out".format(pdb_id_clean))).__str__()]
            handle_debug(print, "out file found for {}".format(pdb_id_clean))
            P = proj.get_protein(pdb_id_clean)
            hydrophobic_count = 0
            layer_ref = {}
            content = _get_filtered_out_lines(Path(out_dir / Path("{}.out".format(pdb_id_clean))))
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
                handle_debug(print, len(layer_ref), len(P.residues), pdb_id_clean)
                cliques = P.centroid_cliques
                # buffer = []
                for clique in cliques:
                    # insert_clique_into_db(clique, P.name, layer_ref, conn, cliques_table)
                    _push_clique_to_buffer(clique, P.name, layer_ref, buffer)
                # _bulk_insert_cliques_into_db(buffer, conn, cliques_table)

            # insert comment 3 here

        else:
            handle_debug(print, "out file for {} does not exist in {}".format(pdb_id_clean, out_dir))
            bad_proteins.append(",".join((pdb_id_clean, "", "missing out file")) + "\n")

    if verbose.VERBOSE:
        with open("bad_proteins_file.txt", "wt") as bp_file:
            bp_file.writelines(bad_proteins)

    _bulk_insert_cliques_into_db(buffer, conn, cliques_table)


