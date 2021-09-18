from TPP.scripts.visualizer import draw_heatmap
from TPP.API.energy import EnergyND2
from pandas import DataFrame
import os
from pathlib import Path

def generate_heatmap(name, M, conn, heatmap_dir, layer="ALL"):
    stmt = None
    HYDROPHOBIC_DIFF = ["1", "2", "5", "6"]
    INTERFACE_DIFF = ["1", "3", "4", "6"]
    WATER_DIFF = ["2", "3", "4", "5"]
    diff = None
    if layer == "ALL":
        stmt = "SELECT clique FROM CLIQUES WHERE size={}".format(M)
    else:
        stmt = "SELECT clique, layerinfo FROM CLIQUES WHERE size={}".format(M)
        if layer == "HYDROPHOBIC": diff = HYDROPHOBIC_DIFF
        elif layer == "INTERFACE": diff = INTERFACE_DIFF
        elif layer == "WATER": diff = WATER_DIFF


    res = conn.execute(stmt)
    rows = None
    if layer != "ALL":
        rows = []
        for r in res:
            layers = set(r[1].split(";"))
            if layers.isdisjoint(diff):
                rows.append(r[0].split(";"))
    else:
        rows = list(map(lambda x: x[0].split(";"), list(res)))

    E_test = EnergyND2(M=M, cliques=rows)
    E_test.update_epair_table()

    AAs = ["G", "P", "D", "E", "K", "R", "H", "S", "T", "N", "Q", "A", "M", "Y", "W", "V", "I", "L", "F", "C"]

    ref = {
        "GLY": 0,
        "PRO": 1,
        "ASP": 2,
        "GLU": 3,
        "LYS": 4,
        "ARG": 5,
        "HIS": 6,
        "SER": 7,
        "THR": 8,
        "ASN": 9,
        "GLN": 10,
        "ALA": 11,
        "MET": 12,
        "TYR": 13,
        "TRP": 14,
        "VAL": 15,
        "ILE": 16,
        "LEU": 17,
        "PHE": 18,
        "CYS": 19
    }

    if M <= 2:
        draw_heatmap(name, DataFrame(E_test.STATIC_EPAIR_TABLE), AAs, AAs, "gist_rainbow_r", path_to_dir=os.path.abspath(heatmap_dir))
    elif M == 3:
        for i in ref:
            draw_heatmap(name + "_{}".format(i), DataFrame(E_test.STATIC_EPAIR_TABLE[ref[i]]), AAs, AAs, "gist_rainbow_r", path_to_dir=os.path.abspath(heatmap_dir))
    else:
        print("Higher order cliques beyond M=3 not yet supported")


def generate_all_2d_3d_heatmaps(path_to_heatmaps_dir, conn, date_created):

    ### ALL_LAYERS ###

    all_layers_2d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "all_layers_plots", "2d"))
    all_layers_3d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "all_layers_plots", "3d"))

    if not all_layers_2d_path.is_dir():
        all_layers_2d_path.mkdir(parents=True)
    if not all_layers_3d_path.is_dir():
        all_layers_3d_path.mkdir(parents=True)

    generate_heatmap(f"ALL_LAYERS_E_test_M2_{date_created}", 2, conn, all_layers_2d_path)
    generate_heatmap(f"ALL_LAYERS_E_test_M3_{date_created}", 3, conn, all_layers_3d_path)


    ### HYDROPHOBIC_LAYER ###

    hydrophobic_2d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "hydrophobic_plots", "2d"))
    hydrophobic_3d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "hydrophobic_plots", "3d"))

    if not hydrophobic_2d_path.is_dir():
        hydrophobic_2d_path.mkdir(parents=True)
    if not hydrophobic_3d_path.is_dir():
        hydrophobic_3d_path.mkdir(parents=True)

    generate_heatmap(f'HYDROPHOBIC_E_test_M2_{date_created}', 2, conn, hydrophobic_2d_path, layer="HYDROPHOBIC")
    generate_heatmap(f'HYDROPHOBIC_E_test_M3_{date_created}', 3, conn, hydrophobic_3d_path, layer="HYDROPHOBIC")


    ### INTERFACE_LAYER ###

    interface_2d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "interface_plots", "2d"))
    interface_3d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "interface_plots", "3d"))

    if not interface_2d_path.is_dir():
        interface_2d_path.mkdir(parents=True)
    if not interface_3d_path.is_dir():
        interface_3d_path.mkdir(parents=True)

    generate_heatmap(f'INTERFACE_E_test_M2_{date_created}', 2, conn, interface_2d_path, layer="INTERFACE")
    generate_heatmap(f'INTERFACE_E_test_M3_{date_created}', 3, conn, interface_3d_path, layer="INTERFACE")


    ### WATER_LAYER ###

    water_2d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "water_plots", "2d"))
    water_3d_path = Path(os.path.join(os.path.abspath(path_to_heatmaps_dir), "water_plots", "3d"))

    if not water_2d_path.is_dir():
        water_2d_path.mkdir(parents=True)
    if not water_3d_path.is_dir():
        water_3d_path.mkdir(parents=True)

    generate_heatmap(f'WATER_E_test_M2_{date_created}', 2, conn, water_2d_path, layer="WATER")
    generate_heatmap(f'WATER_E_test_M3_{date_created}', 3, conn, water_3d_path, layer="WATER")




#generate_all_2d_3d_heatmaps(path_to_heatmaps_dir, conn, "9_1_2021")
