from TPP.scripts.visualizer import draw_heatmap
from TPP.API.energy import EnergyND2
from TPP.db.query import get_rows
from pandas import DataFrame
import os
from pathlib import Path

def generate_heatmap(name, db_path, M, heatmap_dir, layer="ALL"):
    HYDROPHOBIC_DIFF = ["1", "2", "5", "6"]
    INTERFACE_DIFF = ["1", "3", "4", "6"]
    WATER_DIFF = ["2", "3", "4", "5"]

    res, rows, diff = None, list(), None
    if layer == "ALL":
        res = get_rows(db_path, ['clique'], size=M)
        rows = list(map(lambda x: x[0].split(";"), res))
    else:
        if layer == "HYDROPHOBIC":
            diff = HYDROPHOBIC_DIFF
        elif layer == "INTERFACE":
            diff = INTERFACE_DIFF
        elif layer == "WATER":
            diff = WATER_DIFF
        res = get_rows(db_path, ['clique', 'layerinfo'], size=M)
        for r in res:
            layers = set(r[1].split(";"))
            if layers.isdisjoint(diff):
                rows.append(r[0].split(";"))

    E_test = EnergyND2(M=M, cliques=rows)
    E_test.update_epair_table()

    AAs = [
        "G",
        "P",
        "D",
        "E",
        "K",
        "R",
        "H",
        "S",
        "T",
        "N",
        "Q",
        "A",
        "M",
        "Y",
        "W",
        "V",
        "I",
        "L",
        "F",
        "C",
    ]

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
        "CYS": 19,
    }

    if M <= 2:
        draw_heatmap(
            name,
            DataFrame(E_test.STATIC_EPAIR_TABLE),
            AAs,
            AAs,
            "gist_rainbow_r",
            path_to_dir=os.path.abspath(heatmap_dir),
        )
    elif M == 3:
        for i in ref:
            draw_heatmap(
                name + "_{}".format(i),
                DataFrame(E_test.STATIC_EPAIR_TABLE[ref[i]]),
                AAs,
                AAs,
                "gist_rainbow_r",
                path_to_dir=os.path.abspath(heatmap_dir),
            )
    elif M == 4:
        for i in ref:
            for j in ref:
                draw_heatmap(
                    name + "_{}_{}".format(i, j),
                    DataFrame(E_test.STATIC_EPAIR_TABLE[(ref[i], ref[j])]),
                    AAs,
                    AAs,
                    "gist_rainbow_r",
                    path_to_dir=os.path.abspath(heatmap_dir)
                )
    else:
        print("Higher order cliques beyond M=4 not yet supported")

def generate_all_2d_3d_heatmaps(path_to_heatmaps_dir, db_path, date_created, layers=["ALL", "HYDROPHOBIC", "INTERFACE", "WATER"]):
    if not Path(os.path.abspath(path_to_heatmaps_dir)).is_dir():
        Path(os.path.abspath(path_to_heatmaps_dir)).mkdir(parents=True)

    ### ALL_LAYERS ###
    if "ALL" in layers:
        all_layers_2d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "all_layers_plots", "2d")
        )
        all_layers_3d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "all_layers_plots", "3d")
        )
        all_layers_4d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "all_layers_plots", "4d")
        )

        if not all_layers_2d_path.is_dir():
            all_layers_2d_path.mkdir(parents=True)
        if not all_layers_3d_path.is_dir():
            all_layers_3d_path.mkdir(parents=True)
        if not all_layers_4d_path.is_dir():
            all_layers_4d_path.mkdir(parents=True)

        generate_heatmap(
            f"ALL_LAYERS_E_test_M2_{date_created}", db_path, 2, all_layers_2d_path
        )
        generate_heatmap(
            f"ALL_LAYERS_E_test_M3_{date_created}", db_path, 3, all_layers_3d_path
        )
        generate_heatmap(
            f"ALL_LAYERS_E_test_M4_{date_created}", db_path, 4, all_layers_4d_path
        )

    ### HYDROPHOBIC_LAYER ###
    if "HYDROPHOBIC" in layers:
        hydrophobic_2d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "hydrophobic_plots", "2d")
        )
        hydrophobic_3d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "hydrophobic_plots", "3d")
        )
        hydrophobic_4d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "hydrophobic_plots", "4d")
        )

        if not hydrophobic_2d_path.is_dir():
            hydrophobic_2d_path.mkdir(parents=True)
        if not hydrophobic_3d_path.is_dir():
            hydrophobic_3d_path.mkdir(parents=True)
        if not hydrophobic_4d_path.is_dir():
            hydrophobic_4d_path.mkdir(parents=True)

        generate_heatmap(
            f"HYDROPHOBIC_E_test_M2_{date_created}",
            db_path,
            2,
            hydrophobic_2d_path,
            layer="HYDROPHOBIC",
        )
        generate_heatmap(
            f"HYDROPHOBIC_E_test_M3_{date_created}",
            db_path,
            3,
            hydrophobic_3d_path,
            layer="HYDROPHOBIC",
        )
        generate_heatmap(
            f"HYDROPHOBIC_E_test_M4_{date_created}",
            db_path,
            4,
            hydrophobic_4d_path,
            layer="HYDROPHOBIC",
        )

    ### INTERFACE_LAYER ###
    if "INTERFACE" in layers:
        interface_2d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "interface_plots", "2d")
        )
        interface_3d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "interface_plots", "3d")
        )
        interface_4d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "interface_plots", "4d")
        )

        if not interface_2d_path.is_dir():
            interface_2d_path.mkdir(parents=True)
        if not interface_3d_path.is_dir():
            interface_3d_path.mkdir(parents=True)
        if not interface_4d_path.is_dir():
            interface_4d_path.mkdir(parents=True)

        generate_heatmap(
            f"INTERFACE_E_test_M2_{date_created}",
            db_path,
            2,
            interface_2d_path,
            layer="INTERFACE",
        )
        generate_heatmap(
            f"INTERFACE_E_test_M3_{date_created}",
            db_path,
            3,
            interface_3d_path,
            layer="INTERFACE",
        )
        generate_heatmap(
            f"INTERFACE_E_test_M4_{date_created}",
            db_path,
            4,
            interface_4d_path,
            layer="INTERFACE",
        )

    ### WATER_LAYER ###
    if "WATER" in layers:
        water_2d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "water_plots", "2d")
        )
        water_3d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "water_plots", "3d")
        )
        water_4d_path = Path(
            os.path.join(os.path.abspath(path_to_heatmaps_dir), "water_plots", "4d")
        )

        if not water_2d_path.is_dir():
            water_2d_path.mkdir(parents=True)
        if not water_3d_path.is_dir():
            water_3d_path.mkdir(parents=True)
        if not water_4d_path.is_dir():
            water_4d_path.mkdir(parents=True)

        generate_heatmap(
            f"WATER_E_test_M2_{date_created}", db_path, 2, water_2d_path, layer="WATER"
        )
        generate_heatmap(
            f"WATER_E_test_M3_{date_created}", db_path, 3, water_3d_path, layer="WATER"
        )
        generate_heatmap(
            f"WATER_E_test_M4_{date_created}", db_path, 4, water_4d_path, layer="WATER"
        )



