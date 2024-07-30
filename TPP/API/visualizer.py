import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from TPP.API.energy import EnergyMD
from pandas import DataFrame
import os
from pathlib import Path
from TPP.API.constants import AA_MINIs, AA_REF


def draw_heatmap(
    name, heatmap_data, x_labels, y_labels, cmap, path_to_dir=r""
):
    plot = sns.heatmap(
        heatmap_data,
        xticklabels=x_labels,
        yticklabels=y_labels,
        center=0,
        vmin=-1,
        vmax=1,
        robust=True,
        cmap=cmap,
    )
    plt.savefig(Path(Path(path_to_dir) / Path("{}.png".format(name))))
    plt.clf()


def generate_heatmap(name, project, heatmap_dir, M=2, L="ALL"):
    E_group = EnergyMD(project, M=M, L=L)
    if M <= 2:
        draw_heatmap(
            name,
            DataFrame(E_group.get_result()),
            AA_MINIs,
            AA_MINIs,
            "gist_rainbow_r",
            path_to_dir=os.path.abspath(heatmap_dir),
        )
    elif M == 3:
        for i in AA_REF:
            draw_heatmap(
                name + "_{}".format(i),
                DataFrame(E_group.get_result()[AA_REF[i]]),
                AA_MINIs,
                AA_MINIs,
                "gist_rainbow_r",
                path_to_dir=os.path.abspath(heatmap_dir),
            )
    elif M == 4:
        for i in AA_REF:
            for j in AA_REF:
                draw_heatmap(
                    name + "_{}_{}".format(i, j),
                    DataFrame(E_group.get_result()[(AA_REF[i], AA_REF[j])]),
                    AA_MINIs,
                    AA_MINIs,
                    "gist_rainbow_r",
                    path_to_dir=os.path.abspath(heatmap_dir)
                )
    else:
        print("Higher order cliques beyond M=4 not yet supported")


def generate_all_2d_3d_heatmaps(project, path_to_heatmaps_dir, date_created, layers=("ALL", "HYDROPHOBIC", "INTERFACE", "WATER")):
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
            f"ALL_LAYERS_E_group_M2_{date_created}", project, all_layers_2d_path, M=2, L="ALL"
        )
        generate_heatmap(
            f"ALL_LAYERS_E_group_M3_{date_created}", project, all_layers_3d_path, M=3, L="ALL"
        )
        generate_heatmap(
            f"ALL_LAYERS_E_group_M4_{date_created}", project, all_layers_4d_path, M=4, L="ALL"
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
            f"HYDROPHOBIC_E_group_M2_{date_created}", project, hydrophobic_2d_path, M=2, L="HYDROPHOBIC"
        )
        generate_heatmap(
            f"HYDROPHOBIC_E_group_M3_{date_created}", project, hydrophobic_3d_path, M=3, L="HYDROPHOBIC"
        )
        generate_heatmap(
            f"HYDROPHOBIC_E_group_M4_{date_created}", project, hydrophobic_4d_path, M=4, L="HYDROPHOBIC"
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
            f"INTERFACE_E_group_M2_{date_created}", project, interface_2d_path, M=2, L="INTERFACE"
        )
        generate_heatmap(
            f"INTERFACE_E_group_M3_{date_created}", project, interface_3d_path, M=3, L="INTERFACE"
        )
        generate_heatmap(
            f"INTERFACE_E_group_M4_{date_created}", project, interface_4d_path, M=4, L="INTERFACE"
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
            f"WATER_E_group_M2_{date_created}", project, water_2d_path, M=2, L="WATER"
        )
        generate_heatmap(
            f"WATER_E_group_M3_{date_created}", project, water_3d_path, M=3, L="WATER"
        )
        generate_heatmap(
            f"WATER_E_group_M4_{date_created}", project, water_4d_path, M=4, L="WATER"
        )


def plot_E_env(aa_i, result_res, plots_dir, layers=(1, 2, 3, 4, 5, 6)):
    mpl.style.use("seaborn")
    fig, ax = plt.subplots()
    L_color_map = {
        "W_OUT": "deepskyblue",
        "I_OUT": "lightgreen",
        "H_OUT": "lightsalmon",
        "H_IN": "red",
        "I_IN": "green",
        "W_IN": "blue",
        "ALL": "purple"
    }
    layer_id_L_map = {
        1: "W_IN",
        2: "I_IN",
        3: "H_IN",
        4: "H_OUT",
        5: "I_OUT",
        6: "W_OUT",
        7: "ALL"
    }
    for layer in layers:
        items = sorted(result_res[layer].items())
        sorted_keys, sorted_values = zip(*items)
        ll = layer_id_L_map[layer]
        ax.scatter(sorted_keys, sorted_values, color=L_color_map[ll])
        ax.plot(sorted_keys, sorted_values, color=L_color_map[ll], label=ll)
    ax.set_title(f"{aa_i}_E_env_plot")
    ax.set_xlabel("burial count")
    ax.set_ylabel('E_env')
    ax.set_ylim(-2, 5)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.90, box.height])
    # Put a legend to the right of the current axis
    # ax.legend()
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(Path(plots_dir) / Path(f"{aa_i}_E_env_plot.png"))
    plt.clf()

