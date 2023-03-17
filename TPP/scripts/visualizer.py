import matplotlib as plt
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import plotly.express as px
import numpy as np
import py3Dmol


def view_pdb(pdb_path):
  with open(pdb_path, "rt") as file:
    view = py3Dmol.view(width=800, height=600)
    view.addModel(file.read(), "pdb")
    style = {'cartoon': {'colorscheme': 'chain'}}
    view.setStyle({'model': -1}, style)
    view.zoomTo()
    view.show()

def view_pdb(pdb_path, protein_obj):
    with open(pdb_path, "rt") as file:
        view = py3Dmol.view(width=800, height=600)

def draw_countplot(data, x, title):
    df = pd.DataFrame()
    df[x] = data
    ax = sns.countplot(x=x, data=df, order=list(range(max(data) + 1)))
    plt.title("Clique Sizes Histogram")
    for p in ax.patches:
        ax.annotate("{}".format(p.get_height()), (p.get_x() + 0.1, p.get_height()))
    plt.show()


def draw_histogram(data, title, normalized=False):
    ax = None
    if not normalized:
        ax = sns.distplot(data, kde=False, norm_hist=False)
    else:
        ax = sns.distplot(data)
    plt.show()


def draw_heatmap(
    name, heatmap_data, x_labels, y_labels, cmap, center=0, path_to_dir=r""
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


# plt.show()


def draw_3d_heatmap():  # TODO: Currently just a template, need to load/display actual data
    x = []
    y = []
    z = []
    res = [
        "GLY",
        "PRO",
        "ASP",
        "GLU",
        "LYS",
        "ARG",
        "HIS",
        "SER",
        "THR",
        "ASN",
        "GLN",
        "ALA",
        "MET",
        "TYR",
        "TRP",
        "VAL",
        "ILE",
        "LEU",
        "PHE",
        "CYS",
    ]
    frame = []
    for k in range(20):
        for i in range(20):
            for j in range(20):
                x.append(i + 1)
                y.append(j + 1)
                z.append(k)
                frame.append(res[k])
    heat = np.random.randint(100, size=(20, 20, 20)).flatten()
    df = pd.DataFrame({"x": x, "y": y, "z": z, "frame": frame, "heat": heat})

    fig = px.scatter_3d(
        df,
        x="x",
        y="y",
        z="z",
        animation_frame="frame",
        color="heat",
        range_x=[0, 21],
        range_y=[0, 21],
        range_z=[0, 21],
    )
    fig.update_layout(
        scene=dict(
            xaxis=dict(tickmode="array", ticktext=res, tickvals=list(range(1, 21))),
            yaxis=dict(tickmode="array", ticktext=res, tickvals=list(range(1, 21))),
            zaxis=dict(tickmode="array", ticktext=res, tickvals=list(range(0, 20))),
        )
    )
    fig["layout"].pop("updatemenus")

    fig.show()


# draw_countplot([len(s) for s in centroid_protein.P.centroid_cliques], "clique_sizes", "Clique Size Frequencies")
# draw_histogram(centroid_protein.P.centroid_clique_distances, "Centroid Clique Distances", normalized=True)
