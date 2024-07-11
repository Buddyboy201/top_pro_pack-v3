from bs4 import BeautifulSoup
import numpy as np
from decimal import Decimal

# TODO: create separate paths for tmaf/tmdet and std/out cases to apply tmaf coord adjustment
# TODO: use thresholds Vladimir suggested for inferring tmaf after coord adjustment mem layer info
# TODO: set up test to compare results of tmaf mem layer analysis vs golden set
# TODO: cmdline tool that runs whole pipeline without code
# TODO: update docs (both in code and in README/collab)

# need to create a way to adjust the coordinates of every atom such that membrane normal plane is flat halfway
# across structure z-axis


# layer cutoffs for tmaf layer analysis
# Hydrophobic: +/- 11 A
# Interface: from -11 to -22 and +11 to +22
# Water: from -22 to -30 and +22 to +30

# 1) parse rotation matrix and translation vector
# 2) update all atom coords with rotation and translation
# 3) all layers have fixed ranges from the normal, and some atoms lie in between these ranges
# we want to make it so ranges are as consistent as possible with Vladimir's version of this approach
# limitation: this only works for AF predicted structures
# we need a dataset to compare my results against golden (Vladimir's results for similar dataset)


# 1)

def get_transformation(tmdet_file):
    with open(tmdet_file, "rt") as tmdet_file:
        R = np.empty(9).reshape(3, 3)
        T = np.empty(3).reshape(3, 1)
        tmdet = BeautifulSoup(tmdet_file, "xml")
        x_row = tmdet.find("ROWX")
        y_row = tmdet.find("ROWY")
        z_row = tmdet.find("ROWZ")
        R[0, 0] = x_row.get("X")
        R[0, 1] = x_row.get("Y")
        R[0, 2] = x_row.get("Z")
        R[1, 0] = y_row.get("X")
        R[1, 1] = y_row.get("Y")
        R[1, 2] = y_row.get("Z")
        R[2, 0] = z_row.get("X")
        R[2, 1] = z_row.get("Y")
        R[2, 2] = z_row.get("Z")
        T[0] = x_row.get("T")
        T[1] = y_row.get("T")
        T[2] = z_row.get("T")
        return lambda v: R.dot(v) + T

def get_layer(z):
    if z >= 0 and z < 11:
        return 3
    elif z < 0 and z > -11:
        return 4
    elif z >= 11 and z < 22:
        return 2
    elif z <= -11 and z > -22:
        return 5
    elif z >= 22 and z <= 30:
        return 1
    elif z <= -22 and z >= -30:
        return 6
    else:
        return -1

def gen_pdb(pdb_file, tmdet_file, new_pdb_file, out_file):
    T = get_transformation(tmdet_file)

    with open(out_file, "wt") as out:
        with open(new_pdb_file, "wt") as new_pdb:
            with open(pdb_file) as pdb:
                for line in pdb:
                    if line[0:4] == "ATOM":
                        res_name = line[17:20].strip(" ")
                        res_id = int(line[22:26].strip(" "))
                        coordx = float(line[30:38].strip(" "))
                        coordy = float(line[38:46].strip(" "))
                        coordz = float(line[46:54].strip(" "))
                        T_coords = np.round(T(np.array([coordx, coordy, coordz]).reshape(3, 1)), 3)
                        x, y, z = T_coords[0][0], T_coords[1][0], T_coords[2][0]
                        #print(f"({coordx}, {coordy}, {coordz}) -> ({x}, {y}, {z})") # len of 8
                        x_str = " "*(8-len(f"{x:.3f}"))+f"{x:.3f}"
                        y_str = " " * (8 - len(f"{y:.3f}")) + f"{y:.3f}"
                        z_str = " " * (8 - len(f"{z:.3f}")) + f"{z:.3f}"
                        new_pdb.write(line[0:30]+f"{x_str}{y_str}{z_str}"+line[54:])
                        layer = get_layer(z)
                        out.write(f"2016Menv    {'' if res_id > 9 else ' '}{res_id} {res_name} 15   {layer}    3.0    8.0   11.1   52.5   -0.8   -2.0    3.1    3.4   -0.7   -0.3\n")
                    else:
                        new_pdb.write(line)




