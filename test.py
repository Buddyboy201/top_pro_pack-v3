import sys
sys.path.insert(1, r'C:\Users\aprak\PycharmProjects\TopProPack_v2_2\API')

from API.atom import *
from API.residue import *
from API.centroid_protein import *
from API.energy import *
from API.visualizer import *
from API.top_pro_pack import *

E = TPP_Engine()
E.create_new_project(name="project_1")
#P = E.init_protein("7prc.pdb", "")
with open(r"C:\alpha\7prc.pdb", "rt") as file:
    data = file.read().split("\n")
    P = E.init_protein("7prc.pdb", "", json_load=False, data_load=False, raw_data=data)

