import sys
sys.path.insert(1, r'C:\Users\aprak\PycharmProjects\TopProPack_v2_2\API')

from boxsdk import JWTAuth, Client, OAuth2
from API.atom import *
from API.residue import *
from API.centroid_protein import *
from API.energy import *
from API.visualizer import *
from API.top_pro_pack import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from API.visualizer import *

# .out -> res => layer mapping
# .out -> .pdb => cliques
# (cliques + (res => mapping)) => cliques filtered by layer
# manually build pairwise counts from cliques and P.E.ref
# update epair table for each protein using associated pairwise counts

# step 1: .out -> res => layer mapping
auth = JWTAuth.from_settings_file('config.json')
access_token = auth.authenticate_instance()
client_auth = OAuth2(client_id="ukaag2sjuuko9fjei5n9q6fufg4ca6b2",
                     client_secret="a87GrrszqV6LthoU8e2CwkC8Tec9DYQN",
                     access_token=access_token)
client = Client(client_auth)
current_user = client.user(user_id='me').get()
#print(current_user.name)
app_user = client.user(user_id="8678622593")
app_client = client.as_user(app_user)
items_iter = app_client.folder(folder_id="115006420152").get_items()


items = []
while True:
    try:
        file = items_iter.next()
        if ".out" in file.name and "id" not in file.name: items.append(file)
    except: break

pdb_iter = app_client.folder(folder_id="115006398552").get_items()

pdb_items = []
while True:
    try:
        file = pdb_iter.next()
        if ".pdb" in file.name:
            pdb_items.append(file)
    except: break


for i in range(len(items)):
    items[i].name = items[i].name[:4] + ".out"
    #print(items[i].name)

for i in range(len(pdb_items)):
    pdb_items[i].name = pdb_items[i].name[24:28] + ".pdb"
    #print(pdb_items[i].name)

items.sort(key=lambda item: item.name)
pdb_items.sort(key=lambda item: item.name)

#print(len(items), len(pdb_items))

missing_protein = "6idp"
bad_pos = 390
#for i in range(len(pdb_items)):
#    if pdb_items[i].name[:4] == missing_protein:
        #pdb_items.pop(i)
#        bad_pos = i
#print(bad_pos)
pdb_items.pop(bad_pos)

#for out, pdb in list(zip(items, pdb_items)):
#    if out.name[:4] != pdb.name[:4]:
#        print(out.name, pdb.name)

def get_filtered_lines(file_id, app_client):
    content = app_client.file(file_id).content().decode("utf-8").split("\n")
    return [[i for i in line.split(" ") if i != ""] for line in content if line.split(" ")[0].strip(" ") == "2016Menv"]

def get_pdb_lines(file_id, app_client):
    return app_client.file(file_id).content().decode("utf-8").split("\n")

#print(get_pdb_lines(pdb_items[0].id, app_client))

min_hydrophobic_residues = 14

water_cliques = []
interface_cliques = []
hydrophobic_cliques = []

def get_layer(res, layer_ref):
    return layer_ref[res.resid+1]

def filter_water(cliques, layer_ref):
    new_cliques = []
    for clique in cliques:
        if len([layer for layer in [get_layer(res, layer_ref) for res in clique] if layer not in [1, 6]]) == 0:
            new_cliques.append(clique)
    return new_cliques


def filter_interface(cliques, layer_ref):
    new_cliques = []
    for clique in cliques:
        if len([layer for layer in [get_layer(res, layer_ref) for res in clique] if layer not in [2, 5]]) == 0:
            new_cliques.append(clique)
    return new_cliques

def filter_hydrophobic(cliques, layer_ref):
    new_cliques = []
    for clique in cliques:
        if len([layer for layer in [get_layer(res, layer_ref) for res in clique] if layer not in [3, 4]]) == 0:
            new_cliques.append(clique)
    return new_cliques


good_proteins = []
bad_proteins = []
exceptions = []
E_total = Energy()
E_water = Energy()
E_interface = Energy()
E_hydrophobic = Energy()
Engine = TPP_Engine()
log_num = 0
file = open("layer_heatmap_log_{}.txt".format(log_num), "w")
for out, pdb in list(zip(items, pdb_items)):
    pdb_name = pdb.name
    P = Engine.init_protein(pdb_name, "", json_load=False, data_load=False, raw_data=get_pdb_lines(pdb.id, app_client))
    content = get_filtered_lines(out.id, app_client)
    hydrophobic_count = 0
    layer_ref = {}
    for line in content:
        res = line[2].strip(" ")
        id = int(line[1].strip(" "))
        layer = int(line[4].strip(" "))
        layer_ref[id] = layer
        if layer == 3 or layer == 4:
            hydrophobic_count += 1
    if type(P) is Exception:
        exceptions.append(pdb_name + " " + str(P))
    elif hydrophobic_count >= min_hydrophobic_residues:
        good_proteins.append(P)
        E_total.update_static_total_pairs_table(get_protein_pairs_matrix(P.centroid_cliques))
        E_water.update_static_total_pairs_table(get_protein_pairs_matrix(filter_water(P.centroid_cliques, layer_ref)))
        E_interface.update_static_total_pairs_table(get_protein_pairs_matrix(filter_interface(P.centroid_cliques, layer_ref)))
        E_hydrophobic.update_static_total_pairs_table(get_protein_pairs_matrix(filter_hydrophobic(P.centroid_cliques, layer_ref)))
        file.write("{} {} {}\n".format(pdb_name, len(content), len(P.residues)))
    else:
        bad_proteins.append(P)

file.write("bad_proteins:\n")
for i in bad_proteins:
    file.write("\t"+i.name+"\n")
file.write("exceptions:\n")
for i in exceptions:
    file.write("\t"+i+"\n")
file.close()
E_total.update_epair_values()
E_water.update_epair_values()
E_interface.update_epair_values()
E_hydrophobic.update_epair_values()
print(E_total.get_static_epair_table())
print(E_water.get_static_epair_table())
print(E_interface.get_static_epair_table())
print(E_hydrophobic.get_static_epair_table())
labels = "G P D E K R H S T N Q A M Y W V I L F C".split(" ")
draw_heatmap("E_total", E_total.get_static_epair_table(), labels, labels, "gist_rainbow_r")
draw_heatmap("E_water", E_water.get_static_epair_table(), labels, labels, "gist_rainbow_r")
draw_heatmap("E_interface", E_interface.get_static_epair_table(), labels, labels, "gist_rainbow_r")
draw_heatmap("E_hydrophobic", E_hydrophobic.get_static_epair_table(), labels, labels, "gist_rainbow_r")


test = '''
good_proteins = []
bad_proteins = []
E = Energy()
Engine = TPP_Engine()
short, long = 0, 0
for i in items[:30]:
    pdb_name = i.name[:4] + ".pdb"
    out = Engine.init_protein(pdb_name, "")
    content = get_filtered_lines(i.id, app_client)
    hydrophobic_count = 0
    layer_ref = {}
    for line in content:
        res = line[2].strip(" ")
        id = int(line[1].strip(" "))
        #print(res, out.residues[id-1].name)
        layer = int(line[4].strip(" "))
        if layer == 3 or layer == 4:
            hydrophobic_count += 1
    if type(out) is Exception: print(out)
    elif hydrophobic_count >= min_hydrophobic_residues:
        good_proteins.append(out)
    else:
        bad_proteins.append(out)
    print(len(content), len(out.residues))
    if len(content) > len(out.residues): long += 1
    elif len(content) < len(out.residues): short += 1
print(len(items), len(good_proteins), len(bad_proteins))
#print(short, long)

'''

