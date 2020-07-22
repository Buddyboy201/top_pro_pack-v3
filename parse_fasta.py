import os
import requests

pdb_list = []

with open("fasta_file.txt", "rt") as file:
    lines = file.readlines()
    lines = lines[1:]
    for line in lines:
        pdb_list.append(line[0:4].lower() + ".pdb")

file_path = r"C:\Users\aprak\PycharmProjects\TopProPack_v2_2\{}"
file = open(file_path.format("list_of_pdbs.txt"), "wt")
file.write("{}\n".format(len(pdb_list)))
for f in pdb_list:
    file.write(f+"\n")
file.close()


print(len(pdb_list))

base_url = "https://files.rcsb.org/download/"
full_url = base_url+pdb_list[0]
data = requests.get(full_url).text
print(data)
print(data.split('\n'))