# top_pro_pack-v3
 
## Getting Started

### 1. Creating Projects

TPP_v3 uses json config files to configure and initialize various projects. Config files will contain a name,
a path to a directory containing the project's pdb files, a path to a directory containing the project's json files,
a list of paths that should be ignored in any of the specified directories, and two variables, *exclude_backbone* and 
*distance_cutoff*, which modify the parameters used when generating centroid-based cliques.

Config files can be easily generated with the *TPP.API.top_pro_pack.create_project* function, as shown in this example:

~~~python
from TPP.API.top_pro_pack import create_project
from pathlib import Path

config_path = Path("<insert_config_path_here>")
create_project(config_path, "<insert_project_name_here>", "<insert_json_dir_path_here>", "<insert_pdb_dir_path_here>")
~~~

Optionally, addition parameters can be added to *create_project* to change the *exclude_backbone*, *distance_cutoff*,
and *ignored_paths* parameters from their default values.

Once *create_project* has been used to generate a config file for a project, it does not need to be called again until
additional projects need to be configured.

To use to config file to initialize a project, the path to the config file can be provided to an instance of the 
*TPP.API.top_pro_pack.Project* class to create a *Project* object with the settings provided in the config file like so:

~~~python
from TPP.API.top_pro_pack import Project
from pathlib import Path

config_path = "<insert_config_file_path_here>"
example_proj = Project(config_path)
~~~

### 2. Loading CentroidProtein objects through the Project class

After a project has been initialized with a config file, protein data saved within either the pdb_dir or json_dir specified
can be used to generate *TPP.API.centroid_protein.CentroidProtein* objects for that data, allowing for a variety of more
complex operations to be run on the protein described once it has been processed. **NOTE: functionality for json protein
data has NOT been implemented in its entirety yet - it should not be used unless you are prepared for various errors**
Pdb files can be processed by a *Project* object like follows:

~~~python
from TPP.API.top_pro_pack import Project
from pathlib import Path

config_path = "<insert_config_file_path_here>"
example_proj = Project(config_path)
example_proj.load_protein("<insert_protein_id_here>", "<insert_pdb_file_name__from_pdb_dir_here (e.g. 'test.pdb')>")
test_Protein = example_proj.get_protein("<insert_protein_id_here>") # this command can be used to receive a CentroidProtein object for the pdb data once it has been loaded
~~~

If a use case arises where the whole dataset of pdb files in the specified pdb_dir needs to be loaded, all pdb_files can be
loaded like so:

~~~python
from TPP.API.top_pro_pack import Project
from pathlib import Path

config_path = "<insert_config_file_path_here>"
example_proj = Project(config_path)

example_proj.load_all_pdbs(id_list) # list of ids for each pdb in pdb_dir in order of the files as listed in pdb_dir
~~~

The method *TPP.API.top_pro_pack.Project.generate_default_ids* can be used to create a list ids based of the file names of
the pdb files in pdb_dir while taking into account any ignored files in the directory. It is reccomended this method be used 
when loading the entire pdb_dir since it eliminates any errors relating to pdb_id misatrribution unless a custom set of ids is desired.
An example of this methods use can be seen as follows:

~~~python
from TPP.API.top_pro_pack import Project
from pathlib import Path

config_path = "<insert_config_file_path_here>"
example_proj = Project(config_path)

example_proj.load_all_pdbs(example_proj.generate_default_ids())
~~~

### 3. Creating a Database for Clique Info

*top_pro_pack-v3* makes use of sqlite databases to effectively store clique data and other associated 
pdb information. Each project can generate a table with the following structure:

~~~python
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
~~~
(Above code uses sqlalchemy library; does not run)

In order to generate useful membrane layer information in relation to clique data, it is necessary to also provide
a text file for every pdb file in a project that associates every residue with one of 6 layers (water, interface, and hydrophobic with
each major layer having an "outer" and "inner" membrane variant). This api makes use of special ".out" files that pair with
and individual pdb file to produce complete clique and layer information about a membrane structure. Note that cliques from structures
not meeting minimum baselines for number of residues overall and hydrophobic residues are excluded from analysis as they quite often fail to
provide useful data. Additionally, cliques with residues crossing between layers are also excluded for simplicity.

Sqlite "clique tables" can be generated in a 2-step-process with the functions
*TPP.scripts.db_gen.gen_table_and_connect_engine* and *TPP.scripts.db_gen.generate_clique_db* as seen below:

~~~python
from TPP.scripts.db_gen import gen_table_and_connect_engine, generate_clique_db
proj = "<...>" # project placeholder
db_path = "<insert sqlite db path here>"
out_dir = "<insert '.out' files directory path here>"
conn, cliques_table, meta, engine = gen_table_and_connect_engine(db_path)
generate_clique_db(proj, cliques_table, conn, out_dir)
~~~

### 3.5 Querying SQL Clique Tables

As a pre-step to generating potentials for cliques, queries must first be made to a clique table to isolate the desired data. This can be done through standard
sqlite queries, python scripts using a sql library to parse the database, a combination of the prior two methods for complex queries (recommended), or any method
that results in a python "2d-array" of cliques.

### 4. Calculating Pairwise / "Groupwise" Energies

Once a set of *TPP.API.centroid_protein.CentroidProtein* objects have been loaded into a Project and a sqlite db has been generated, the 
*TPP.API.energy.EnergyND2* class can be used to compute pairwise or "groupwise" (Formula given later in the section) energies within cliques of a given size 
for the structures in the project. 

~~~python
from TPP.API.energy import EnergyND2
cliques = "<...>" # cliques generated from query on cliques table
M = "<insert INTEGER size of cliques here>"
energy = EnergyND2(M, cliques)
energy.update_epair_table()


# Updated matrix of potentials can now be accessed here:
energy.STATIC_EPAIR_TABLE
~~~

The formula used to derive the "groupwise" potentials for a given residue triplet A-B-C is: *-ln(P(A)*P(BC)/P(ABC))*


### 5. Heatmap Visualizations for Potentials

The TPP API makes it extremely easy to generate heatmaps in bulk for analyzing both pairwise potentials as well as residue triplets (currently only M=2 or 3
supported for generating potentials) through the function *TPP.scripts.heatmap_gen.generate_all_2d_3d_heatmaps*


## Structure Representation






