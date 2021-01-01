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


