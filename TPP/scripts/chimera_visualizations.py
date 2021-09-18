from pathlib import Path
from subprocess import check_call
from threading import Thread


def _display_chimera(row, pdb_directory, chimera_path):
    #chimera_path = Path(r"C:\Program Files\Chimera 1.15rc\bin\chimera.exe")
    chimera_path = Path(chimera_path)
    residues = ",".join(row[1].split(";"))
    path_to_cmd = Path.cwd() / Path("{}.cmd".format(row[0]))
    path_to_pdb = Path(pdb_directory) / Path("{}.pdb".format(row[0]))
    with open(path_to_cmd, "w") as file:
        file.writelines(["select: {}\n".format(residues), "display: {}\n".format(residues),
                         "focus: {}\n".format(residues)])
    if not path_to_pdb.exists():
        return Exception("Invalid pdb file path/file already exists")
    p = check_call([str(chimera_path), str(path_to_pdb), str(path_to_cmd)], shell=True)
    path_to_cmd.unlink()


def display_chimera(conn, sql_id, pdb_directory, chimera_path):
    stmt_id = "SELECT pdbname, oldresid FROM cliques WHERE id={}".format(sql_id)
    row = list(conn.execute(stmt_id))[0]
    t = Thread(target=_display_chimera, args=(row, pdb_directory, chimera_path,))
    t.start()
