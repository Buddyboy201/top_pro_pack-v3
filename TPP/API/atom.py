
from mendeleev import element
from Bio.PDB import Vector

# TODOs:
# chain handling - handle by seperating out
# pychimera

element_mass = {}

def deserialize_json(data):
    pass

class Atom:
    def __init__(self, symbol, name, atomid, coords, bfactor, load_json=False):
        if not load_json:
            self.symbol = symbol.capitalize()
            self.name = name
            self.atomid = atomid
            self.coords = coords  # (), for consistency save everything as np.array()
            self.mc_sc = False
            if self.name == "CA" or self.name == "C" or self.name == "N" or self.name == "O":
                self.mc_sc = True
            if element_mass.get(self.symbol) is None:
                element_mass[self.symbol] = element(self.symbol).atomic_weight
            self.atomic_mass = element_mass[self.symbol]
            self.vector = Vector(x=self.coords[0], y=self.coords[1], z=self.coords[2])
            self.bfactor = bfactor
        else:
            self.symbol = None
            self.name = None
            self.atomid = None
            self.coords = None
            self.mc_sc = None
            self.atomic_mass = None
            self.vector = None
            self.bfactor = None

    def get_json_dict(self):
        return {
            "symbol": self.symbol,
            "name": self.name,
            "atomid": self.atomid,
            "coords": tuple(self.coords),
            "mc_sc": self.mc_sc,
            "atomic_mass": self.atomic_mass
        }

    def is_mainchain(self):
        return self.mc_sc

    def get_symbol(self):
        return self.symbol

    def get_atomid(self):
        return self.atomid

    def get_name(self):
        return self.name

    def get_coords(self):
        return self.coords

    def get_mass(self):
        return self.atomic_mass

    def get_bfactor(self):
        return self.bfactor

    def get_vector(self):
        return self.vector









