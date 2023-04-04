from mendeleev import element

# TODOs:
# chain handling - handle by seperating out

element_mass = {}


class Atom:
    def __init__(self, symbol, name, atomid, coords):
        self.symbol = symbol.capitalize()
        self.name = name
        self.atomid = atomid
        self.coords = coords  # (), for consistency save everything as np.array()
        self.mc_sc = False
        if (
            self.name == "CA"
            or self.name == "C"
            or self.name == "N"
            or self.name == "O"
        ):
            self.mc_sc = True
        if element_mass.get(self.symbol) is None:
            element_mass[self.symbol] = element(self.symbol).atomic_weight
        self.atomic_mass = element_mass[self.symbol]

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
