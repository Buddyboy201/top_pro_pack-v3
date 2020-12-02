

import TPP.API.atom as atom


class Residue:
    def __init__(self, name, resid, atoms, chain, load_json=False, old_resid=None):
        if not load_json:
            self.name = name
            self.resid = resid
            self.atoms = atoms
            self.centroid = None
            self.chain = chain
            self.old_resid = old_resid
        else:
            self.name = None
            self.resid = None
            self.atoms = None
            self.centroid = None
            self.chain = None

    def get_json_dict(self):
        return {
            "name": self.name,
            "resid": self.resid,
            "atoms": tuple((atm.get_json_dict() for atm in self.atoms)),
            "centroid": self.centroid
        }

    def deserialize_json(self, data):
        pass

    def get_name(self):
        return self.name

    def get_atoms(self):
        return self.atoms

    def get_resid(self):
        return self.resid

    def get_old_resid(self):
        return self.old_resid

    def get_COM(self, exclude_backbone=False):
        if self.name == "GLY":
            for atm in self.atoms:
                if atm.get_name() == "CA":
                    self.centroid = atm.get_coords()
                    return self.centroid
            if self.centroid is None:
                print("No CA atom found for GLY molecule at {}".format(self.resid))
        COM = [0.0, 0.0, 0.0]
        mass_sum = 0
        for atm in self.atoms:
            if exclude_backbone:
                if (atm.get_name() == "CA" or atm.get_name() == "C" or atm.get_name() == "N" or atm.get_name() == "O"):
                    continue
            COM[0] += atm.get_mass() * atm.get_coords()[0]
            COM[1] += atm.get_mass() * atm.get_coords()[1]
            COM[2] += atm.get_mass() * atm.get_coords()[2]
            mass_sum += atm.get_mass()
        COM[0] /= float(mass_sum)
        COM[1] /= float(mass_sum)
        COM[2] /= float(mass_sum)
        self.centroid = tuple([round(i, 3) for i in COM])
        return self.centroid

    def get_centroid(self, exclude_backbone=False):
        if self.centroid == None: self.update_COM(exclude_backbone=exclude_backbone)
        return self.centroid

    def add_atom(self, atm):
        self.atoms.append(atm)

    def update_COM(self, exclude_backbone=False):
        self.centroid = self.get_COM(exclude_backbone=exclude_backbone)

    def get_chain(self):
        return self.chain


