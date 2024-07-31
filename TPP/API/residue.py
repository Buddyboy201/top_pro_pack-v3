

class Residue:
    def __init__(self, name, resid, atoms, chain, conf_score, old_resid=None):
        self.name = name
        self.resid = resid
        self.atoms = atoms
        self.centroid = None
        self.chain = chain
        self.old_resid = old_resid
        self.conf_score = conf_score
        self.layerinfo = None

    def _get_COM(self, exclude_backbone=False):
        if self.name == "GLY":
            centroid_tmp = None
            for atm in self.atoms:
                if atm.get_name() == "CA":
                    centroid_tmp = atm.get_coords()
            # self.centroid = centroid_tmp
            if centroid_tmp is None:
                print("No CA atom found for GLY molecule at {}".format(self.resid))
            return centroid_tmp
        COM = [0.0, 0.0, 0.0]
        mass_sum = 0
        for atm in self.atoms:
            if exclude_backbone:
                if atm.is_mainchain():
                    continue
            COM[0] += atm.get_mass() * atm.get_coords()[0]
            COM[1] += atm.get_mass() * atm.get_coords()[1]
            COM[2] += atm.get_mass() * atm.get_coords()[2]
            mass_sum += atm.get_mass()
        if mass_sum <= 0:
            return None
        COM[0] /= float(mass_sum)
        COM[1] /= float(mass_sum)
        COM[2] /= float(mass_sum)
        # self.centroid = tuple([round(i, 3) for i in COM])
        return tuple([round(i, 3) for i in COM])

    def get_centroid(self, exclude_backbone=False):
        if self.centroid is None:
            self.update_COM(exclude_backbone=exclude_backbone)
        return self.centroid

    def add_atom(self, atm):
        self.atoms.append(atm)

    def update_COM(self, exclude_backbone=False):
        self.centroid = self._get_COM(exclude_backbone=exclude_backbone)


