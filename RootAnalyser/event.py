################################################################################
# Event class with containers for different final states
################################################################################
class Event(): # Basic struct
    def __init__(self):
        # Containers
        self.leptons, self.photons, self.jets = [],[],[]
        self.electrons, self.muons, self.taus = [],[],[]
        self.ljets, self.bjets, self.tops = [],[],[]
        self.higgs, self.ws, self.zs = [],[],[]
        self.nus, self.exotics = [], []
        self.visible = []

        # Counters
        self.npho, self.nlep, self.njet = 0, 0, 0
        self.nele, self.nmu, self.ntau =0, 0, 0
        self.nljet, self.nbjet, self.ntop = 0, 0, 0
        self.nhiggs, self.nz, self.nw = 0, 0, 0
        self.nnu, self.nexo = 0, 0
        # Global event info
        self.ht_tot, self.ht_jet, self.ee_jet, self.MET, self.MET_phi = (0.,)*5
        self.weight=1.
    def __repr__(self):
        # print self.__dict__
        return '''
        npho = {npho}, nlep = {nlep}, njet = {njet}
        nele = {nele}, nmu = {nmu}, ntau = {ntau}
        nljet = {nljet}, nbjet = {nbjet}, ntop = {ntop}
        nhiggs = {nhiggs}, nz = {nz}, nw = {nw}
        nnu = {nnu}, nexo = {nexo}
        '''.format(**self.__dict__)
################################################################################
