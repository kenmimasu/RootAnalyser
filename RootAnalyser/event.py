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
    
    def smear_jets(self, smear_function, **kwargs):
        newjets = smear_function(*self.jets, **kwargs)
        self.jets = newjets
        self.bjets = [j for j in self.jets if j.btag]
        self.ljets = [j for j in self.jets if not j.btag]
    
    def smear_MET(self, smear_function, **kwargs):
        self.MET, self.MET_phi = smear_function(self.MET, self.MET_phi, **kwargs)
        
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
