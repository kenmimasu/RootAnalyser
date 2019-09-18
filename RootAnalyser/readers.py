from particles import Photon, Electron, Tau, Muon, Jet, Top, Particle
from event import Event
import numpy as np
import sys
################################################################################
'''
########################################################################
LHCO, LHEF selector
Available data members of Event instance, referenced by e.g. evt.photons
Containers:
    leptons - List of all lepton objects
    electrons, muons, taus, photons - Lists of each particle object
    jets - List of all jets
    bjets, ljets - Lists of b-tagged and non b-tagged jets
    exotics - Lists of all particles not identified above 
             (ONLY for LHE event format & they also contribute to MET)
Counters:
    npho, nlep, nele, nmu, ntau, njet, nljet, nbjet, nexo -
    numbers of each particle type i.e. length of corresponding list
Global event info:
    ht_tot - scalar sum of pT of all visible particles
    ht_jet - scalar sum of pT of all jets
    ee_jet - sum of all jet energies
    weight - evt weight, information only present in LHEF
Missing energy:
    MET - Missing transverse energy
    MET_phi - azimuthal direction of MET
########################################################################
'''
################################################################################
# No acceptance requirements
default_acceptance ={
    'pt_gam_min' : 0.,
    'pt_ele_min' : 0.,
    'pt_mu_min'  : 0.,
    'pt_tau_min' : 0.,
    'pt_jet_min' : 0.,
    'pt_bjet_min' : 0.,
    'eta_gam_max': 9999999.,
    'eta_ele_max': 9999999.,
    'eta_mu_max': 9999999.,
    'eta_tau_max': 9999999.,
    'eta_jet_max': 9999999.,
    'eta_bjet_max': 9999999.
}
# Reader functions for different TTree structures, LHCO and LHEF only.
################################################################################
from collections import Callable

class RAReader(Callable):
    '''
    Callable object that returns Event instance.
    '''
    def __init__(self, acceptance=None):
        self.acc = acceptance if acceptance is not None else default_acceptance
        for k in default_acceptance.keys():
            if k not in self.acc: 
                self.acc[k] = default_acceptance[k]
        
        self.event = Event()
    @staticmethod
    def within_acceptance(particle, pt_min, abs_eta_max):
        return  ( ( particle.pt > self.acc[pt_min] ) and 
                  ( abs(particle.eta) < self.acc['abs_eta_max'] ) )
    
    def add_MET(self, met, met_phi):
        self.event.MET = met
        self.event.MET_phi = met_phi
        if (evt.MET_phi < 0.0): evt.MET_phi += 2.0*np.pi
        
    def add_photon(self, phot):
        if self.within_acceptance(phot,'pt_gam_min','eta_gam_max'): 
            self.event.photons.append(phot)
            self.event.ht_tot+=phot.pt
            self.event.npho+=1
    
    def add_lepton(self, lep):
        self.event.leptons.append(lep)
        self.event.ht_tot+=lep.pt
        self.event.nlep+=1
        
    def add_electron(self, elec):
        if self.within_acceptance(elec,'pt_ele_min','eta_ele_max'): 
            self.add_lepton(elec)
            self.event.electrons.append(elec)
            self.event.nele+=1
    
    def add_muon(self, muon):
        if self.within_acceptance(muon,'pt_mu_min','eta_mu_max'): 
            self.add_lepton(muon)
            self.event.muons.append(muon)
            self.event.nmu+=1
    
    def add_tau(self, tau):
        if self.within_acceptance(tau,'pt_tau_min','eta_tau_max'): 
            self.taus.append(tau)
            self.ht_tot+=tau.pt
            self.ntau+=1
            
    def add_jet(self, jet):
        if jet.btag: # collect b-tagged jets
            acceptance = self.within_acceptance(jet,'pt_bjet_min','eta_bjet_max')
        else:
            acceptance = self.within_acceptance(jet,'pt_jet_min','eta_jet_max')

        if acceptance:
            self.event.jets.append(jet)
            self.event.ht_tot+=jet.pt
            self.event.ht_jet+=jet.pt
            self.event.ee_jet+=jet.ee
            self.event.njet+=1
            if jet.btag: # collect b-tagged jets
                self.event.bjets.append(jet)
                self.event.nbjet+=1
            else: # collect non b-tagged jets
                self.event.ljets.append(jet)
                self.event.nljet+=1
    
    def read_tree(self):
        raise NotImplementedError
        
    def __call__(self, tree):
        try:
            return self.read_tree(tree)
        except Exception as e:
            import traceback
            raise Exception(traceback.format_exc())

class LHCO_reader(RAReader):
    def read_tree(self, tree):
        # MET (Sometimes LHCO tree has no attribute MissingET)
        try:
            MET = tree.MissingET[0].MET
            MET_phi = tree.MissingET[0].Phi
            self.add_MET(MET, MET_phi)
        except AttributeError:
            self.add_MET(0., 0.)
        
        # Photons
        for i in xrange(tree.Photon_size):
            phot = Photon.TRoot(tree.Photon[i])
            self.add_photon(phot)
        
        # Leptons
        for i in xrange(tree.Electron_size):
            elec = Electron.TRoot(tree.Electron[i])
            self.add_electron(elec)
        
        for i in xrange(tree.Muon_size):
            muon = Muon.TRoot(tree.Muon[i])
            self.add_muon(muon)
        
        for i in xrange(tree.Tau_size):
            tau = Tau.TRoot(tree.Tau[i])
            self.add_tau(tau)
            
        # Jets
        for i in xrange(tree.Jet_size):
            jet = Jet.TRoot(tree.Jet[i])
            self.add_jet(jet)
        
        # Empty exotics stuff
        self.event.exotics=[]
        self.event.nexo = 0
        self.event.weight = 1.
    
        return self.event

class LHEF_reader(RAReader):
    def __init__(self, *args, **kwargs):
        self.MET_px, self.MET_py = 0., 0.
        super(LHEF_reader, self).__init__(*args, **kwargs)
        
    def add_neutrino(self, nu):
        self.event.nus.append(nu)
        self.MET_px += part.Px
        self.MET_py += part.Py
        self.event.nnu+=1
    
    def add_top(self, top):
        # no acceptance requirements
        self.event.tops.append(top)
        self.event.ht_tot+=top.pt
        self.event.ntop+=1
    
    def add_Z(self, zbos):
        # no acceptance requirements
        self.event.zs.append(zbos)
        self.event.ht_tot+=zbos.pt
        self.event.nz+=1
    
    def add_W(self, wbos):
        # no acceptance requirements
        self.event.ws.append(wbos)
        self.event.ht_tot+=wbos.pt
        self.event.nw+=1
    
    def add_Higgs(self, higg):
        # no acceptance requirements
        self.event.higgs.append(higg)
        self.event.ht_tot+=higg.pt
        self.event.nhiggs+=1
    
    def add_exotic(self, exo):
        # no acceptance requirements
        self.event.exotics.append(exo)
        self.event.nexo += 1 
        self.MET_px += part.Px
        self.MET_py += part.Py
    
    def add_MET_from_px_py(self):
        self.event.MET = np.sqrt( self.MET_px**2 + self.MET_py**2 )
        try:
            self.event.MET_phi = np.arctan2(MET_py,MET_px)
            if (self.event.MET_phi < 0.0): 
                self.event.MET_phi += 2.0*np.pi
        except ZeroDivisionError:
            self.event.MET_phi = 0.
                    
    def read_tree(self, tree):
        MET_px, MET_py = 0., 0.
        # select only final state particles
        self.event.ECM = np.sqrt( (tree.Particle[0].E+tree.Particle[1].E)**2 - 
                           (tree.Particle[0].Pz+tree.Particle[1].Pz)**2 )
        final_state = filter(lambda p : p.Status == 1 , tree.Particle)
        for part in final_state:
            # Photons
            if abs(part.PID) == 22:
                phot = Photon.LHEF(part)
                self.add_photon(phot)
                
            # Charged leptons                   
            elif abs(part.PID) == 11:
                elec = Electron.LHEF(part)
                self.add_electron(elec)
                
            elif abs(part.PID) == 13:
                muon = Muon.LHEF(part)
                self.add_muon(muon)
                
            elif abs(part.PID) == 15:
                tau = Tau.LHEF(part)
                self.add_tau(tau)
                
            # Neutrinos
            elif abs(part.PID) in (12,14,16):
                nu = Particle.LHEF(part)
                self.add_neutrino(nu)
                
            # Jets (quarks & gluons)                 
            elif abs(part.PID) in (1,2,3,4,5,21):
                jet = Jet.LHEF(part)
                self.add_jet(jet)
                    
            elif abs(part.PID) == 6:
                # check if t or tbar
                anti = True if part.PID==-6 else False
                top = Top.LHEF(part, anti=anti)
                self.add_top(top)
        
            elif abs(part.PID) == 23:
                zbos = Particle.LHEF(part)
                self.add_Z(zbos)
        
            elif abs(part.PID) == 24:
                wbos = Particle.LHEF(part)
                self.add_W(wbos)
        
            elif abs(part.PID) == 25:
                higg = Particle.LHEF(part)
                self.add_Higgs(higg)
                
            # All other PIDs contribute to MET and put in exotics container
            else:
                exo = Particle.LHEF(part)
                exo.PID = part.PID
                self.add_exotic(exo)
       
        # MET
        self.add_MET_from_px_py()
        
        # Weight
        try:
            self.event.weight = tree.Event[0].Weight
        except AttributeError:
            self.event.weight = 1.
    
        self.event.rwgt = []
        try:
            for w in tree.Rwgt:
                self.event.rwgt.append(w.Weight)
        except AttributeError:
            pass
        
        return self.event

class Delphes_reader(RAReader):
    def __init__(self, *args, **kwargs):
        super(LHEF_reader, self).__init__(*args, **kwargs)
        
    def read_tree(self, tree):
        # MET (Sometimes LHCO tree has no attribute MissingET)
        # MET (Sometimes LHCO tree has no attribute MissingET)
        try:
            MET = tree.MissingET[0].MET
            MET_phi = tree.MissingET[0].Phi
            self.add_MET(MET, MET_phi)
        except AttributeError:
            self.add_MET(0., 0.)
            
        # Photons
        for i in xrange(tree.Photon_size):
            phot = Photon.Delphes(tree.Photon[i])
            self.add_photon(phot)
            
        # Leptons
        for i in xrange(tree.Electron_size):
            elec = Electron.Delphes(tree.Electron[i])
            self.add_electron(elec)
            
        for i in xrange(tree.Muon_size):
            muon = Muon.Delphes(tree.Muon[i])
            self.add_muon(muon)
        # Jets
        # check Jet container name
        try:
            jet_size, jet_array = tree.KTjet_size, tree.KTjet
        except AttributeError:
            jet_size, jet_array = tree.Jet_size, tree.Jet
    
        for i in xrange(jet_size):
            # check Tau tag
            if not jet_array[i].TauTag:
                jet = Jet.Delphes(jet_array[i])
                self.add_jet(jet)
            else:
                tau = Tau.Delphes(jet_array[i])
                self.add_tau(tau)
    
        # Empty exotics stuff
        self.event.exotics=[]
        self.event.nexo = 0
        
        try:
            self.event.weight = tree.Event[0].Weight
        except AttributeError:
            self.event.weight = 1.
        
        return self.event
        
################################################################################
def read_tree(TTree, acceptance=None):
    tree_name = TTree.GetName()
    if tree_name=='LHCO':
        reader = LHCO_reader(acceptance=acceptance)
    elif tree_name=='LHEF':
        reader = LHEF_reader(acceptance=acceptance)
    elif tree_name=='Delphes':
        reader = Delphes_reader(acceptance=acceptance)
    else:
        print 'Only LHCO, LHEF and Delphes structures are implemented!'
        sys.exit()
        
    return reader(TTree)

def read_tree_old(TTree, acceptance=None):
    tree_name = TTree.GetName()
    if tree_name=='LHCO':
        reader = _read_LHCO
    elif tree_name=='LHEF':
        reader = _read_LHEF
    elif tree_name=='Delphes':
        reader = _read_Delphes
    else:
        print 'Only LHCO, LHEF and Delphes structures are implemented!'
        sys.exit()
    return reader(TTree, acceptance=acceptance)

################################################################################
def _read_LHCO(tree, acceptance=None):
    '''Read in event ROOT tree of LHCO format imposing acceptance cuts specified by "acceptance"
       returning an Event instance'''
    acc = acceptance if acceptance is not None else default_acceptance
    for k in default_acceptance.keys():
        if k not in acc: acc[k]=default_acceptance[k]
        
    evt = Event() # initialise
    # MET (Sometimes LHCO tree has no attribute MissingET)
    try:
        evt.MET = tree.MissingET[0].MET
        evt.MET_phi = tree.MissingET[0].Phi
    except AttributeError:
        evt.MET = 0.
        evt.MET_phi=0.
        
    if (evt.MET_phi < 0.0): evt.MET_phi += 2.0*np.pi
    # Photons
    for i in xrange(tree.Photon_size):
        phot = Photon.TRoot(tree.Photon[i])
        acceptance = ( ( phot.pt > acc['pt_gam_min'] ) and ( abs(phot.eta) < acc['eta_gam_max'] ) )
        if acceptance: 
            evt.photons.append(phot)
            evt.ht_tot+=phot.pt
            evt.npho+=1
        
    # Leptons
    for i in xrange(tree.Electron_size):
        elec = Electron.TRoot(tree.Electron[i])
        acceptance = ( ( elec.pt > acc['pt_ele_min'] ) and ( abs(elec.eta) < acc['eta_ele_max'] ) )
        if acceptance: 
            evt.electrons.append(elec)
            evt.leptons.append(elec)
            evt.ht_tot+=elec.pt
            evt.nlep+=1
            evt.nele+=1
        
    for i in xrange(tree.Muon_size):
        muon = Muon.TRoot(tree.Muon[i])
        acceptance = ( ( muon.pt > acc['pt_mu_min'] ) and ( abs(muon.eta) < acc['eta_mu_max'] ) )
        if acceptance: 
            evt.muons.append(muon)
            evt.leptons.append(muon)
            evt.ht_tot+=muon.pt
            evt.nlep+=1
            evt.nmu+=1
            
    for i in xrange(tree.Tau_size):
        tau = Tau.TRoot(tree.Tau[i])
        acceptance = ( ( tau.pt > acc['pt_tau_min'] ) and ( abs(tau.eta) < acc['eta_tau_max'] ) )
        if acceptance: 
            evt.taus.append(tau)
            evt.ht_tot+=tau.pt
            evt.ntau+=1

    # Jets
    for i in xrange(tree.Jet_size):
        jet = Jet.TRoot(tree.Jet[i])
        acceptance = ( ( jet.pt > acc['pt_jet_min'] ) and ( abs(jet.eta) < acc['eta_jet_max'] ) )
        if acceptance:
            evt.jets.append(jet)
            evt.ht_tot+=jet.pt
            evt.ht_jet+=jet.pt
            evt.ee_jet+=jet.ee
            evt.njet+=1
            if jet.btag: # collect b-tagged jets
                evt.bjets.append(jet)
                evt.nbjet+=1
            else: # collect non b-tagged jets
                evt.ljets.append(jet)
                evt.nljet+=1
    # Empty exotics stuff
    evt.exotics=[]
    evt.nexo = 0
    evt.weight = 1.
    
    return evt
    
################################################################################
def _read_LHEF(tree, acceptance=None):
    '''Read in event ROOT tree of LHEF format imposing acceptance cuts specified by "acceptance"
       returning an Event instance'''
    acc = acceptance if acceptance is not None else default_acceptance
    for k in default_acceptance.keys():
        if k not in acc: acc[k]=default_acceptance[k]
    evt=Event()
    MET_px, MET_py = 0., 0.
    # select only final state particles
    evt.ECM = np.sqrt( (tree.Particle[0].E+tree.Particle[1].E)**2 - 
                       (tree.Particle[0].Pz+tree.Particle[1].Pz)**2 )
    final_state = filter(lambda p : p.Status == 1 , tree.Particle)
    for part in final_state:
        # Photons
        
        if abs(part.PID) == 22:
            phot = Photon.LHEF(part)
            acceptance = ( ( phot.pt > acc['pt_gam_min'] ) and ( abs(phot.eta) < acc['eta_gam_max'] ) )
            if acceptance: 
                evt.photons.append(phot)
                evt.ht_tot+=phot.pt
                evt.npho+=1 
        # Charged leptons                   
        elif abs(part.PID) == 11:
            elec = Electron.LHEF(part)
            acceptance = ( ( elec.pt > acc['pt_ele_min'] ) and ( abs(elec.eta) < acc['eta_ele_max'] ) )
            if acceptance: 
                evt.electrons.append(elec)
                evt.leptons.append(elec)
                evt.ht_tot+=elec.pt
                evt.nlep+=1
                evt.nele+=1
        elif abs(part.PID) == 13:
            muon = Muon.LHEF(part)
            acceptance = ( ( muon.pt > acc['pt_mu_min'] ) and ( abs(muon.eta) < acc['eta_mu_max'] ) )
            if acceptance: 
                evt.muons.append(muon)
                evt.leptons.append(muon)
                evt.ht_tot+=muon.pt
                evt.nlep+=1
                evt.nmu+=1
        elif abs(part.PID) == 15:
            tau = Tau.LHEF(part)
            acceptance = ( ( tau.pt > acc['pt_tau_min'] ) and ( abs(tau.eta) < acc['eta_tau_max'] ) )
            if acceptance: 
                evt.taus.append(tau)
                evt.ht_tot+=tau.pt
                evt.ntau+=1
        # Neutrinos
        elif abs(part.PID) in (12,14,16):
            nu = Particle.LHEF(part)
            evt.nus.append(nu)
            MET_px+=part.Px; MET_py+=part.Py
            evt.nnu+=1
        # Jets (quarks & gluons)                 
        elif abs(part.PID) in (1,2,3,4,5,21):
            jet = Jet.LHEF(part)
            if jet.btag: # collect b-tagged jets
                acceptance = ( ( jet.pt > acc['pt_bjet_min'] ) and 
                               ( abs(jet.eta) < acc['eta_bjet_max'] ) )
            else:
                acceptance = ( ( jet.pt > acc['pt_jet_min'] ) and 
                               ( abs(jet.eta) < acc['eta_jet_max'] ) )
            if acceptance:
                evt.jets.append(jet)
                evt.ht_tot+=jet.pt
                evt.ht_jet+=jet.pt
                evt.ee_jet+=jet.ee
                evt.njet+=1
                if jet.btag: # collect b-tagged jets
                    evt.bjets.append(jet)
                    evt.nbjet+=1
                else: # collect non b-tagged jets
                    evt.ljets.append(jet)
                    evt.nljet+=1
                    
        elif abs(part.PID) == 6:
            # check if t or tbar
            anti = True if part.PID==-6 else False
            top = Top.LHEF(part, anti=anti)
            # no acceptance requirements
            acceptance = True
            if acceptance:
                evt.tops.append(top)
                evt.ht_tot+=top.pt
                evt.ntop+=1
        
        elif abs(part.PID) == 23:
            zbos = Particle.LHEF(part)
            evt.zs.append(zbos)
            evt.ht_tot+=zbos.pt
            evt.nz+=1
        
        elif abs(part.PID) == 24:
            wbos = Particle.LHEF(part)
            evt.ws.append(wbos)
            evt.ht_tot+=wbos.pt
            evt.nw+=1
        
        elif abs(part.PID) == 25:
            higg = Particle.LHEF(part)
            evt.higgs.append(higg)
            evt.ht_tot+=higg.pt
            evt.nhiggs+=1
                
        # All other PIDs contribute to MET and grouped into exotics container
        else:
            exo = Particle.LHEF(part)
            exo.PID = part.PID
            evt.exotics.append(exo)
            evt.nexo += 1 
            MET_px+=part.Px; MET_py+=part.Py
    # MET
    evt.MET = np.sqrt(MET_px**2+MET_py**2)
    try:
        evt.MET_phi = np.arctan2(MET_py,MET_px)
        if (evt.MET_phi < 0.0): evt.MET_phi += 2.0*np.pi
    except ZeroDivisionError:
        evt.MET_phi = 0.
    # Weight
    try:
        evt.weight = tree.Event[0].Weight
    except AttributeError:
        evt.weight = 1.
    
    evt.rwgt = []
    try:
        for w in tree.Rwgt:
            evt.rwgt.append(w.Weight)
    except AttributeError:
        pass
        
    return evt
################################################################################
def _read_Delphes(tree, acceptance=None):
    '''Read in event ROOT tree of Delphes format imposing acceptance cuts specified by "acceptance"
       returning an Event instance'''
    # print sorted([ d for d in dir(tree) if not d.startswith('_')])
    # print tree
    # print dir(tree.GetListOfBranches())
    # for i in xrange(tree.GetListOfBranches().GetEntries()):
    #     print tree.GetListOfBranches()[i]
    # l1 = sorted([ d for d in dir(tree.Muon[1]) if not d.startswith('_')])
    # # l2 = sorted([ d for d in dir(tree.EFlowElectron[1]) if not d.startswith('_')])
    #
    # print l1
    # sys.exit()
    
    
    acc = acceptance if acceptance is not None else default_acceptance
    for k in default_acceptance.keys():
        if k not in acc: acc[k]=default_acceptance[k]
    evt = Event() # initialise
    # MET (Sometimes LHCO tree has no attribute MissingET)
    try:
        evt.MET = tree.MissingET[0].MET
        evt.MET_phi = tree.MissingET[0].Phi
    except AttributeError:
        evt.MET = 0.
        evt.MET_phi=0.
        
    if (evt.MET_phi < 0.0): evt.MET_phi += 2.0*np.pi
    # Photons
    for i in xrange(tree.Photon_size):
        phot = Photon.Delphes(tree.Photon[i])
        acceptance = ( ( phot.pt > acc['pt_gam_min'] ) and ( abs(phot.eta) < acc['eta_gam_max'] ) )
        if acceptance: 
            evt.photons.append(phot)
            evt.ht_tot+=phot.pt
            evt.npho+=1
        
    # Leptons
    for i in xrange(tree.Electron_size):
        elec = Electron.Delphes(tree.Electron[i])
        acceptance = ( ( elec.pt > acc['pt_ele_min'] ) and ( abs(elec.eta) < acc['eta_ele_max'] ) )
        if acceptance: 
            evt.electrons.append(elec)
            evt.leptons.append(elec)
            evt.ht_tot+=elec.pt
            evt.nlep+=1
            evt.nele+=1
        
    for i in xrange(tree.Muon_size):
        muon = Muon.Delphes(tree.Muon[i])
        acceptance = ( ( muon.pt > acc['pt_mu_min'] ) and ( abs(muon.eta) < acc['eta_mu_max'] ) )
        if acceptance: 
            evt.muons.append(muon)
            evt.leptons.append(muon)
            evt.ht_tot+=muon.pt
            evt.nlep+=1
            evt.nmu+=1
            
    # for i in xrange(tree.Tau_size):
    #     tau = Tau.Delphes(tree.Tau[i])
    #     acceptance = ( ( tau.pt > acc['pt_tau_min'] ) and ( abs(tau.eta) < acc['eta_tau_max'] ) )
    #     if acceptance:
    #         evt.taus.append(tau)
    #         evt.ht_tot+=tau.pt
    #         evt.ntau+=1

    # Jets
    # check Jet container name
    try:
        jet_size, jet_array = tree.KTjet_size, tree.KTjet
    except AttributeError:
        jet_size, jet_array = tree.Jet_size, tree.Jet
    
    for i in xrange(jet_size):
        # check Tau tag
        if not jet_array[i].TauTag:
            jet = Jet.Delphes(jet_array[i])
            acceptance = ( ( jet.pt > acc['pt_jet_min'] ) and ( abs(jet.eta) < acc['eta_jet_max'] ) )
            if acceptance:
                evt.jets.append(jet)
                evt.ht_tot+=jet.pt
                evt.ht_jet+=jet.pt
                evt.ee_jet+=jet.ee
                evt.njet+=1
                if jet.btag: # collect b-tagged jets
                    evt.bjets.append(jet)
                    evt.nbjet+=1
                else: # collect non b-tagged jets
                    evt.ljets.append(jet)
                    evt.nljet+=1
        else:
            tau = Tau.Delphes(jet_array[i])
            acceptance = ( ( tau.pt > acc['pt_tau_min'] ) and ( abs(tau.eta) < acc['eta_tau_max'] ) )
            if acceptance: 
                evt.taus.append(tau)
                evt.ht_tot+=tau.pt
                evt.ntau+=1
    
    # Empty exotics stuff
    evt.exotics=[]
    evt.nexo = 0
    
        # try:
        #     for i in xrange(tree.KTjet_size):
        #         jet = Jet.Delphes(tree.KTjet[i])
        #         acceptance = ( ( jet.pt > acc['pt_jet_min'] ) and ( abs(jet.eta) < acc['eta_jet_max'] ) )
        #         if acceptance:
        #             evt.jets.append(jet)
        #             evt.ht_tot+=jet.pt
        #             evt.ht_jet+=jet.pt
        #             evt.ee_jet+=jet.ee
        #             evt.njet+=1
        #             if jet.btag: # collect b-tagged jets
        #                 evt.bjets.append(jet)
        #                 evt.nbjet+=1
        #             else: # collect non b-tagged jets
        #                 evt.ljets.append(jet)
        #                 evt.nljet+=1
        #
        # except AttributeError:
        #     for i in xrange(tree.Jet_size):
        #         jet = Jet.Delphes(tree.Jet[i])
        #         acceptance = ( ( jet.pt > acc['pt_jet_min'] ) and ( abs(jet.eta) < acc['eta_jet_max'] ) )
        #         if acceptance and not :
        #             evt.jets.append(jet)
        #             evt.ht_tot+=jet.pt
        #             evt.ht_jet+=jet.pt
        #             evt.ee_jet+=jet.ee
        #             evt.njet+=1
        #             if jet.btag: # collect b-tagged jets
        #                 evt.bjets.append(jet)
        #                 evt.nbjet+=1
        #             else: # collect non b-tagged jets
        #                 evt.ljets.append(jet)
        #                 evt.nljet+=1
        #     # Empty exotics stuff
        #     evt.exotics=[]
        #     evt.nexo = 0
        
    try:
        evt.weight = tree.Event[0].Weight
    except AttributeError:
        evt.weight = 1.
        
    return evt
################################################################################