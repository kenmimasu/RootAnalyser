from particles import Photon, Electron, Tau, Muon, Jet, Top, Particle
from event import Event
import numpy as np
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
    'eta_gam_max': 9999999.,
    'eta_ele_max': 9999999.,
    'eta_mu_max': 9999999.,
    'eta_tau_max': 9999999.,
    'eta_jet_max': 9999999.
}
# Reader functions for different TTree structures, LHCO and LHEF only.
################################################################################
def read_tree(TTree, acceptance=None):
    tree_name = TTree.GetName()
    if tree_name=='LHCO':
        reader = _read_LHCO
    elif tree_name=='LHEF':
        reader = _read_LHEF
    else:
        print 'Only LHCO and LHEF structures are implemented!'
        sys.exit()
    return reader(TTree, acceptance)
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
    '''Read in event ROOT tree of LHCO format imposing acceptance cuts specified by "acceptance"
       returning an Event instance'''
    acc = acceptance if acceptance is not None else default_acceptance
    for k in default_acceptance.keys():
        if k not in acc: acc[k]=default_acceptance[k]
        
    evt=Event()
    MET_px, MET_py = 0., 0.
    # select only final state particles
    final_state = filter(lambda p : p.Status > 0, tree.Particle)
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
            acceptance = ( ( tau.pt > acc['pt_mu_min'] ) and ( abs(tau.eta) < acc['eta_mu_max'] ) )
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
        elif abs(part.PID) == 6:
            top = Top.LHEF(part)
            acceptance = ( ( top.pt > acc['pt_jet_min'] ) and ( abs(top.eta) < acc['eta_jet_max'] ) )
            if acceptance:
                evt.tops.append(top)
                evt.ht_tot+=top.pt
                evt.ntop+=1
                
        # All other PIDs contribute to MET and grouped into exotics container
        elif abs(part.PID) not in (23,24,25):
            exo = Particle.LHEF(part)
            exo.PID = part.PID
            evt.exotics.append(exo)
            evt.nexo += 1 
            MET_px+=part.Px; MET_py+=part.Py
    # MET
    evt.MET = np.sqrt(MET_px**2+MET_py**2)
    try:
        evt.MET_phi = np.arctan(MET_py/MET_px)
        if (evt.MET_phi < 0.0): evt.MET_phi += 2.0*np.pi
    except ZeroDivisionError:
        evt.MET_phi = 0.
    try:
        evt.weight = tree.Event[0].Weight
    except AttributeError:
        evt.weight = 1.
        
    return evt
################################################################################