########################################################################
# Template selector for LHCO formatted root files used for detector level studies
# Requires pyROOT, numpy and selector_tools.py
# Python version 2.7 to be safe 
from ROOT import TPySelector,TH1F
import ROOT
import random, os
import numpy as np
from selector_tools import *
from itertools import product, combinations
########################################################################
avail_objects='''
########################################################################
LHCO selector
Available objects
Containers:
    leptons - List of all lepton objects
    electrons,muons,taus,photons - Lists of each particle object
    jets - List of all jets
    bjets,ljets - Lists of b-tagged and non b-tagged jets
Counters:
    npho, nlep, nele, nmu, ntau, njet, nljet, nbjet -
    numbers of each particle type i.e. length of each relevant list
Global event info:
    ht_tot - scalar sum of pT of all visible particles
    ht_jet - scalar sum of pT of all jets
    ee_jet - sum of all jet energies
Missing energy:
    MET - Missing transverse energy
    MET_phi, MET_px, MET_px - azimuthal direction, x component and y component of MET
########################################################################
'''
########################################################################
# Acceptance
# pt_gam_min, eta_gam_max = 10., 2.5
# pt_ele_min, eta_ele_max = 10., 2.5
# pt_mu_min,  eta_mu_max  = 10., 2.7
# pt_tau_min, eta_tau_max = 25., 2.5
# pt_jet_min, eta_jet_max = 20., 3.5
tag=''
pt_gam_min, eta_gam_max = 0., 9999999.
pt_ele_min, eta_ele_max = 0., 9999999.
pt_mu_min,  eta_mu_max  = 0., 9999999.
pt_tau_min, eta_tau_max = 0., 9999999.
pt_jet_min, eta_jet_max = 0., 9999999.
########################################################################
write_histos=True
# write_histos=False
plot_dir='/Users/Ken/GoogleDrive/Desk/HELNLO/MG5/pp_zh_eebb_nogam/Analysis/Plots' # where to write out plot data
assert os.path.exists(plot_dir),"'plot_dir path doesn't exist!"
########################################################################
# Tuple of histogram init arguments
histargs = (
            ('pT_elec','electron pT',40,0.,200.),
            ('pT_posi','positron pT',40,0.,200.),
            ('pT_bjet1','Leading b-jet pT',40,0.,200.),
            ('pT_bjet2','SubLeading b-jet pT',40,0.,200.),
            ('pT_Z','Z boson pT',100,0.,2000.),
            ('eta_elec','electron rapidity',50,-5.,5.),
            ('eta_posi','positron rapidity',50,-5.,5.),
            ('eta_bjet1','Leading b-jet rapidity',50,-5.,5.),
            ('eta_bjet2','SubLeading b-jet rapidity',50,-5.,5.),
            ('eta_Z','Z boson pseudo-rapidity',100,-5.,5.),
            ('y_Z','Z boson rapidity',50,-5.,5.),
            ('m_bb','bbbar invariant mass',80,0.,800.),
            ('m_ll','dilepton invariant mass',40,0.,200.),
            ('m_HV','bbbar-dilepton invariant mass',100,0.,2000.),
            ('R_bb','bbbar separation',60,0.,6.),
            ('R_ll','dilepton separation',60,0.,6.),
            )
########################################################################
def try_except(fn):
    """decorator for extra debugging output"""
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except:
            import traceback
            traceback.print_exc()
            assert(0)
    return wrapped
########################################################################    
class MyPySelector( TPySelector ):
    def Begin( self ):
        return 1
    def SlaveBegin( self, tree ):
        print 'Processing {} events...'.format(tree.GetEntries())
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        [hist.Sumw2() for hist in self.histos.values()]
        return 1
        
    @try_except
    def Process( self, entry ):
        tree = self.fChain
        if tree.GetEntry( entry ) <= 0:
            return 0
        ################################################################
        # Event info containers and counters, made global for access by reader functions
        global leptons, electrons, muons, taus, photons, jets, ljets, bjets
        global npho, nlep, nele, nmu, ntau, njet, nljet, nbjet
        global ht_tot, ht_jet, ee_jet, MET, MET_phi
        leptons, electrons, muons, taus, photons, jets, ljets, bjets = [],[],[],[],[],[],[],[]
        npho, nlep, nele, nmu, ntau, njet, nljet, nbjet = (0,)*8
        ht_tot, ht_jet, ee_jet, MET, MET_phi = (0.,)*5
        ################################################################
        # Read tree and fill info
        if not entry % 10000: print entry
        evt_wgt = read_tree(tree) # returns event weight, set to 1 if no such information is present        
        ################################################################
        # Begin Analysis
        # Book particles
        bjet1, bjet2 = pt_sort(bjets)
        for lep in leptons:
            if lep.charge==-1: elec = lep
            if lep.charge==1: posi = lep
        ################################################################
        # Calculate kinematic quantities
        m_bb, m_ll, m_bbll = Minv(bjet1,bjet2), Minv(elec,posi), Minv(bjet1,bjet2,elec,posi)
        pT_bb, pT_ll = pT(bjet1,bjet2), pT(elec,posi)
        dR_bb, dR_ll = deltaR(bjet1,bjet2), deltaR(elec,posi)
        y_Z, eta_Z = rap(elec,posi), psrap(elec,posi)
        ################################################################
        # if entry ==0: self.Abort('ABORT!')
        #Fill histos
        self.histos['pT_elec'].Fill(elec.pt, evt_wgt) 
        self.histos['pT_posi'].Fill(posi.pt, evt_wgt) 
        self.histos['pT_bjet1'].Fill(bjet1.pt, evt_wgt) 
        self.histos['pT_bjet2'].Fill(bjet2.pt, evt_wgt) 
        self.histos['pT_Z' ].Fill(pT_ll, evt_wgt) 
        self.histos['eta_elec'].Fill(elec.eta, evt_wgt) 
        self.histos['eta_posi'].Fill(posi.eta, evt_wgt) 
        self.histos['eta_bjet1'].Fill(bjet1.eta, evt_wgt) 
        self.histos['eta_bjet2'].Fill(bjet2.eta, evt_wgt) 
        self.histos['eta_Z'].Fill(eta_Z, evt_wgt) 
        self.histos['y_Z'].Fill(y_Z, evt_wgt) 
        self.histos['m_bb'].Fill(m_bb, evt_wgt) 
        self.histos['m_ll'].Fill(m_ll, evt_wgt) 
        self.histos['m_HV'].Fill(m_bbll, evt_wgt) 
        self.histos['R_bb'].Fill(dR_bb, evt_wgt) 
        self.histos['R_ll'].Fill(dR_ll, evt_wgt) 
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
        fname = self.fChain.GetCurrentFile().GetName()
        file_name= os.path.basename(fname).replace('.root','_analysis{}.root'.format(tag))
        cur_dir = ROOT.gDirectory
        out_file = ROOT.TFile(file_name,'RECREATE','test')
        out_file.cd()
        for hist in self.histos.values():
            hist.Write()
            if write_histos: 
                hist_data = '{}/{}_{}.dat'.format( plot_dir, file_name.replace('.root',''), hist.GetName() )
                write_hist(hist, filename = hist_data)
        cur_dir.cd()
        out_file.Close()
        del out_file
        print ''
        print 'wrote root output to {}'.format(os.path.abspath(file_name))
        if write_histos: print 'wrote histogram data to {}'.format(plot_dir)
        print '########################################################################'
        return 1

########################################################################################################
def try_except(fn):
    """decorator for extra debugging output"""
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except:
            import traceback
            traceback.print_exc()
            assert(0)
    return wrapped
########################################################################    
# Reader functions for different TTree structures, LHCO and LHEF only.
def read_tree(TTree):
    tree_name = TTree.GetName()
    if tree_name=='LHCO':
        reader = _read_LHCO
    elif tree_name=='LHEF':
        reader = _read_LHEF
    else:
        print 'Only LHCO and LHEF structures are implemented!'
        sys.exit()
    return reader(TTree)
    
def _read_LHCO(tree):
    global npho, nlep, nele, nmu, ntau, njet, nljet, nbjet
    global ht_tot, ht_jet, ee_jet, MET, MET_phi
    # MET
    MET = tree.MissingET[0].MET
    MET_phi = tree.MissingET[0].Phi
    if (MET_phi < 0.0): MET_phi += 2.0*np.pi
    # Photons
    for i in xrange(tree.Photon_size):
        phot = Photon.TRoot(tree.Photon[i])
        acceptance = ( ( phot.pt > pt_gam_min ) and ( abs(phot.eta) < eta_gam_max ) )
        if acceptance: photons.append(phot); ht_tot+=phot.pt; npho+=1
        
    # Leptons
    for i in xrange(tree.Electron_size):
        elec = Electron.TRoot(tree.Electron[i])
        acceptance = ( ( elec.pt > pt_ele_min ) and ( abs(elec.eta) < eta_ele_max ) )
        if acceptance: electrons.append(elec); leptons.append(elec); ht_tot+=elec.pt; nlep+=1; nele+=1
        
    for i in xrange(tree.Muon_size):
        muon = Muon.TRoot(tree.Muon[i])
        acceptance = ( ( muon.pt > pt_mu_min ) and ( abs(muon.eta) < eta_mu_max ) )
        if acceptance: muons.append(muon); leptons.append(muon); ht_tot+=muon.pt; nlep+=1; nmu+=1
            
    for i in xrange(tree.Tau_size):
        tau = Tau.TRoot(tree.Tau[i])
        acceptance = ( ( tau.pt > pt_tau_min ) and ( abs(tau.eta) < eta_tau_max ) )
        if acceptance: taus.append(tau); ht_tot+=tau.pt; ntau+=1

    # Jets
    for i in xrange(tree.Jet_size):
        jet = Jet.TRoot(tree.Jet[i])
        acceptance = ( ( jet.pt > pt_jet_min ) and ( abs(jet.eta) < eta_jet_max ) )
        if acceptance:
            jets.append(jet)
            ht_tot+=jet.pt; ht_jet+=jet.pt; ee_jet+=jet.ee; njet+=1
            if jet.btag: # collect b-tagged jets
                bjets.append(jet)
                nbjet+=1
            else: # collect non b-tagged jets
                ljets.append(jet)
                nljet+=1
    return 1.
    
def _read_LHEF(tree):
    global npho, nlep, nele, nmu, ntau, njet, nljet, nbjet
    global ht_tot, ht_jet, ee_jet, MET, MET_phi
    MET_px, MET_py = 0., 0.
    non_final_state,final_state=[],[]
    # select only final state particles
    for part in tree.Particle: 
        non_final_state.append(part.Mother1)
        non_final_state.append(part.Mother2)
    non_final_state=set(non_final_state)
    for i,part in enumerate(tree.Particle): 
        if i+1 not in non_final_state:
            final_state.append(part)
    for part in final_state:
        # Photons
        if abs(part.PID) == 22:
            phot = Photon.LHEF(part)
            acceptance = ( ( phot.pt > pt_gam_min ) and ( abs(phot.eta) < eta_gam_max ) )
            if acceptance: photons.append(phot); ht_tot+=phot.pt; npho+=1 
        # Leptons                   
        elif abs(part.PID) == 11:
            elec = Electron.LHEF(part)
            acceptance = ( ( elec.pt > pt_ele_min ) and ( abs(elec.eta) < eta_ele_max ) )
            if acceptance: leptons.append(elec); electrons.append(elec); ht_tot+=elec.pt; nlep+=1; nele+=1
        elif abs(part.PID) == 13:
            muon = Muon.LHEF(part)
            acceptance = ( ( muon.pt > pt_mu_min ) and ( abs(muon.eta) < eta_mu_max ) )
            if acceptance: leptons.append(muon); muons.append(muon); ht_tot+=muon.pt; nlep+=1; nmu+=1
        elif abs(part.PID) == 15:
            tau = Tau.LHEF(part)
            acceptance = ( ( tau.pt > pt_mu_min ) and ( abs(tau.eta) < eta_mu_max ) )
            if acceptance: leptons.append(tau); ht_tot+=tau.pt; ntau+=1
        # Jets (quarks)                 
        elif abs(part.PID) in (1,2,3,4,5):
            jet = Jet.LHEF(part)
            acceptance = ( ( jet.pt > pt_jet_min ) and ( abs(jet.eta) < eta_jet_max ) )
            if acceptance:
                jets.append(jet); ht_tot+=jet.pt; ht_jet+=jet.pt; ee_jet+=jet.ee; njet+=1
                if jet.btag: # collect b-tagged jets
                    bjets.append(jet)
                    nbjet+=1
                else: # collect non b-tagged jets
                    ljets.append(jet)
                    nljet+=1
        # All other PIDs contribute to MET
        else:
            MET_px+=part.Px; MET_py+=part.Py
    # MET
    MET = np.sqrt(MET_px**2+MET_py**2)
    try:
        MET_phi = np.arctan(MET_py/MET_px)
        if (MET_phi < 0.0): MET_phi += 2.0*np.pi
    except ZeroDivisionError:
        MET_phi = 0.
    try:
        weight = tree.Event[0].Weight
    except AttributeError:
        weight = 1.
    return weight
########################################################################################################

