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
LHCO,LHEF selector
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
    MET_phi - azimuthal direction of MET
########################################################################
'''
########################################################################
# Acceptance
# acceptance = {
#     'pt_gam_min ': 10.,
#     'pt_ele_min ': 10.,
#     'pt_mu_min  ': 10.,
#     'pt_tau_min ': 25.,
#     'pt_jet_min ': 20.,
#     'eta_gam_max': 2.5,
#     'eta_ele_max': 2.5,
#     'eta_mu_max ': 2.7,
#     'eta_tau_max': 2.5,
#     'eta_jet_max': 3.5
# }
acceptance=None # no acceptance requirements
########################################################################
tag=''
########################################################################
write_histos=True
plot_dir='/Users/Ken/GoogleDrive/Desk/ALPs/Analysis/plots' # where to write out plot data
assert os.path.exists(plot_dir),"'plot_dir path doesn't exist!"
########################################################################
# Tuple of histogram init arguments
histargs = (('pTgamma','photon pT',100,100.,1000.),
            ('eta_gamma_100','photon pseudorapidity',50,-2.5,2.5),
            ('eta_gamma_150','photon pseudorapidity',50,-2.5,2.5),
            ('eta_gamma_300','photon pseudorapidity',50,-2.5,2.5),
            ('eta_gamma_500','photon pseudorapidity',50,-2.5,2.5),
            ('eta_gamma_800','photon pseudorapidity',50,-2.5,2.5),
            ('ct_gamma','photon polar angle',50,-1.,1.),
            ('phi_gamma','photon azimuthal angle',50,-3.,3.),
            )
########################################################################
class MyPySelector( TPySelector ):
    def Begin( self ):
        return 1
    def SlaveBegin( self, tree ):
        self.nEvents = tree.GetEntries()
        print 'Processing {} events...'.format(self.nEvents)
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        self.pTcuts = np.arange(150.,850.,50.)
        self.nPass = [0 for p in self.pTcuts]
        self.weights = [0. for p in self.pTcuts]
        self.tot_wgt = 0.
        return 1
        
    @try_except
    def Process( self, entry ):
        tree = self.fChain
        if tree.GetEntry( entry ) <= 0:
            return 0
        ################################################################
        # Event info containers and counters, made global for access by reader functions
        # global leptons, electrons, muons, taus, photons, jets, ljets, bjets
        # global npho, nlep, nele, nmu, ntau, njet, nljet, nbjet
        # global ht_tot, ht_jet, ee_jet, MET, MET_phi
        # leptons, electrons, muons, taus, photons, jets, ljets, bjets = [],[],[],[],[],[],[],[]
        # npho, nlep, nele, nmu, ntau, njet, nljet, nbjet = (0,)*8
        # ht_tot, ht_jet, ee_jet, MET, MET_phi = (0.,)*5
        ################################################################
        # Read tree and fill info
        # if not entry % 10000: print entry
        evt = read_tree(tree) # returns selector_tools.Event instance
        evt_wgt = evt.weight
        self.tot_wgt += evt_wgt
        ################################################################
        # Begin Analysis
        # if entry ==10: self.Abort('ABORT!')
        gamma = evt.photons[0]
        pta = gamma.pt
        self.histos['pTgamma'].Fill(pta,evt_wgt)
        self.histos['eta_gamma_100'].Fill(gamma.eta,evt_wgt)
        for cut in [100.,150.,300.,500.,800.]:
            if pta > cut:
                hist_name = 'eta_gamma_{}'.format(int(cut))
                self.histos[hist_name].Fill(gamma.eta,evt_wgt)
        
        self.histos['ct_gamma'].Fill(gamma.costheta,evt_wgt)
        self.histos['phi_gamma'].Fill(gamma.phi,evt_wgt)
        for i,cut in enumerate(self.pTcuts):
            if pta > cut: 
                self.nPass[i]+=1
                self.weights[i]+=evt_wgt
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
        fname = self.fChain.GetCurrentFile().GetName()
        file_name= os.path.basename(fname).replace('.root','_analysis_cuts{}.root'.format(tag))
        cur_dir = ROOT.gDirectory
        out_file = ROOT.TFile(file_name,'RECREATE','test')
        out_file.cd()
        for h in histargs:
            hist= self.histos[h[0]]
            hist.Write()
            if write_histos: write_hist(hist,'{}/{}_{}.dat'.format(plot_dir,file_name.replace('.root',''),hist.GetName()))
        cur_dir.cd()
        out_file.Close()
        del out_file
        print ''
        print 'wrote root output to {}'.format(os.path.abspath(file_name))
        if write_histos: print 'wrote histogram data to {}'.format(plot_dir)
        print '########################################################################'
        print 'Total cross section (pb): ', self.tot_wgt
        print 'pT Cut Efficiencies:'
        # print self.pTcuts,self.nPass,self.weights
        for cut,count,wgt in zip(self.pTcuts,self.nPass,self.weights):
            print '{}: {} (unweighted); {} (weighted) => {} pb'.format(cut, float(count)/float(self.nEvents),wgt/self.tot_wgt,wgt )
        return 1

########################################################################################################
# Reader functions for different TTree structures, LHCO and LHEF only.
# def read_tree(TTree):
#     tree_name = TTree.GetName()
#     if tree_name=='LHCO':
#         reader = _read_LHCO
#     elif tree_name=='LHEF':
#         reader = _read_LHEF
#     else:
#         print 'Only LHCO and LHEF structures are implemented!'
#         sys.exit()
#     return reader(TTree)
#
# def _read_LHCO(tree):
#     global npho, nlep, nele, nmu, ntau, njet, nljet, nbjet
#     global ht_tot, ht_jet, ee_jet, MET, MET_phi
#     # MET
#     MET = tree.MissingET[0].MET
#     MET_phi = tree.MissingET[0].Phi
#     if (MET_phi < 0.0): MET_phi += 2.0*np.pi
#     # Photons
#     for i in xrange(tree.Photon_size):
#         phot = Photon.TRoot(tree.Photon[i])
#         acceptance = ( ( phot.pt > pt_gam_min ) and ( abs(phot.eta) < eta_gam_max ) )
#         if acceptance: photons.append(phot); ht_tot+=phot.pt; npho+=1
#
#     # Leptons
#     for i in xrange(tree.Electron_size):
#         elec = Electron.TRoot(tree.Electron[i])
#         acceptance = ( ( elec.pt > pt_ele_min ) and ( abs(elec.eta) < eta_ele_max ) )
#         if acceptance: electrons.append(elec); leptons.append(elec); ht_tot+=elec.pt; nlep+=1; nele+=1
#
#     for i in xrange(tree.Muon_size):
#         muon = Muon.TRoot(tree.Muon[i])
#         acceptance = ( ( muon.pt > pt_mu_min ) and ( abs(muon.eta) < eta_mu_max ) )
#         if acceptance: muons.append(muon); leptons.append(muon); ht_tot+=muon.pt; nlep+=1; nmu+=1
#
#     for i in xrange(tree.Tau_size):
#         tau = Tau.TRoot(tree.Tau[i])
#         acceptance = ( ( tau.pt > pt_tau_min ) and ( abs(tau.eta) < eta_tau_max ) )
#         if acceptance: taus.append(tau); ht_tot+=tau.pt; ntau+=1
#
#     # Jets
#     for i in xrange(tree.Jet_size):
#         jet = Jet.TRoot(tree.Jet[i])
#         acceptance = ( ( jet.pt > pt_jet_min ) and ( abs(jet.eta) < eta_jet_max ) )
#         if acceptance:
#             jets.append(jet)
#             ht_tot+=jet.pt; ht_jet+=jet.pt; ee_jet+=jet.ee; njet+=1
#             if jet.btag: # collect b-tagged jets
#                 bjets.append(jet)
#                 nbjet+=1
#             else: # collect non b-tagged jets
#                 ljets.append(jet)
#                 nljet+=1
#     return 1.
#
# def _read_LHEF(tree):
#     global npho, nlep, nele, nmu, ntau, njet, nljet, nbjet
#     global ht_tot, ht_jet, ee_jet, MET, MET_phi
#     MET_px, MET_py = 0., 0.
#     non_final_state,final_state=[],[]
#     # select only final state particles
#     for part in tree.Particle:
#         non_final_state.append(part.Mother1)
#         non_final_state.append(part.Mother2)
#     non_final_state=set(non_final_state)
#     for i,part in enumerate(tree.Particle):
#         if i+1 not in non_final_state:
#             final_state.append(part)
#     for part in final_state:
#         # Photons
#         if abs(part.PID) == 22:
#             phot = Photon.LHEF(part)
#             acceptance = ( ( phot.pt > pt_gam_min ) and ( abs(phot.eta) < eta_gam_max ) )
#             if acceptance: photons.append(phot); ht_tot+=phot.pt; npho+=1
#         # Leptons
#         elif abs(part.PID) == 11:
#             elec = Electron.LHEF(part)
#             acceptance = ( ( elec.pt > pt_ele_min ) and ( abs(elec.eta) < eta_ele_max ) )
#             if acceptance: leptons.append(elec); electrons.append(elec); ht_tot+=elec.pt; nlep+=1; nele+=1
#         elif abs(part.PID) == 13:
#             muon = Muon.LHEF(part)
#             acceptance = ( ( muon.pt > pt_mu_min ) and ( abs(muon.eta) < eta_mu_max ) )
#             if acceptance: leptons.append(muon); muons.append(muon); ht_tot+=muon.pt; nlep+=1; nmu+=1
#         elif abs(part.PID) == 15:
#             tau = Tau.LHEF(part)
#             acceptance = ( ( tau.pt > pt_mu_min ) and ( abs(tau.eta) < eta_mu_max ) )
#             if acceptance: leptons.append(tau); ht_tot+=tau.pt; ntau+=1
#         # Jets (quarks)
#         elif abs(part.PID) in (1,2,3,4,5):
#             jet = Jet.LHEF(part)
#             acceptance = ( ( jet.pt > pt_jet_min ) and ( abs(jet.eta) < eta_jet_max ) )
#             if acceptance:
#                 jets.append(jet); ht_tot+=jet.pt; ht_jet+=jet.pt; ee_jet+=jet.ee; njet+=1
#                 if jet.btag: # collect b-tagged jets
#                     bjets.append(jet)
#                     nbjet+=1
#                 else: # collect non b-tagged jets
#                     ljets.append(jet)
#                     nljet+=1
#         # All other PIDs contribute to MET
#         else:
#             MET_px+=part.Px; MET_py+=part.Py
#     # MET
#     MET = np.sqrt(MET_px**2+MET_py**2)
#     try:
#         MET_phi = np.arctan(MET_py/MET_px)
#         if (MET_phi < 0.0): MET_phi += 2.0*np.pi
#     except ZeroDivisionError:
#         MET_phi = 0.
#     try:
#         weight = tree.Event[0].Weight
#     except AttributeError:
#         weight = 1.
#     return weight
########################################################################################################

