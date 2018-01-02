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

########################################################################
write_histos=True
plot_dir='/Users/Ken/GoogleDrive/Desk/ALPs/Analysis/plots' # where to write out plot data
assert os.path.exists(plot_dir),"'plot_dir path doesn't exist!"
########################################################################
# Tuple of histogram init arguments
histargs = (('pTgamma','photon pT',50,330.,500.),
            ('eta_gamma','photon pseudorapidity',50,-2.5,2.5),
            ('ct_gamma','photon polar angle',50,-1.,1.),
            ('th_gamma','photon polar angle',50,0.,180.),
            ('e_gamma','photon energy',30,495.,500.),
            )
########################################################################
    
# def try_except(fn):
#     """decorator for extra debugging output"""
#     def wrapped(*args, **kwargs):
#         try:
#             return fn(*args, **kwargs)
#         except:
#             import traceback
#             traceback.print_exc()
#             assert(0)
#     return wrapped
    
class MyPySelector( TPySelector ):
    def Begin( self ):
        return 1
    def SlaveBegin( self, tree ):
        self.nEvents = tree.GetEntries()
        print 'Processing {} events...'.format(self.nEvents)
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        self.pTcuts = np.arange(0.,110.,10.)
        self.nPass = [0 for p in self.pTcuts]
        self.weights = [0. for p in self.pTcuts]
        self.tot_wgt = 0.
        return 1
        
    @try_except
    def Process( self, entry ):
        # print self.__dict__
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
        self.histos['eta_gamma'].Fill(gamma.eta,evt_wgt)
        self.histos['ct_gamma'].Fill(gamma.costheta,evt_wgt)
        self.histos['th_gamma'].Fill(360./2/np.pi*np.arccos(gamma.costheta),evt_wgt)
        self.histos['e_gamma'].Fill(gamma.ee,evt_wgt)
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
####################################################################################################

