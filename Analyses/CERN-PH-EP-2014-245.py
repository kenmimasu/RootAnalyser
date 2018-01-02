########################################################################
# Selector for LHEF/LHCO formatted root files used for parton/
# detector level studies.
# Requires pyROOT, numpy and selector_tools.py
# Python version 2.7 to be safe 
from ROOT import TPySelector,TH1F
import ROOT
import random, os
import numpy as np
from RootAnalyser.event import Event
from RootAnalyser.root import *
from RootAnalyser.kinematics import *
from RootAnalyser.readers import read_tree
from itertools import product, combinations
########################################################################
# CERN-PH-EP-2014-245: ATLAS 8 TeV monophoton + MET analysis
########################################################################
avail_objects='''
########################################################################
LHCO,LHEF selector
Available data members of Event instance
Containers:
    leptons - List of all lepton objects
    electrons,muons,taus,photons - Lists of each particle object
    jets - List of all jets
    bjets,ljets - Lists of b-tagged and non b-tagged jets
Counters:
    npho, nlep, nele, nmu, ntau, njet, nljet, nbjet -
    numbers of each particle type i.e. length of corresponding list
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
acceptance={
    'pt_gam_min' : 0.,
    'pt_ele_min' : 7.,
    'pt_mu_min'  : 6.,
    'pt_tau_min' : 0.,
    'pt_jet_min' : 30.,
    'eta_gam_max': 9999999.,
    'eta_ele_max': 2.47,
    'eta_mu_max' : 2.5,
    'eta_tau_max': 9999999.,
    'eta_jet_max': 9999999.}
tag='test'
########################################################################
MET_cut = 150.
etaa_cut = 1.37
pta_trig = 125.
# hadem_cut = 0.05
dphi_MET_cut = 0.4
jet_dR_thres = 0.4
# lep_dR_thres = 0.5

########################################################################
write_histos=True
plot_dir='/Users/Ken/GoogleDrive/python/Analyses/plots' # where to write out plot data
assert os.path.exists(plot_dir),"'plot_dir path doesn't exist!"
########################################################################
# Tuple of histogram init arguments
histargs = (
            ('presel','MET',50,0.,600.),
            # ('hadem','Photon H/E',50,0.,0.05),
            ('pTa','Photon pT',50,0.,600.),
            ('etaa','Photon pseudorapidity',50,-2.,2.),
            ('dphi','dPhi(a,MET)',50,0.,3.14),
            ('pTa_veto','Photon pT',50,0.,600.),
            ('MET_veto','MET',50,0.,600.),
            ('veto_pTj','jet pT',50,0.,100.),
            ('veto_dRaj','dR(a,j)',50,0.,3.),
            ('veto_etaj','jet rapidity',50,-5.,5.),
            ('veto_pTl','jet pT',50,0.,100.),
            ('veto_dRal','dR(a,l)',50,0.,3.),
            ('veto_etal','lepton rapidity',50,-5.,5.),
            )
########################################################################
class MyPySelector( TPySelector ):
    def Begin( self ):
        return 1
    def SlaveBegin( self, tree ):
        self.nEvents = tree.GetEntries()
        print 'Processing {} events...'.format(self.nEvents)
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        self.pTcuts = [125.]
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
        # Read tree and fill info
        evt = read_tree(tree, acceptance=acceptance) # returns selector_tools.Event instance
        evt_wgt = evt.weight
        self.tot_wgt += evt_wgt
        if entry ==10: self.Abort('ABORT!')
        ################################################################
        # Begin Analysis
        # Begin Analysis
        self.histos['presel'].Fill(evt.MET)
        veto=False
        if npho>=1 and MET > MET_cut and nljets < 2:
            photon = pt_sort(evt.photons)[0]
            # self.histos['hadem'].Fill(photon.hadem)
            if (abs(photon.eta) < etaa_cut) and(photon.pt > pta_trig):
                self.histos['pTa'].Fill(photon.pt)
                self.histos['etaa'].Fill(photon.eta)
                dphi_MET = dphi(photon,ref_phi=evt.MET_phi)
                self.histos['dphi'].Fill(dphi_MET)
                veto_ctr=0
                for jet in evt.jets:
                    dRaj = deltaR(photon,jet)
                    veto = ((jet.pt > jet_pt_thres) and (dRaj > jet_dR_thres))
                    if veto: 
                        veto_ctr+=1
                    if veto_ctr==2:
                        self.histos['veto_pTj'].Fill(jet.pt)
                        self.histos['veto_dRaj'].Fill(dRaj)
                        self.histos['veto_etaj'].Fill(jet.eta)
                        break
                for lep in evt.leptons:
                    dRaj = deltaR(photon,jet)
                    veto = ((lep.pt > lep_pt_thres) and
                            (dRaj > lep_dR_thres))
                    if veto: 
                        self.histos['veto_pTl'].Fill(lep.pt)
                        self.histos['veto_dRal'].Fill(dRaj)
                        self.histos['veto_etal'].Fill(lep.eta)
                        break
                if not veto and dphi_MET > dphi_MET_cut:
                    self.histos['pTa_veto'].Fill(photon.pt)
                    self.histos['MET_veto'].Fill(MET)
                    for i,cut in enumerate(self.pTcuts):
                        if photon.pt > cut:self.npass[i]+=1
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

