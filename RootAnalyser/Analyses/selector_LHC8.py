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
pt_gam_min, eta_gam_max = 10., 2.5
pt_ele_min, eta_ele_max = 10., 2.5
pt_mu_min,  eta_mu_max  = 10., 2.7
pt_tau_min, eta_tau_max = 25., 2.5
pt_jet_min, eta_jet_max = 20., 3.5
tag=''
# pt_gam_min, eta_gam_max = 0., 9999999.
# pt_ele_min, eta_ele_max = 0., 9999999.
# pt_mu_min,  eta_mu_max  = 0., 9999999.
# pt_tau_min, eta_tau_max = 0., 9999999.
# pt_jet_min, eta_jet_max = 0., 9999999.
# cuts
########################################################################
MET_cut = 140.
etaa_cut = 1.4442
pta_trig = 145.
hadem_cut = 0.05
dphi_MET_cut = 2.
ptcuts=[145.,160.,190.,250.,400.,700.]
jet_pt_thres, jet_dR_thres = 30., 0.5
lep_pt_thres, lep_dR_thres = 10., 0.5
########################################################################
write_histos=True
plotdir,resdir = '{}/Plots'.format(os.getcwd()),'{}/Results'.format(os.getcwd())
if not os.path.exists(plotdir): os.mkdir(plotdir)
if not os.path.exists(resdir): os.mkdir(resdir)
########################################################################
# Tuple of histogram init arguments
histargs = (
            ('presel','MET',50,0.,600.),
            ('hadem','Photon H/E',50,0.,0.05),
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
    
class MyPySelector( TPySelector ):
    def Begin( self ):
        return 1
    def SlaveBegin( self, tree ):
        self.ntot=tree.GetEntries()
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        self.npass=[0,0,0,0,0,0]
        print 'Processing {} events...'.format(self.ntot)
        return 1
        
    @try_except
    def Process( self, entry ):
        tree = self.fChain
        if tree.GetEntry( entry ) <= 0:
            return 0
        ################################################################
        leptons, electrons, muons, taus, photons, jets, ljets, bjets = [],[],[],[],[],[],[],[]
        npho, nlep, nele, nmu, ntau, njet, nljet, nbjet = (0,)*8
        ht_tot, ht_jet, ee_jet = (0.,)*3
        ################################################################
        # if entry ==20: self.Abort('ABORT!')
        try:
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
        except AttributeError as e:
            return 1
        ################################################################                    
        if not entry % 10000: print entry                    
        ################################################################
        # Begin Analysis
        self.histos['presel'].Fill(MET)
        veto=False
        if npho>=1 and MET > MET_cut:
            photon = pt_sort(photons)[0]
            self.histos['hadem'].Fill(photon.hadem)
            if (abs(photon.eta) < etaa_cut) and (photon.hadem < hadem_cut) and(photon.pt > pta_trig):
                self.histos['pTa'].Fill(photon.pt)
                self.histos['etaa'].Fill(photon.eta)
                dphi_MET = dphi(photon,ref_phi=MET_phi)
                self.histos['dphi'].Fill(dphi_MET)
                veto_ctr=0
                for jet in jets:
                    dRaj = deltaR(photon,jet)
                    veto = ((jet.pt > jet_pt_thres) and (dRaj > jet_dR_thres))
                    if veto: 
                        veto_ctr+=1
                    if veto_ctr==2:
                        self.histos['veto_pTj'].Fill(jet.pt)
                        self.histos['veto_dRaj'].Fill(dRaj)
                        self.histos['veto_etaj'].Fill(jet.eta)
                        break
                for lep in leptons:
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
                    for i,cut in enumerate(ptcuts):
                        if photon.pt > cut:self.npass[i]+=1
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
        fname = os.path.basename(self.fChain.GetCurrentFile().GetName())
        file_name= '{}/{}'.format(resdir,fname.replace('.root','_analysis_cuts{}.root'.format(tag)))
        cur_dir = ROOT.gDirectory
        out_file = ROOT.TFile(file_name,'RECREATE','test')
        out_file.cd()
        for h in histargs:
            hist= self.histos[h[0]]
            hist.Write()
            if write_histos: write_hist(hist,'{}/{}_{}.dat'.format(plotdir,fname.replace('.root',''),hist.GetName()))
        cur_dir.cd()
        out_file.Close()
        del out_file
        print ''
        for i,cut in enumerate(ptcuts):
            print '{} events with pT > {} ({:.2f}%)'.format(self.npass[i],cut, 100.*self.npass[i]/self.ntot)
        print 'wrote root output to {}'.format(os.path.abspath(file_name))
        if write_histos: print 'wrote histogram data to {}'.format(plotdir)
        print '########################################################################'
        return 1
         