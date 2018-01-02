from ROOT import TPySelector,TH1F
import ROOT
import random, os
import numpy as np
from LHCO import *
from itertools import product, combinations
########################################################################
# Create Plots directory
if not os.path.exists('{}/Plots'.format(os.getcwd())): os.mkdir('Plots')
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
########################################################################
# cuts
MET_cut = 130.
etaa_cut = 1.48
jet_pt_thres, jet_dR_thres, jet_eta_cut = 40., 0.5, 3.
########################################################################
# Tuple of histogram init arguments
histargs = (
            ('presel','MET',50,0.,600.),
            ('pTa','Photon pT',50,0.,600.),
            ('pTa_veto','Photon pT',50,0.,600.),
            ('MET_veto','MET',50,0.,600.),
            ('veto_pTj','jet pT',50,0.,100.),
            ('veto_dRaj','dR(a,j)',50,0.,3.),
            ('veto_etaj','jet rapidity',50,-5.,5.)
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
        print 'Processing {} events...'.format(tree.GetEntries())
        self.histos = {arg[0]:TH1F(*arg) for arg in histargs}
        self.npass = 0
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
        # MET
        MET = tree.MissingET[0].MET
        MET_phi = tree.MissingET[0].Phi
        if (MET_phi < 0.0): MET_phi += 2.0*np.pi
        # Photons
        for i in xrange(tree.Photon_size):
            phot = Photon.TRoot(tree.Photon[i])
            acceptance = ( ( phot.pt > pt_gam_min ) and ( abs(phot.eta) < eta_gam_max ) )
            if acceptance: 
                photons.append(phot)
                ht_tot+=phot.pt; npho+=1
        
        # Leptons
        for i in xrange(tree.Electron_size):
            elec = Electron.TRoot(tree.Electron[i])
            acceptance = ( ( elec.pt > pt_ele_min ) and ( abs(elec.eta) < eta_ele_max ) )
            if acceptance: 
                electrons.append(elec)
                leptons.append(elec)
                ht_tot+=elec.pt; nlep+=1; nele+=1
        for i in xrange(tree.Muon_size):
            muon = Muon.TRoot(tree.Muon[i])
            acceptance = ( ( muon.pt > pt_mu_min ) and ( abs(muon.eta) < eta_mu_max ) )
            if acceptance: 
                muons.append(muon)
                leptons.append(muon)
                ht_tot+=muon.pt; nlep+=1; nmu+=1
        
        for i in xrange(tree.Tau_size):
            tau = Tau.TRoot(tree.Tau[i])
            acceptance = ( ( tau.pt > pt_tau_min ) and ( abs(tau.eta) < eta_tau_max ) )
            if acceptance: 
                taus.append(tau)
                ht_tot+=tau.pt; ntau+=1
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
                else: # collect non b-taged jets
                    ljets.append(jet)
                    nljet+=1
                    
                    
        ################################################################
        # Begin Analysis
        if not entry % 10000: print entry
        self.histos['presel'].Fill(MET)
        veto=False
        if npho==1 and MET > MET_cut:
            photon = photons[0]
            if abs(photon.eta) < 1.48 and (photon.hadem < 0.05):
                self.histos['pTa'].Fill(photon.pt)
                for jet in jets:
                    dRaj = deltaR(photon,jet)
                    veto = ((jet.pt > jet_pt_thres) and
                            (abs(jet.eta) < jet_eta_cut) and
                            (dRaj > jet_dR_thres))
                    if veto: 
                        self.histos['veto_pTj'].Fill(jet.pt)
                        self.histos['veto_dRaj'].Fill(dRaj)
                        self.histos['veto_etaj'].Fill(jet.eta)
                        break
                if not veto:
                    self.histos['pTa_veto'].Fill(photon.pt)
                    self.histos['MET_veto'].Fill(MET)
                    self.npass+=1
                
        # if entry ==10: self.Abort('ABORT!')
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
        fname = self.fChain.GetCurrentFile().GetName()
        file_name= os.path.basename(fname).replace('.root','_analysis{}.root'.format(tag))
        cur_dir = ROOT.gDirectory
        out_file = ROOT.TFile(file_name,'RECREATE','test')
        out_file.cd()
        for h in histargs:
            hist= self.histos[h[0]]
            hist.Write()
            write_hist(hist,'{}/Plots/{}_{}.dat'.format(os.getcwd(),file_name.replace('.root',''),hist.GetName()))
        cur_dir.cd()
        out_file.Close()
        del out_file
        print '{} events passed selection'.format(self.npass)
        print 'wrote output to {}'.format(os.path.abspath(file_name))
        return 1
         