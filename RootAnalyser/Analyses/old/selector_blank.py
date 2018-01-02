from ROOT import TPySelector,TH1F
import ROOT
import random, os
import numpy as np
from LHCO import *
from itertools import product, combinations
########################################################################
# Acceptance
# pt_gam_min, eta_gam_max = 10., 2.5
# pt_ele_min, eta_ele_max = 10., 2.5
# pt_mu_min,  eta_mu_max  = 10., 2.5
# pt_tau_min, eta_tau_max = 25., 2.5
# pt_jet_min, eta_jet_max = 20., 3.5
pt_gam_min, eta_gam_max = 0., 9999999.
pt_ele_min, eta_ele_max = 0., 9999999.
pt_mu_min,  eta_mu_max  = 0., 9999999.
pt_tau_min, eta_tau_max = 0., 9999999.
pt_jet_min, eta_jet_max = 0., 9999999.
########################################################################
# Tuple of histogram init arguments
histargs = (
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
        return 1
        
    @try_except
    def Process( self, entry ):
        tree = self.fChain
        if tree.GetEntry( entry ) <= 0:
            return 0
        # ################################################################
        # leptons, electrons, muons, taus, photons, jets, ljets, bjets = [],[],[],[],[],[],[],[]
        # npho, nlep, nele, nmu, ntau, njet, nljet, nbjet = (0,)*8
        # ht_tot, ht_jet, ee_jet = (0.,)*3
        # ################################################################
        # # MET
        # MET = tree.MissingET[0].MET
        # MET_phi = tree.MissingET[0].Phi
        # if (MET_phi < 0.0): MET_phi += 2.0*np.pi
        # # Photons
        # for i in xrange(tree.Photon_size):
        #     phot = Photon.TRoot(tree.Photon[i])
        #     acceptance = ( ( phot.pt > pt_gam_min ) and ( abs(phot.eta) < eta_gam_max ) )
        #     if acceptance: 
        #         photons.append(phot)
        #         ht_tot+=phot.pt; npho+=1
        # 
        # # Leptons
        # for i in xrange(tree.Electron_size):
        #     elec = Electron.TRoot(tree.Electron[i])
        #     acceptance = ( ( elec.pt > pt_ele_min ) and ( abs(elec.eta) < eta_ele_max ) )
        #     if acceptance: 
        #         electrons.append(elec)
        #         leptons.append(elec)
        #         ht_tot+=elec.pt; nlep+=1; nele+=1
        # for i in xrange(tree.Muon_size):
        #     muon = Muon.TRoot(tree.Muon[i])
        #     acceptance = ( ( muon.pt > pt_mu_min ) and ( abs(muon.eta) < eta_mu_max ) )
        #     if acceptance: 
        #         muons.append(muon)
        #         leptons.append(muon)
        #         ht_tot+=muon.pt; nlep+=1; nmu+=1
        # 
        # for i in xrange(tree.Tau_size):
        #     tau = Tau.TRoot(tree.Tau[i])
        #     acceptance = ( ( tau.pt > pt_tau_min ) and ( abs(tau.eta) < eta_tau_max ) )
        #     if acceptance: 
        #         taus.append(tau)
        #         ht_tot+=tau.pt; ntau+=1
        # # Jets
        # for i in xrange(tree.Jet_size):
        #     jet = Jet.TRoot(tree.Jet[i])
        #     acceptance = ( ( jet.pt > pt_jet_min ) and ( abs(jet.eta) < eta_jet_max ) )
        #     if acceptance:
        #         jets.append(jet)
        #         ht_tot+=jet.pt; ht_jet+=jet.pt; ee_jet+=jet.ee; njet+=1
        #         if jet.btag: # collect b-tagged jets
        #             bjets.append(jet)
        #             nbjet+=1
        #         else: # collect non b-taged jets
        #             ljets.append(jet)
        #             nljet+=1
        #             
        #             
        ################################################################
        # Begin Analysis
        if not entry % 10000: print entry

        # if entry ==10: self.Abort('ABORT!')
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
        fname = self.fChain.GetCurrentFile().GetName()
        print os.path.basename(fname).split('.')[0]
        
        # cur_dir = ROOT.gDirectory
        # out_file = ROOT.TFile(file_name,'RECREATE','test')
        # out_file.cd()
        # for h in histargs:
        #     self.histos[h[0]].Write()
        #     # write_hist(h,'{}/{}.dat'.format(os.getcwd(),h.GetName()))
        # cur_dir.cd()
        # out_file.Close()
        # del out_file
        
        
        # print 'wrote output to {}'.format(os.path.abspath(file_name))
        return 1
         