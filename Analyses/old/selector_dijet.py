########################################################################
# Template selector for LHEF formatted root files
from ROOT import TPySelector,TH1F
import ROOT
import random, os,sys
import numpy as np
from selector_tools import *
from itertools import product, combinations
########################################################################
# Acceptance
# pt_gam_min, eta_gam_max = 10., 2.5
# pt_ele_min, eta_ele_max = 10., 2.5
# pt_mu_min,  eta_mu_max  = 10., 2.7
# pt_tau_min, eta_tau_max = 25., 2.5
# pt_jet_min, eta_jet_max = 20., 3.5
pt_gam_min, eta_gam_max = 0., 9999999.
pt_ele_min, eta_ele_max = 0., 9999999.
pt_mu_min,  eta_mu_max  = 0., 9999999.
pt_tau_min, eta_tau_max = 0., 9999999.
pt_jet_min, eta_jet_max = 0., 9999999.
########################################################################
# Define any other global cut variables here
########################################################################
tag='' # string appended to output file name to potentially distinguish different runs
write_histos=False # will write histograms out as text files if you want to use other plotting software
plot_dir='/Users/Ken/Documents/Work/Computing/Axions/axion_dijet/Analysis/Plots' # where to write out plot data
########################################################################
# Tuple of histogram init arguments.
# Entry should be ('hist_name','hist_title',n_bins,x_min,x_max)
histargs = (('theta_j','polar angle',50,0.,3.2),
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
        ################################################################
        non_final_state,final_state=[],[]
        leptons, electrons, muons, taus, photons, jets, ljets, bjets = [],[],[],[],[],[],[],[]
        npho, nlep, nele, nmu, ntau, njet, nljet, nbjet = (0,)*8
        ht_tot, ht_jet, ee_jet, MET_px, MET_py = (0.,)*5
        ################################################################
        current_event = tree.Event
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
        ################################################################
        # Begin Analysis
        if not entry % 10000: print entry
        if entry ==0: self.Abort('ABORT!')
        jet1,jet2 = pt_sort(jet)
        self.histos['pTa1'].Fill(phot1.pt) # This is how you fill a histogram, referring by name
        
        return 1

    def SlaveTerminate( self ):
        return 1
        
    def Terminate( self ):
        fname = self.fChain.GetCurrentFile().GetName()
        # output file will have the same name as input but with '_analysis' appended
        file_name= os.path.basename(fname).replace('.root','_analysis{}.root'.format(tag)) 
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
        print 'wrote output to {}'.format(os.path.abspath(file_name))
        return 1
         