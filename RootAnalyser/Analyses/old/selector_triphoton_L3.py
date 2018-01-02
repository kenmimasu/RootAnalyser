########################################################################
# Template selector for LHEF formatted root files
from ROOT import TPySelector,TH1F
import ROOT
import random, os,sys
import numpy as np
from selector_tools import *
from itertools import product, combinations
########################################################################
# Plot data
tag=''
write_histos=False
plot_dir='/Users/Ken/Desktop' # where to write out plot data
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
# Additional cuts
deg2rad = 2*np.pi/360
theta_min, theta_max, dtheta_min = 16.1*deg2rad, 163.9*deg2rad, 20.*deg2rad
########################################################################
# Tuple of histogram init arguments
histargs = (('pTa1','Leading Photon pT',50,0.,50.),
            ('pTa2','2nd Photon pT',50,0.,50.),
            ('pTa3','3rd Photon pT',50,0.,50.),
            ('Ea1','Leading Photon E',50,0.,50.),
            ('Ea2','2nd Photon E',50,0.,50.),
            ('Ea3','3rd Photon E',50,0.,50.),
            ('dRa_h','dR between 2 hardest photons',50,2.,5.),
            ('dRa_s','dR between 2 softest photons',50,0.,3.),
            ('tha1','polar angle of photon 1',60,0.,180.),
            ('tha2','polar angle of photon 2',60,0.,180.),
            ('tha3','polar angle of photon 3',60,0.,180.),
            ('tha1_theta','polar angle of photon 1',60,0.,180.),
            ('tha2_theta','polar angle of photon 2',60,0.,180.),
            ('tha3_theta','polar angle of photon 3',60,0.,180.),
            ('dtheta12','polar angle between photons 1 and 2',50,0.,3.14),
            ('dang12','angle between photons 1 and 2',50,0.,3.14),
            ('dang13','angle between photons 1 and 3',50,0.,3.14),
            ('dang23','angle between photons 2 and 3',50,0.,3.14),
            ('dtheta12_theta','polar angle between photons 1 and 2',50,0.,3.14),
            ('Ea1_theta','Leading Photon E',50,0.,50.),
            ('Ea2_theta','2nd Photon E',50,0.,50.),
            ('Ea3_theta','3rd Photon E',50,0.,50.),
            ('Ea1_dtheta','Leading Photon E',50,0.,50.),
            ('Ea2_dtheta','2nd Photon E',50,0.,50.),
            ('Ea3_dtheta','3rd Photon E',50,0.,50.),
            ('dtheta12_sel','polar angle between photons 1 and 2',50,0.,3.14),
            ('dtheta13_sel','polar angle between photons 1 and 3',50,0.,3.14),
            ('dtheta23_sel','polar angle between photons 2 and 3',50,0.,3.14),
            ('Ea3_sel','Leading Photon E',50,0.,50.),
            ('dtheta13_E3','polar angle between photons 1 and 3',50,0.,3.14),
            ('dtheta23_E3','polar angle between photons 2 and 3',50,0.,3.14),
            ('ct3_E3','polar angle between photons 2 and 3',20,-1.,1.),
            ('dtheta13_ct3','polar angle between photons 1 and 3',50,0.,3.14),
            ('dtheta23_ct3','polar angle between photons 2 and 3',50,0.,3.14),
            ('ct3_ct3','polar angle between photons 2 and 3',20,-1.,1.),
            
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
        self.npass1=0
        self.npass2=0
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
        phot1,phot2,phot3 = pt_sort(photons)
        self.histos['pTa1'].Fill(phot1.pt)
        self.histos['pTa2'].Fill(phot2.pt)
        self.histos['pTa3'].Fill(phot3.pt)
        self.histos['dRa_h'].Fill(deltaR(phot1,phot2))
        self.histos['dRa_s'].Fill(deltaR(phot2,phot3))
        phot4,phot5,phot6 = ee_sort(photons)
        ee1,ee2,ee3 = phot4.ee,phot5.ee,phot6.ee
        ct1, ct2, ct3 = phot4.pz/ee1, phot5.pz/ee2, phot6.pz/ee3
        pol1, pol2, pol3 = np.arccos(ct1), np.arccos(ct2), np.arccos(ct3)
        self.histos['Ea1'].Fill(ee1)
        self.histos['Ea2'].Fill(ee2)
        self.histos['Ea3'].Fill(ee3)        
        self.histos['tha1'].Fill(pol1/deg2rad)
        self.histos['tha2'].Fill(pol2/deg2rad)
        self.histos['tha3'].Fill(pol3/deg2rad)
        dtheta12,dtheta13,dtheta23 = abs(pol1-pol2), abs(pol1-pol3), abs(pol2-pol3)
        dang12, dang13, dang23 = np.arccos(pdot(phot4,phot5)/ee1/ee2), np.arccos(pdot(phot4,phot6)/ee1/ee3), np.arccos(pdot(phot5,phot6)/ee2/ee3)
        self.histos['dang12'].Fill(dang12)
        self.histos['dang13'].Fill(dang13)
        self.histos['dang23'].Fill(dang23)
        self.histos['dtheta12'].Fill(dtheta12)
        if ((pol1 > theta_min) and (pol1 < theta_max) 
        and (pol2 > theta_min) and (pol2 < theta_max) ):
            self.histos['dtheta12_theta'].Fill(dtheta12)
            self.histos['Ea1_theta'].Fill(ee1)
            self.histos['Ea2_theta'].Fill(ee2)
            self.histos['Ea3_theta'].Fill(ee3)
            self.histos['tha1_theta'].Fill(pol1/deg2rad)
            self.histos['tha2_theta'].Fill(pol2/deg2rad)
            if dang12 > dtheta_min:
                self.histos['Ea1_dtheta'].Fill(ee1)
                self.histos['Ea2_dtheta'].Fill(ee2)
                self.histos['Ea3_dtheta'].Fill(ee3)
                if ((dang13 > dtheta_min) and (dang23 > dtheta_min) 
                and (pol3 > theta_min) and (pol3 < theta_max) and (phot6.ee>2.) ):
                    self.histos['dtheta12_sel'].Fill(dang12)
                    self.histos['dtheta13_sel'].Fill(dang13)
                    self.histos['dtheta23_sel'].Fill(dang23)
                    self.histos['Ea3_sel'].Fill(phot6.ee)
                    self.histos['tha3_theta'].Fill(pol3/deg2rad)
                    if ee3 > 0.125*91.:
                        self.histos['dtheta13_E3'].Fill(dang13)
                        self.histos['dtheta23_E3'].Fill(dang23)
                        self.histos['ct3_E3'].Fill(ct3)
                        self.npass1+=1
                        if abs(ct3) < 0.75:
                            self.histos['dtheta13_ct3'].Fill(dang13)
                            self.histos['dtheta23_ct3'].Fill(dang23)
                            self.histos['ct3_ct3'].Fill(ct3)
                        
                            self.npass2+=1
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
            if write_histos: write_hist(hist,'{}/{}_{}.dat'.format(plot_dir,file_name.replace('.root',''),hist.GetName()))
        cur_dir.cd()
        out_file.Close()
        del out_file
        print '{} events passed selection 1'.format(self.npass1)
        print '{} events passed selection 2'.format(self.npass2)
        print 'wrote output to {}'.format(os.path.abspath(file_name))
        return 1
         