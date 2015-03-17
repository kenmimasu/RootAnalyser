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
tag='_DELPHI'
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
theta_min, theta_max, isol_cone = 16.1*deg2rad, 163.9*deg2rad, 30.*deg2rad
isol_cut= 15.*deg2rad 
########################################################################
# Tuple of histogram init arguments
histargs = (
            ('Ea1','Leading Photon E',50,50.,100.),
            ('Ea2','2nd Photon E',50,50.,100.),
            ('Ea3','3rd Photon E',50,0.,100.),
            ('dRa_h','dR between 2 hardest photons',50,2.,5.),
            ('dRa_s','dR between 2 softest photons',50,0.,3.),
            ('tha1','polar angle of photon 1',60,0.,180.),
            ('tha2','polar angle of photon 2',60,0.,180.),
            ('tha3','polar angle of photon 3',60,0.,180.),
            ('tha1_theta','polar angle of photon 1',60,0.,180.),
            ('tha2_theta','polar angle of photon 2',60,0.,180.),
            ('tha3_theta','polar angle of photon 3',60,0.,180.),
            ('dang12','polar angle between photons 1 and 2',50,0.,3.14),
            ('dang13','polar angle between photons 1 and 3',50,0.,3.14),
            ('dang23','polar angle between photons 2 and 3',50,0.,3.14),
            ('dang12_theta','polar angle between photons 1 and 2',50,0.,3.14),
            ('dang13_theta','polar angle between photons 1 and 3',50,0.,3.14),
            ('dang23_theta','polar angle between photons 2 and 3',50,0.,3.14),
            ('Ea1_theta','Leading Photon E',50,50.,100.),
            ('Ea2_theta','2nd Photon E',50,50.,100.),
            ('Ea3_theta','3rd Photon E',50,0.,100.),
            ('Ea1_isol','Leading Photon E',50,0.,100.),
            ('Ea2_isol','2nd Photon E',50,50.,100.),
            ('Ea3_isol','3rd Photon E',50,0.,100.),
            ('planar_isol','planar angle of 3 photons',50,0.,6.),
            ('Ea3_planar','3rd Photon E',50,0.,100.),
            )
########################################################################
def theta_acceptance(pol):
    poldeg = pol/deg2rad
    return (25. < poldeg < 35.) or (42. < poldeg < 88.) or (92. < poldeg < 138.) or (145. < poldeg < 155.)
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
        self.npass=0
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
        sqrts = tree.Particle[0].E+tree.Particle[1].E
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
        phot1,phot2,phot3 = ee_sort(photons)
        self.histos['dRa_h'].Fill(deltaR(phot1,phot2))
        self.histos['dRa_s'].Fill(deltaR(phot2,phot3))
        ee1,ee2,ee3 = phot1.ee,phot2.ee,phot3.ee
        ct1, ct2, ct3 = phot1.pz/phot1.ee, phot2.pz/phot2.ee, phot3.pz/phot3.ee
        pol1, pol2, pol3 = np.arccos(ct1), np.arccos(ct2), np.arccos(ct3)
        dang12, dang13, dang23 = np.arccos(pdot(phot1,phot2)/ee1/ee2), np.arccos(pdot(phot1,phot3)/ee1/ee3), np.arccos(pdot(phot2,phot3)/ee2/ee3)
        planar_angle = dang12+dang13+dang23
        self.histos['Ea1'].Fill(ee1)
        self.histos['Ea2'].Fill(ee2)
        self.histos['Ea3'].Fill(ee3)
        self.histos['tha1'].Fill(pol1/deg2rad)
        self.histos['tha2'].Fill(pol2/deg2rad)
        self.histos['tha3'].Fill(pol3/deg2rad)
        self.histos['dang12'].Fill(dang12)
        self.histos['dang13'].Fill(dang13)
        self.histos['dang23'].Fill(dang23)
        if theta_acceptance(pol1) and theta_acceptance(pol2) and theta_acceptance(pol3) and min(ee1,ee2,ee3) > 3.:
            self.histos['dang12_theta'].Fill(dang12)
            self.histos['dang13_theta'].Fill(dang13)
            self.histos['dang23_theta'].Fill(dang23)
            self.histos['Ea1_theta'].Fill(ee1)
            self.histos['Ea2_theta'].Fill(ee2)
            self.histos['Ea3_theta'].Fill(ee3)
            self.histos['tha1_theta'].Fill(pol1/deg2rad)
            self.histos['tha2_theta'].Fill(pol2/deg2rad)
            self.histos['tha3_theta'].Fill(pol3/deg2rad)
            if min(ee1,ee2) > 0.15*sqrts and min(dang12,dang13,dang23) > isol_cut and dang12 > 30.*deg2rad:
                self.histos['Ea1_isol'].Fill(ee1)
                self.histos['Ea2_isol'].Fill(ee2)
                self.histos['Ea3_isol'].Fill(ee3)
                self.histos['planar_isol'].Fill(planar_angle)
                if ee3 > 0.06*sqrts and planar_angle > 350.*deg2rad and min(dang13,dang23) > isol_cone:
                    self.histos['Ea3_planar'].Fill(ee3)
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
            if write_histos: write_hist(hist,'{}/{}_{}.dat'.format(plot_dir,file_name.replace('.root',''),hist.GetName()))
        cur_dir.cd()
        out_file.Close()
        del out_file
        print '{} events passed selection'.format(self.npass)
        print 'wrote output to {}'.format(os.path.abspath(file_name))
        return 1
         