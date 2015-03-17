from ROOT import TPySelector,TH1F
import ROOT
import random, os
import numpy as np
from LHCO import *
from itertools import product, combinations
########################################################################
tag='no_ptCut'
# Acceptance
pt_gam_min, eta_gam_max = 0., 9999999.
pt_ele_min, eta_ele_max = 20., 2.5
pt_mu_min,  eta_mu_max  = 20., 2.7
pt_tau_min, eta_tau_max = 0., 9999999.
pt_jet_min, eta_jet_max = 20., 9999999.
pt_bjet_min, eta_bjet_max = 20., 2.5
# Analysis cuts
pt_b1, pt_b2 = 45.,20.
pt_l1, pt_l2 = 40.,20.
Mll_min, Mll_max = 83.,99.
drbb_min, drbb_max = 0., 2.5
drll_min, drll_max = 0., 1.6
ptbb_min, ptbb_max = 0.,9999999.
# ptb1_min, ptb1_max = 80.,9999999.
ptb1_min, ptb1_max = 0.,9999999.
#ptb2_min, ptb2_max = 40.,9999999.
ptb2_min, ptb2_max = 0.,9999999.
# ptll_min, ptll_max = 110.,9999999.
ptll_min, ptll_max = 0.,9999999.
# ptl1_min, ptl1_max = 70.,9999999.
ptl1_min, ptl1_max = 0.,9999999.
ptl2_min, ptl2_max = 0.,9999999.
htbb_min, htbb_max = 120., 9999999.
htbbll_min, htbbll_max = 260., 9999999.
# signal region
Mbb_min, Mbb_max = 130.,190.
Mbbll_min, Mbbll_max = 340.,420.

########################################################################
# Tuple of histogram init arguments
histargs = (
            ('Mbb_2btag','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_sel','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_0l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_0l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_1l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_1l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_1l_MTW','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_MZ','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ1','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ2','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ3','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ4','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_2j_2l_ptZ5','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_0l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_0l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_1l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_1l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_1l_MTW','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_MET','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_MZ','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ1','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ2','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ3','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ4','Inv. mass of two hardest b jets',12,10.,250.),
            ('Mbb_3j_2l_ptZ5','Inv. mass of two hardest b jets',12,10.,250.),
            
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
                if jet.btag and ( jet.pt > pt_bjet_min ) and ( abs(jet.eta) < eta_bjet_max ): # collect b-tagged jets
                    bjets.append(jet)
                    nbjet+=1
                else: # collect non b-taged jets
                    ljets.append(jet)
                    nljet+=1

        ################################################################
        # Begin Analysis
        two_btag = nbjet == 2
        two_leptons = (nele == 2 and nmu < 1) or (nmu == 2 and nele < 1)
        one_lepton = (nele == 1 and nmu < 1) or (nmu == 1 and nele < 1)
        no_leptons = (nele == 0 and nmu == 0)
        if two_btag:
            # bjets
            sortbjets= pt_sort(bjets)
            bjet1, bjet2 = sortbjets[0], sortbjets[1]
            Mbb = Minv(bjet1,bjet2)
            self.histos['Mbb_2btag'].Fill(Mbb)
            if njet==2: # exactly 2 jet category
                no_additional_jets, one_additional_jet = True,False
            elif nljet==1: # one additional jet category
                if ljets[0].pt < min([b.pt for b in bjets]):
                    no_additional_jets, one_additional_jet = False,True
                else: 
                    no_additional_jets, one_additional_jet = False,False
            else: no_additional_jets, one_additional_jet = False,False 
                     
            if (bjet1.pt > pt_b1) and (bjet2.pt > pt_b2):
                self.histos['Mbb_sel'].Fill(Mbb)
                if no_additional_jets: self.histos['Mbb_2j'].Fill(Mbb)
                if no_additional_jets and no_leptons:
                    self.histos['Mbb_2j_0l'].Fill(Mbb)
                    if MET > 120.: self.histos['Mbb_2j_0l_MET'].Fill(Mbb) 
                if no_additional_jets and one_lepton:
                    self.histos['Mbb_2j_1l'].Fill(Mbb)
                    if MET > 25.: 
                        self.histos['Mbb_2j_1l_MET'].Fill(Mbb)
                        lep = leptons[0]
                        MET_px, MET_py=met_px_py(MET,MET_phi)
                        MTW = (MET + lep.pt) - (lep.px + MET_px)**2- (lep.py + MET_py)**2
                        if MTW < 120.:self.histos['Mbb_2j_1l_MTW'].Fill(Mbb)
                if no_additional_jets and two_leptons:
                    self.histos['Mbb_2j_2l'].Fill(Mbb)
                    if MET < 60.: 
                        self.histos['Mbb_2j_2l_MET'].Fill(Mbb)
                        sortleps=pt_sort(leptons)
                        lep1, lep2 = sortleps[0], sortleps[1]
                        Mll =  Minv(lep1,lep2)
                        if (Mll > Mll_min) and (Mll < Mll_max):
                            self.histos['Mbb_2j_2l_MZ'].Fill(Mbb)
                            ptZ = pT(lep1,lep2)
                            dRbb = deltaR(bjet1,bjet2)
                            if ptZ <= 90.:
                                if dRbb >=0.7 and dRbb <=3.4:
                                     self.histos['Mbb_2j_2l_ptZ1'].Fill(Mbb)
                            if ptZ > 90. and ptZ <= 120.:
                                if dRbb >=0.7 and dRbb <=3.0:
                                     self.histos['Mbb_2j_2l_ptZ2'].Fill(Mbb)
                            if ptZ > 120. and ptZ <= 160.:
                                if dRbb >=0.7 and dRbb <=2.3:
                                     self.histos['Mbb_2j_2l_ptZ3'].Fill(Mbb)
                            if ptZ > 160. and ptZ <= 200.:
                                if dRbb >=0.7 and dRbb <=1.8:
                                     self.histos['Mbb_2j_2l_ptZ4'].Fill(Mbb)
                            if ptZ > 200.:
                                if dRbb < 1.4:
                                     self.histos['Mbb_2j_2l_ptZ5'].Fill(Mbb)
                                     
                if one_additional_jet: self.histos['Mbb_3j'].Fill(Mbb)
                if one_additional_jet and no_leptons:
                    self.histos['Mbb_3j_0l'].Fill(Mbb)
                    if MET > 120.: self.histos['Mbb_3j_0l_MET'].Fill(Mbb)
                if one_additional_jet and one_lepton:
                    self.histos['Mbb_3j_1l'].Fill(Mbb)
                    if MET > 25.: 
                        self.histos['Mbb_3j_1l_MET'].Fill(Mbb)
                        lep = leptons[0]
                        MET_px, MET_py=met_px_py(MET,MET_phi)
                        MTW = (MET + lep.pt) - (lep.px + MET_px)**2- (lep.py + MET_py)**2
                        if MTW < 120.:self.histos['Mbb_2j_1l_MTW'].Fill(Mbb)
                if one_additional_jet and two_leptons:
                    self.histos['Mbb_3j_2l'].Fill(Mbb)
                    if MET < 60.: 
                        self.histos['Mbb_3j_2l_MET'].Fill(Mbb)
                        sortleps=pt_sort(leptons)
                        lep1, lep2 = sortleps[0], sortleps[1]
                        Mll =  Minv(lep1,lep2)
                        if (Mll > Mll_min) and (Mll < Mll_max):
                            self.histos['Mbb_3j_2l_MZ'].Fill(Mbb)
                            ptZ = pT(lep1,lep2)
                            dRbb = deltaR(bjet1,bjet2)
                            if ptZ <= 90.:
                                if dRbb >=0.7 and dRbb <=3.4:
                                     self.histos['Mbb_3j_2l_ptZ1'].Fill(Mbb)
                            if ptZ > 90. and ptZ <= 120.:
                                if dRbb >=0.7 and dRbb <=3.0:
                                     self.histos['Mbb_3j_2l_ptZ2'].Fill(Mbb)
                            if ptZ > 120. and ptZ <= 160.:
                                if dRbb >=0.7 and dRbb <=2.3:
                                     self.histos['Mbb_3j_2l_ptZ3'].Fill(Mbb)
                            if ptZ > 160. and ptZ <= 200.:
                                if dRbb >=0.7 and dRbb <=1.8:
                                     self.histos['Mbb_3j_2l_ptZ4'].Fill(Mbb)
                            if ptZ > 200.:
                                if dRbb < 1.4:
                                     self.histos['Mbb_3j_2l_ptZ5'].Fill(Mbb) 
                                                                
        if nmu+nele!=nlep or nbjet+nljet!=njet: 
            self.Abort('ABORT!' )
        if not entry % 50000: print entry

        return 1

    def SlaveTerminate( self ):
        return 1
    @try_except
    def Terminate( self ):
        fname = self.fChain.GetCurrentFile().GetName()
        file_name= os.path.basename(fname).replace('.root','_analysis_cuts_{}.root'.format(tag))
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
        print 'wrote output to {}'.format(os.path.abspath(file_name))
        return 1
         