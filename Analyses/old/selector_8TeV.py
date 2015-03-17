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
pt_ele_min, eta_ele_max = 0., 2.5
pt_mu_min,  eta_mu_max  = 0., 2.7
pt_tau_min, eta_tau_max = 0., 9999999.
pt_jet_min, eta_jet_max = 0., 9999999.
pt_bjet_min, eta_bjet_max = 0., 2.5
# Analysis cuts
pt_b1, pt_b2 = 40.,20.
pt_l1, pt_l2 = 40.,20.
Mll_min, Mll_max = 80.,100.
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
            ('PtBJet1','Pt of hardest b jet',50,0.,300.),
            ('PtBJet2','Pt of 2nd hardest b jet',50,0.,200.),
            ('Ptbb','Pt of bb system',50,0.,300.),
            ('Ptbb_ht','Pt of bb system',50,0.,300.),
            ('Ptbb_dr','Pt of bb system',50,0.,300.),
            ('Mbb','Inv. mass of two hardest b jets',50,0.,300.),
            ('dRbb','delta R of two hardest b jets',50,0.,4.5),
            ('dRbb_htbb','delta R of two hardest b jets',50,0.,4.5),
            ('dRbb_htbbll','delta R of two hardest b jets',50,0.,4.5),
            ('HTbb','HT of two hardest b jets',50,0.,350.),
            ('PtLep1','Pt of hardest lepton',50,0.,250.),
            ('PtLep2','Pt of 2nd hardest lepton',50,0.,150.),
            ('Ptll','Pt of di-lepton system',50,0.,350.),
            ('Ptll_ht','Pt of di-lepton system',50,0.,350.),
            ('Ptll_dr','Pt of di-lepton system',50,0.,350.),
            ('Ptbbll','Pt of bb system + Pt of di-lepton system',50,0.,600.),
            ('Ptbbll_ht','Pt of bb system + Pt of di-lepton system',50,0.,600.),
            ('Ptbbll_dr','Pt of bb system + Pt of di-lepton system',50,0.,600.),
            ('Mll','Inv. mass of two hardest leptons',50,0.,200.),
            ('dRll','delta R of two hardest leptons',50,0.,4.5),
            ('dRll_htbb','delta R of two hardest leptons',50,0.,4.5),
            ('dRll_htbbll','delta R of two hardest leptons',50,0.,4.5),
            ('HTll','HT of two of hardest leptons',50,0.,300.),
            ('Mbbll','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('HTbbll','HT of two hardest bjets and leptons',50,0.,500.),
            ('MET','Missing Energy',50,0.,200.),
            ('Mbb_pt','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_pt','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_ht','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_ht','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_Mll','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_Mll','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_drll','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_drll','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_drbb','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_drbb','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_ptb1','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_ptb1','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_ptb2','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_ptb2','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_ptl1','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_ptl1','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_ptl2','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_ptl2','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_ptbb','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_ptbb','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_ptll','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_ptll','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_htbb','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_htbb','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_htbbll','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_htbbll','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_mbb','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_mbb','Inv. mass of two hardest bjets and leptons',50,0.,600.),
            ('Mbb_mbbll','Inv. mass of two hardest b jets',50,0.,300.),
            ('Mbbll_mbbll','Inv. mass of two hardest bjets and leptons',50,0.,600.),
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
        one_jet = njet > 0
        two_btag = nbjet > 1
        two_leptons = (nele > 1) or (nmu > 1)
        if two_btag and two_leptons:
            # bjets
            sortbjets= pt_sort(bjets)
            bjet1, bjet2 = sortbjets[0], sortbjets[1]
            # leptons
            sortleps= pt_sort(leptons)
            lep1, lep2 = sortleps[0], sortleps[1]
            bjet1_pt, bjet2_pt, lep1_pt, lep2_pt = bjet1.pt, bjet2.pt, lep1.pt, lep2.pt

            # histos
            if (bjet1_pt > pt_b1) and (bjet2_pt > pt_b2) and (lep1_pt > pt_l1) and (lep2_pt > pt_l2):
                ptbb, ptll, = pT(bjet1,bjet2),pT(lep1,lep2)
                ptbbll=ptbb+ptll
                htbb, htll, htbbll = bjet1_pt+bjet2_pt, lep1_pt+lep2_pt, bjet1_pt+bjet2_pt+lep1_pt+lep2_pt
                Mbb, Mll, Mbbll = Minv(bjet1,bjet2), Minv(lep1,lep2), Minv(bjet1,bjet2,lep1,lep2)
                drbb, drll = deltaR(bjet1,bjet2), deltaR(lep1,lep2)
                self.histos['PtBJet1'].Fill(bjet1_pt)
                self.histos['PtBJet2'].Fill(bjet2_pt)
                self.histos['Mbb'].Fill(Mbb)
                self.histos['dRbb'].Fill(drbb)
                self.histos['Ptbb'].Fill(ptbb)
                self.histos['HTbb'].Fill(htbb)
                self.histos['PtLep1'].Fill(lep1_pt)
                self.histos['PtLep2'].Fill(lep2_pt)
                self.histos['Mll'].Fill(Mll)
                self.histos['dRll'].Fill(drll)
                self.histos['Ptll'].Fill(ptll)
                self.histos['Ptbbll'].Fill(ptbbll)
                self.histos['HTll'].Fill(htll)
                # bbll system
                self.histos['Mbbll'].Fill(Mbbll)
                self.histos['HTbbll'].Fill(htbbll)
                self.histos['MET'].Fill(MET)
                if (bjet1_pt > ptb1_min) and (bjet2_pt > ptb2_min):
                    self.histos['Mbb_pt'].Fill(Mbb)
                    self.histos['Mbbll_pt'].Fill(Mbbll) 
                if (htbb > htbb_min) and (htbb < htbb_max):                    
                    self.histos['Mbb_ht'].Fill(Mbb)
                    self.histos['Mbbll_ht'].Fill(Mbbll)               
                # perform cuts one by one
                if (Mll > Mll_min) and (Mll < Mll_max):
                    self.histos['Mbb_Mll'].Fill(Mbb)
                    self.histos['Mbbll_Mll'].Fill(Mbbll)
                    if (bjet1_pt > ptb1_min) and (bjet1_pt < ptb1_max):
                        self.histos['Mbb_ptb1'].Fill(Mbb)
                        self.histos['Mbbll_ptb1'].Fill(Mbbll)
                        if (bjet1_pt > ptb1_min) and (bjet1_pt < ptb1_max):
                            self.histos['Mbb_ptb2'].Fill(Mbb)
                            self.histos['Mbbll_ptb2'].Fill(Mbbll)
                            if (lep1_pt > ptl1_min) and (lep1_pt < ptl1_max):
                                self.histos['Mbb_ptl1'].Fill(Mbb)
                                self.histos['Mbbll_ptl1'].Fill(Mbbll)
                                if (lep2_pt > ptl1_min) and (lep2_pt < ptl1_max):
                                    self.histos['Mbb_ptl2'].Fill(Mbb)
                                    self.histos['Mbbll_ptl2'].Fill(Mbbll)
                                    if (htbb > htbb_min) and (htbb < htbb_max):
                                        self.histos['Mbb_htbb'].Fill(Mbb)
                                        self.histos['Mbbll_htbb'].Fill(Mbbll)
                                        self.histos['dRbb_htbb'].Fill(drbb)
                                        self.histos['dRll_htbb'].Fill(drll)
                                        if (htbbll > htbbll_min) and (htbbll < htbbll_max):
                                            self.histos['Mbb_htbbll'].Fill(Mbb)
                                            self.histos['Mbbll_htbbll'].Fill(Mbbll)
                                            self.histos['dRbb_htbbll'].Fill(drbb)
                                            self.histos['dRll_htbbll'].Fill(drll)
                                            self.histos['Ptbb_ht'].Fill(ptbb)
                                            self.histos['Ptll_ht'].Fill(ptll)
                                            self.histos['Ptbbll_ht'].Fill(ptbbll)
                                            if (drll > drll_min) and (drll < drll_max):
                                                self.histos['Mbb_drll'].Fill(Mbb)
                                                self.histos['Mbbll_drll'].Fill(Mbbll)
                                                if (drbb > drbb_min) and (drbb < drbb_max):
                                                    self.histos['Mbb_drbb'].Fill(Mbb)
                                                    self.histos['Mbbll_drbb'].Fill(Mbbll)
                                                    self.histos['Ptbb_dr'].Fill(ptbb)
                                                    self.histos['Ptll_dr'].Fill(ptll)      
                                                    self.histos['Ptbbll_dr'].Fill(ptbbll)
                                                    if (ptbb > ptbb_min) and (ptbb < ptbb_max):
                                                        self.histos['Mbb_ptbb'].Fill(Mbb)
                                                        self.histos['Mbbll_ptbb'].Fill(Mbbll)
                                                        if (ptll > ptll_min) and (ptll < ptll_max):
                                                            self.histos['Mbb_ptll'].Fill(Mbb)
                                                            self.histos['Mbbll_ptll'].Fill(Mbbll)
                                                            if (Mbb > Mbb_min) and (Mbb < Mbb_max):
                                                                self.histos['Mbb_mbb'].Fill(Mbb)
                                                                self.histos['Mbbll_mbb'].Fill(Mbbll)
                                                                if (Mbbll > Mbbll_min) and (Mbbll < Mbbll_max):
                                                                    self.histos['Mbb_mbbll'].Fill(Mbb)
                                                                    self.histos['Mbbll_mbbll'].Fill(Mbbll)
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
         